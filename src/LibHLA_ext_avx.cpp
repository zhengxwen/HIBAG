// ===============================================================
//
// HIBAG R package (HLA Genotype Imputation with Attribute Bagging)
// Copyright (C) 2020   Xiuwen Zheng (zhengx@u.washington.edu)
// All rights reserved.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// ===============================================================
// Name           : LibHLA_ext_avx
// Author         : Xiuwen Zheng
// Kernel Version : 1.5
// Copyright      : Xiuwen Zheng (GPL v3)
// Description    : HLA imputation C++ library
// ===============================================================


// Optimization level
#ifndef HIBAG_NO_COMPILER_OPTIMIZE
#if defined(__clang__) && !defined(__APPLE__)
    #pragma clang optimize on
#endif
#if defined(__GNUC__) && ((__GNUC__>4) || (__GNUC__==4 && __GNUC_MINOR__>=4))
    #pragma GCC optimize("O3")
#endif
#endif


#include "LibHLA_ext.h"
// need a patch for gcc_v4.8
#if defined(HIBAG_CPU_ARCH_X86) && defined(__GNUC__) && (__GNUC__==4) && (__GNUC_MINOR__==8)
#   pragma GCC target("avx")
#   define __AVX__       1
#   define __SSE4_1__    1
#   define __SSE4_2__    1
#   define __SSE3__      1
#   define __SSSE3__     1
#   define __POPCNT__    1
#   include <xmmintrin.h>  // SSE
#   include <emmintrin.h>  // SSE2
#   include <immintrin.h>  // AVX
#   define TARGET_AVX
#endif

#include "LibHLA.h"
#include <R.h>

using namespace std;
using namespace HLA_LIB;


#ifdef HIBAG_CPU_ARCH_X86_AVX
extern const bool HIBAG_ALGORITHM_AVX = true;
#else
extern const bool HIBAG_ALGORITHM_AVX = false;
#endif


#define SIMD_NAME(NAME)  NAME ## _avx
#define THROW_ERROR      throw ErrHLA("No AVX support!")


#ifdef HIBAG_CPU_ARCH_X86_AVX

#ifdef __ICC
	#pragma intel optimization_parameter target_arch=AVX
#elif !defined(__AVX__) && !defined(__clang__)
	#pragma GCC target("avx")
#endif

#include <xmmintrin.h>  // SSE
#include <emmintrin.h>  // SSE2
#include <immintrin.h>  // AVX, AVX2

#ifndef TARGET_AVX
#   define TARGET_AVX    __attribute__((target("avx")))
#endif
#undef SIMD_NAME
#define SIMD_NAME(NAME)  TARGET_AVX NAME ## _avx

#if defined(__ICC)
#   define SIMD_ANDNOT_I256(x1, x2)    \
		(__m256i)_mm256_andnot_pd((__m256d)(x1), (__m256d)(x2))
#else
#   define SIMD_ANDNOT_I256(x1, x2)    (x2) & ~(x1)
#endif


#define U_POPCOUNT    __builtin_popcountll


/// Prepare the internal genotype structure
struct TGenoStruct_avx
{
public:
	__m128i S1, S2;  ///< packed genotypes
	__m256i S1_0, S2_0, S1_1, S2_1;
	bool Low64b;     ///< whether length <= 64 or not
	int64_t *p_H_0, *p_H_1;
	double *p_Freq;
	/// constructor
#ifndef __ICC
	TARGET_AVX
#endif
	inline TGenoStruct_avx(const CHaplotypeList &Haplo, const TGenotype &G)
	{
		const INT64 *s1 = G.PackedSNP1, *s2 = G.PackedSNP2;
		Low64b = (Haplo.Num_SNP <= 64);
		if (Low64b)
		{
			__m128i I1 = { s1[0], s2[0] }, I2 = { s2[0], s1[0] };  // genotypes
			S1 = I1; S2 = I2;
			S1_0 = _mm256_set1_epi64x(s1[0]);
			S2_0 = _mm256_set1_epi64x(s2[0]);
		} else {
			__m128i I1 = { s1[0], s1[1] }, I2 = { s2[0], s2[1] };  // genotypes
			S1 = I1; S2 = I2;
			S1_0 = _mm256_set1_epi64x(s1[0]);
			S2_0 = _mm256_set1_epi64x(s2[0]);
			S1_1 = _mm256_set1_epi64x(s1[1]);
			S2_1 = _mm256_set1_epi64x(s2[1]);
		}
		p_H_0 = Haplo.aux_haplo;
		p_H_1 = Haplo.aux_haplo + Haplo.Num_Haplo;
		p_Freq = Haplo.aux_freq;
	}
};


static ALWAYS_INLINE TARGET_AVX
	int hamm_d(const TGenoStruct_avx &G, const THaplotype &H1, const THaplotype &H2)
{
	const INT64 *h1 = H1.PackedHaplo, *h2 = H2.PackedHaplo;
	// here, UTYPE = int64_t
	if (G.Low64b)
	{
		__m128i H  = { h1[0], h2[0] };  // two haplotypes
		__m128i S1 = G.S1, S2 = G.S2;  // genotypes
		__m128i m1 = H ^ S2, m2 = { m1[1], m1[0] };
		// worry about n < UTYPE_BIT_NUM? unused bits are set to be a missing flag
		__m128i M = _mm_andnot_si128(S1, S2);  // missing value, 1 is missing
		__m128i M2 = { M[0], M[0] };
		__m128i MASK = _mm_andnot_si128(M2, m1 | m2);
		__m128i v = (H ^ S1) & MASK;  // (H1 ^ S1) & MASK, (H2 ^ S2) & MASK
		// popcount
		return U_POPCOUNT(v[0]) + U_POPCOUNT(v[1]);
	} else {
		// here, __m128i is faster than using __m256i
		// since HIBAG_MAXNUM_SNP_IN_CLASSIFIER = 128
		__m128i H1 = { h1[0], h1[1] }, H2 = { h2[0], h2[1] };  // two haplotypes
		__m128i S1 = G.S1, S2 = G.S2;  // genotypes
		// worry about n < UTYPE_BIT_NUM? unused bits are set to be a missing flag
		__m128i M = _mm_andnot_si128(S1, S2);  // missing value, 1 is missing
		__m128i MASK = _mm_andnot_si128(M, (H1 ^ S2) | (H2 ^ S1));
		__m128i va = (H1 ^ S1) & MASK;  // (H1 ^ S1) & MASK, (H2 ^ S2) & MASK
		__m128i vb = (H2 ^ S2) & MASK;
		// popcount
		return U_POPCOUNT(va[0]) + U_POPCOUNT(va[1]) +
			U_POPCOUNT(vb[0]) + U_POPCOUNT(vb[1]);
	}
}


static inline TARGET_AVX
	size_t add_geno_freq4(size_t n, const THaplotype *i1, size_t i2,
		const TGenoStruct_avx &GS, double &prob)
{
	const double ff = 2 * i1->Freq;
	if (GS.Low64b)
	{
		const __m256i H1 = _mm256_set1_epi64x(i1[0].PackedHaplo[0]);
		const __m256i S1 = GS.S1_0, S2 = GS.S2_0;
		const __m256i M = SIMD_ANDNOT_I256(S1, S2);  // missing value, 1 is missing
		for (; n >= 4; n -= 4, i2 += 4)
		{
			__m256i H2 = _mm256_loadu_si256((__m256i*)(GS.p_H_0 + i2));
			__m256i MASK = SIMD_ANDNOT_I256(M, (H1 ^ S2) | (H2 ^ S1));
			__m256i va = (H1 ^ S1) & MASK, vb = (H2 ^ S2) & MASK;
			// popcount for 64b integers
			__m256d f2 = _mm256_set1_pd(ff) * _mm256_loadu_pd(GS.p_Freq + i2);
			prob += f2[0] *
				EXP_LOG_MIN_RARE_FREQ[U_POPCOUNT(va[0]) + U_POPCOUNT(vb[0])];
			prob += f2[1] *
				EXP_LOG_MIN_RARE_FREQ[U_POPCOUNT(va[1]) + U_POPCOUNT(vb[1])];
			prob += f2[2] *
				EXP_LOG_MIN_RARE_FREQ[U_POPCOUNT(va[2]) + U_POPCOUNT(vb[2])];
			prob += f2[3] *
				EXP_LOG_MIN_RARE_FREQ[U_POPCOUNT(va[3]) + U_POPCOUNT(vb[3])];
		}
	} else {
		const __m256i H1_0 = _mm256_set1_epi64x(i1[0].PackedHaplo[0]);
		const __m256i H1_1 = _mm256_set1_epi64x(i1[0].PackedHaplo[1]);
		const __m256i S1_0 = GS.S1_0, S2_0 = GS.S2_0;
		const __m256i S1_1 = GS.S1_1, S2_1 = GS.S2_1;
		const __m256i M_0 = SIMD_ANDNOT_I256(S1_0, S2_0);  // missing value, 1 is missing
		const __m256i M_1 = SIMD_ANDNOT_I256(S1_1, S2_1);  // missing value, 1 is missing
		for (; n >= 4; n -= 4, i2 += 4)
		{
			__m256i H2_0 = _mm256_loadu_si256((__m256i*)(GS.p_H_0 + i2));
			__m256i H2_1 = _mm256_loadu_si256((__m256i*)(GS.p_H_1 + i2));
			__m256i MASK_0 = SIMD_ANDNOT_I256(M_0, (H1_0 ^ S2_0) | (H2_0 ^ S1_0));
			__m256i MASK_1 = SIMD_ANDNOT_I256(M_1, (H1_1 ^ S2_1) | (H2_1 ^ S1_1));
			__m256i va_0 = (H1_0 ^ S1_0) & MASK_0, vb_0 = (H2_0 ^ S2_0) & MASK_0;
			__m256i va_1 = (H1_1 ^ S1_1) & MASK_1, vb_1 = (H2_1 ^ S2_1) & MASK_1;
			// popcount for 64b integers
			__m256d f2 = _mm256_set1_pd(ff) * _mm256_loadu_pd(GS.p_Freq + i2);
			prob += f2[0] * EXP_LOG_MIN_RARE_FREQ[ U_POPCOUNT(va_0[0]) +
				U_POPCOUNT(vb_0[0]) + U_POPCOUNT(va_1[0]) + U_POPCOUNT(vb_1[0]) ];
			prob += f2[1] * EXP_LOG_MIN_RARE_FREQ[ U_POPCOUNT(va_0[1]) +
				U_POPCOUNT(vb_0[1]) + U_POPCOUNT(va_1[1]) + U_POPCOUNT(vb_1[1]) ];
			prob += f2[2] * EXP_LOG_MIN_RARE_FREQ[ U_POPCOUNT(va_0[2]) +
				U_POPCOUNT(vb_0[2]) + U_POPCOUNT(va_1[2]) + U_POPCOUNT(vb_1[2]) ];
			prob += f2[3] * EXP_LOG_MIN_RARE_FREQ[ U_POPCOUNT(va_0[3]) +
				U_POPCOUNT(vb_0[3]) + U_POPCOUNT(va_1[3]) + U_POPCOUNT(vb_1[3]) ];
		}
	}
	return n;
}


THLAType SIMD_NAME(CAlg_Prediction::_BestGuess)(const CHaplotypeList &Haplo,
	const TGenotype &Geno)
{
	const TGenoStruct_avx GS(Haplo, Geno);
	THLAType rv;
	rv.Allele1 = rv.Allele2 = NA_INTEGER;
	double max=0, prob;
	const int nHLA = Haplo.nHLA();
	THaplotype *base=Haplo.List, *I1=base, *I2;

	for (int h1=0; h1 < nHLA; h1++)
	{
		const size_t n1 = Haplo.LenPerHLA[h1];

		// diagonal
		prob = 0;
		THaplotype *i1 = I1;
		for (size_t m1=n1; m1 > 0; m1--, i1++)
		{
			// i2 = i1
			ADD_FREQ_MUTANT(prob, i1->Freq * i1->Freq, hamm_d(GS, *i1, *i1));
			// i2 > i1
			const double ff = 2 * i1->Freq;
			THaplotype *i2 = i1 + 1;
			size_t m2 = m1 - 1;
			if (m2 >= 4)
			{
				m2 = add_geno_freq4(m2, i1, i2-base, GS, prob);
				i2 += m1 - 1 - m2;
			}
			for (; m2 > 0; m2--, i2++)
				ADD_FREQ_MUTANT(prob, ff * i2->Freq, hamm_d(GS, *i1, *i2));
		}
		I2 = I1 + n1;
		if (max < prob)
		{
			max = prob;
			rv.Allele1 = rv.Allele2 = h1;
		}

		// off-diagonal
		for (int h2=h1+1; h2 < nHLA; h2++)
		{
			const size_t n2 = Haplo.LenPerHLA[h2];
			prob = 0;
			THaplotype *i1 = I1;
			for (size_t m1=n1; m1 > 0; m1--, i1++)
			{
				const double ff = 2 * i1->Freq;
				THaplotype *i2 = I2;
				size_t m2 = n2;
				if (m2 >= 4)
				{
					m2 = add_geno_freq4(m2, i1, i2-base, GS, prob);
					i2 += n2 - m2;
				}
				for (; m2 > 0; m2--, i2++)
					ADD_FREQ_MUTANT(prob, ff * i2->Freq, hamm_d(GS, *i1, *i2));
			}
			I2 += n2;
			if (max < prob)
			{
				max = prob;
				rv.Allele1 = h1; rv.Allele2 = h2;
			}
		}

		I1 += n1;
	}

	return rv;
}


double SIMD_NAME(CAlg_Prediction::_PostProb)(const CHaplotypeList &Haplo,
	const TGenotype &Geno, const THLAType &HLA)
{
	const TGenoStruct_avx GS(Haplo, Geno);
	int H1=HLA.Allele1, H2=HLA.Allele2;
	if (H1 > H2) std::swap(H1, H2);
	const int nHLA = Haplo.nHLA();
	int IxHLA = H2 + H1*(2*nHLA-H1-1)/2;
	int idx = 0;
	double sum=0, hlaProb=0, prob;
	THaplotype *base=Haplo.List, *I1=base, *I2;

	for (int h1=0; h1 < nHLA; h1++)
	{
		const size_t n1 = Haplo.LenPerHLA[h1];

		// diagonal
		prob = 0;
		THaplotype *i1 = I1;
		for (size_t m1=n1; m1 > 0; m1--, i1++)
		{
			// i2 = i1
			ADD_FREQ_MUTANT(prob, i1->Freq * i1->Freq, hamm_d(GS, *i1, *i1));
			// i2 > i1
			const double ff = 2 * i1->Freq;
			THaplotype *i2 = i1 + 1;
			size_t m2 = m1 - 1;
			if (m2 >= 4)
			{
				m2 = add_geno_freq4(m2, i1, i2-base, GS, prob);
				i2 += m1 - 1 - m2;
			}
			for (; m2 > 0; m2--, i2++)
				ADD_FREQ_MUTANT(prob, ff * i2->Freq, hamm_d(GS, *i1, *i2));
		}
		I2 = I1 + n1;
		if (IxHLA == idx) hlaProb = prob;
		idx ++; sum += prob;

		// off-diagonal
		for (int h2=h1+1; h2 < nHLA; h2++)
		{
			const size_t n2 = Haplo.LenPerHLA[h2];
			prob = 0;
			THaplotype *i1 = I1;
			for (size_t m1=n1; m1 > 0; m1--, i1++)
			{
				const double ff = 2 * i1->Freq;
				THaplotype *i2 = I2;
				size_t m2 = n2;
				if (m2 >= 4)
				{
					m2 = add_geno_freq4(m2, i1, i2-base, GS, prob);
					i2 += n2 - m2;
				}
				for (; m2 > 0; m2--, i2++)
					ADD_FREQ_MUTANT(prob, ff * i2->Freq, hamm_d(GS, *i1, *i2));
			}
			I2 += n2;
			if (IxHLA == idx) hlaProb = prob;
			idx ++; sum += prob;
		}

		I1 += n1;
	}

	return hlaProb / sum;
}


double SIMD_NAME(CAlg_Prediction::_PostProb2)(const CHaplotypeList &Haplo,
	const TGenotype &Geno, double Prob[])
{
	const TGenoStruct_avx GS(Haplo, Geno);
	double *p = Prob, sum;
	const int nHLA = Haplo.nHLA();
	THaplotype *base=Haplo.List, *I1=base, *I2;

	for (int h1=0; h1 < nHLA; h1++)
	{
		const size_t n1 = Haplo.LenPerHLA[h1];

		// diagonal
		sum = 0;
		THaplotype *i1 = I1;
		for (size_t m1=n1; m1 > 0; m1--, i1++)
		{
			// i2 = i1
			ADD_FREQ_MUTANT(sum, i1->Freq * i1->Freq, hamm_d(GS, *i1, *i1));
			// i2 > i1
			THaplotype *i2 = i1 + 1;
			size_t m2 = m1 - 1;
			if (m2 >= 4)
			{
				m2 = add_geno_freq4(m2, i1, i2-base, GS, sum);
				i2 += m1 - 1 - m2;
			}
			const double ff = 2 * i1->Freq;
			for (; m2 > 0; m2--, i2++)
				ADD_FREQ_MUTANT(sum, ff * i2->Freq, hamm_d(GS, *i1, *i2));
		}
		*p++ = sum;
		I2 = I1 + n1;

		// off-diagonal
		for (int h2=h1+1; h2 < nHLA; h2++)
		{
			const size_t n2 = Haplo.LenPerHLA[h2];
			sum = 0;
			THaplotype *i1 = I1;
			for (size_t m1=n1; m1 > 0; m1--, i1++)
			{
				const double ff = 2 * i1->Freq;
				THaplotype *i2 = I2;
				size_t m2 = n2;
				if (m2 >= 4)
				{
					m2 = add_geno_freq4(m2, i1, i2-base, GS, sum);
					i2 += n2 - m2;
				}
				for (; m2 > 0; m2--, i2++)
					ADD_FREQ_MUTANT(sum, ff * i2->Freq, hamm_d(GS, *i1, *i2));
			}
			*p++ = sum;
			I2 += n2;
		}

		I1 += n1;
	}

	// normalize
	const size_t n = nHLA*(nHLA+1)/2;
	sum = 0;
	for (size_t i=0; i < n; i++) sum += Prob[i];
	const double ff = 1 / sum;
	for (size_t i=0; i < n; i++) Prob[i] *= ff;
	return sum;
}


#else

THLAType SIMD_NAME(CAlg_Prediction::_BestGuess)(const CHaplotypeList &Haplo,
	const TGenotype &Geno)
{
	THROW_ERROR;
}

double SIMD_NAME(CAlg_Prediction::_PostProb)(const CHaplotypeList &Haplo,
	const TGenotype &Geno, const THLAType &HLA)
{
	THROW_ERROR;
}

double SIMD_NAME(CAlg_Prediction::_PostProb2)(const CHaplotypeList &Haplo,
	const TGenotype &Geno, double Prob[])
{
	THROW_ERROR;
}

#endif
