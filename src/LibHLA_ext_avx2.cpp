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
// Name           : LibHLA
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


#include "LibHLA.h"
#include <R.h>

using namespace std;
using namespace HLA_LIB;


#if defined(HIBAG_CPU_ARCH_X86) && defined(__GNUC__) && ((__GNUC__>4) || (__GNUC__==4 && __GNUC_MINOR__>=8))
#   define HIBAG_CPU_ARCH_X86_AVX2
#endif

#if defined(HIBAG_CPU_ARCH_X86) && !defined(HIBAG_CPU_ARCH_X86_AVX2) && defined(__clang__)
#   define HIBAG_CPU_ARCH_X86_AVX2
#endif

#if defined(__AVX2__) && !defined(HIBAG_CPU_ARCH_X86_AVX2)
#   define HIBAG_CPU_ARCH_X86_AVX2
#endif

#ifdef HIBAG_CPU_ARCH_X86_AVX2
extern const bool HIBAG_ALGORITHM_AVX2 = true;
#else
extern const bool HIBAG_ALGORITHM_AVX2 = false;
#endif


#define SIMD_NAME(NAME)  NAME ## _avx2
#define THROW_ERROR      throw ErrHLA("No AVX2 support!")


#ifdef HIBAG_CPU_ARCH_X86_AVX2

#   include <xmmintrin.h>  // SSE
#   include <emmintrin.h>  // SSE2
#   include <immintrin.h>  // AVX, AVX2
#   if !defined(__AVX2__) && !defined(__clang__)
		#pragma GCC target("avx2")
#   endif

#define TARGET_AVX2    __attribute__((target("avx2")))
#undef SIMD_NAME
#define SIMD_NAME(NAME)  TARGET_AVX2 NAME ## _avx2

#if defined(__MINGW32__) || defined(__MINGW64__)
#   define SIMD_ANDNOT_I256(x1, x2)    _mm256_andnot_si256(x1, x2)
#else
#   define SIMD_ANDNOT_I256(x1, x2)    (x2) & ~(x1)
#endif


typedef int64_t UTYPE;
#define U_POPCOUNT    __builtin_popcountll
#define U_H0(x, i)    ((UTYPE*)&x[i].PackedHaplo[0])[0]
#define U_H1(x, i)    ((UTYPE*)&x[i].PackedHaplo[0])[1]


/// Prepare the internal genotype structure
struct TGenoStruct
{
public:
	__m128i S1, S2;  ///< packed genotypes
	__m256i S1_0, S2_0, S1_1, S2_1;
	bool Low64b;     ///< whether length <= 64 or not
	/// constructor
	TGenoStruct(size_t Length, const TGenotype &G)
	{
		const UTYPE *s1 = (const UTYPE*)&G.PackedSNP1[0];
		const UTYPE *s2 = (const UTYPE*)&G.PackedSNP2[0];
		Low64b = (Length <= 64);
		if (Low64b)
		{
			__m128i I1 = { s1[0], s2[0] }, I2 = { s2[0], s1[0] };  // genotypes
			S1 = I1; S2 = I2;
			__m256i J1 = { s1[0], s1[0], s1[0], s1[0] };
			__m256i J2 = { s2[0], s2[0], s2[0], s2[0] };
			S1_0 = J1; S2_0 = J2;
		} else {
			__m128i I1 = { s1[0], s1[1] }, I2 = { s2[0], s2[1] };  // genotypes
			S1 = I1; S2 = I2;
			__m256i J1_0 = { s1[0], s1[0], s1[0], s1[0] };
			__m256i J2_0 = { s2[0], s2[0], s2[0], s2[0] };
			S1_0 = J1_0; S2_0 = J2_0;
			__m256i J1_1 = { s1[1], s1[1], s1[1], s1[1] };
			__m256i J2_1 = { s2[1], s2[1], s2[1], s2[1] };
			S1_1 = J1_1; S2_1 = J2_1;
		}
	}
};


static ALWAYS_INLINE TARGET_AVX2
	int hamm_d(const TGenoStruct &G, const THaplotype &H1, const THaplotype &H2)
{
	const UTYPE *h1 = (const UTYPE*)&H1.PackedHaplo[0];
	const UTYPE *h2 = (const UTYPE*)&H2.PackedHaplo[0];
	// here, UTYPE = int64_t
	if (G.Low64b)
	{
		__m128i H  = { *h1, *h2 };  // two haplotypes
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


// defined for Wojciech Mula algorithm's popcnt in 64-bit integers
static const __m256i pcnt_lookup = {
		0x0302020102010100LL, 0x0403030203020201LL,
		0x0302020102010100LL, 0x0403030203020201LL };
static const __m256i pcnt_low_mask = {
		0x0F0F0F0F0F0F0F0FLL, 0x0F0F0F0F0F0F0F0FLL,
		0x0F0F0F0F0F0F0F0FLL, 0x0F0F0F0F0F0F0F0FLL };

static ALWAYS_INLINE TARGET_AVX2
	size_t add_geno_freq4(size_t n, const THaplotype *i1, const THaplotype *i2,
		const TGenoStruct &GS, double &prob)
{
	const double ff = 2 * i1->Freq;
	const __m256d ff4 = { ff, ff, ff, ff };
	if (GS.Low64b)
	{
		const __m256i H1 = _mm256_set1_epi64x(U_H0(i1, 0));
		for (; n >= 4; n -= 4, i2 += 4)
		{
			__m256i H2 = { U_H0(i2,0), U_H0(i2,1), U_H0(i2,2), U_H0(i2,3) };
			__m256i S1 = GS.S1_0, S2 = GS.S2_0;
			__m256i M = SIMD_ANDNOT_I256(S1, S2);  // missing value, 1 is missing
			__m256i MASK = SIMD_ANDNOT_I256(M, (H1 ^ S2) | (H2 ^ S1));
			__m256i va = (H1 ^ S1) & MASK, vb = (H2 ^ S2) & MASK;
			// popcount for 64b integers
			__m256i lo_a = va & pcnt_low_mask;
			__m256i lo_b = vb & pcnt_low_mask;
			__m256i hi_a = _mm256_srli_epi32(va, 4) & pcnt_low_mask;
			__m256i hi_b = _mm256_srli_epi32(vb, 4) & pcnt_low_mask;
			__m256i pc1_a = _mm256_shuffle_epi8(pcnt_lookup, lo_a);
			__m256i pc1_b = _mm256_shuffle_epi8(pcnt_lookup, lo_b);
			__m256i pc2_a = _mm256_shuffle_epi8(pcnt_lookup, hi_a);
			__m256i pc2_b = _mm256_shuffle_epi8(pcnt_lookup, hi_b);
			__m256i tot_a = _mm256_add_epi8(pc1_a, pc2_a);
			__m256i tot_b = _mm256_add_epi8(pc1_b, pc2_b);
			__m256i total = _mm256_add_epi8(tot_a, tot_b);
			__m256i ii4 = _mm256_sad_epu8(total, _mm256_setzero_si256());
			// four frequencies
			__m256d f = _mm256_i64gather_pd(EXP_LOG_MIN_RARE_FREQ, ii4, 8);
			__m256d f2 = { i2[0].Freq, i2[1].Freq, i2[2].Freq, i2[3].Freq };
			f = ff4 * f2 * f;
			prob += f[0] + f[1] + f[2] + f[3];
		}
	} else {
		const __m256i H1_0 = _mm256_set1_epi64x(U_H0(i1, 0));
		const __m256i H1_1 = _mm256_set1_epi64x(U_H1(i1, 0));
		for (; n >= 4; n -= 4, i2 += 4)
		{
			__m256i H2_0 = { U_H0(i2,0), U_H0(i2,1), U_H0(i2,2), U_H0(i2,3) };
			__m256i H2_1 = { U_H1(i2,0), U_H1(i2,1), U_H1(i2,2), U_H1(i2,3) };
			__m256i S1_0 = GS.S1_0, S2_0 = GS.S2_0;
			__m256i S1_1 = GS.S1_1, S2_1 = GS.S2_1;
			__m256i M_0 = SIMD_ANDNOT_I256(S1_0, S2_0);  // missing value, 1 is missing
			__m256i M_1 = SIMD_ANDNOT_I256(S1_1, S2_1);  // missing value, 1 is missing
			__m256i MASK_0 = SIMD_ANDNOT_I256(M_0, (H1_0 ^ S2_0) | (H2_0 ^ S1_0));
			__m256i MASK_1 = SIMD_ANDNOT_I256(M_1, (H1_1 ^ S2_1) | (H2_1 ^ S1_1));
			__m256i va_0 = (H1_0 ^ S1_0) & MASK_0, vb_0 = (H2_0 ^ S2_0) & MASK_0;
			__m256i va_1 = (H1_1 ^ S1_1) & MASK_1, vb_1 = (H2_1 ^ S2_1) & MASK_1;
			// popcount for 64b integers
			__m256i lo_a_0 = va_0 & pcnt_low_mask;
			__m256i lo_a_1 = va_1 & pcnt_low_mask;
			__m256i lo_b_0 = vb_0 & pcnt_low_mask;
			__m256i lo_b_1 = vb_1 & pcnt_low_mask;
			__m256i hi_a_0 = _mm256_srli_epi32(va_0, 4) & pcnt_low_mask;
			__m256i hi_a_1 = _mm256_srli_epi32(va_1, 4) & pcnt_low_mask;
			__m256i hi_b_0 = _mm256_srli_epi32(vb_0, 4) & pcnt_low_mask;
			__m256i hi_b_1 = _mm256_srli_epi32(vb_1, 4) & pcnt_low_mask;
			__m256i pc1_a_0 = _mm256_shuffle_epi8(pcnt_lookup, lo_a_0);
			__m256i pc1_a_1 = _mm256_shuffle_epi8(pcnt_lookup, lo_a_1);
			__m256i pc1_b_0 = _mm256_shuffle_epi8(pcnt_lookup, lo_b_0);
			__m256i pc1_b_1 = _mm256_shuffle_epi8(pcnt_lookup, lo_b_1);
			__m256i pc2_a_0 = _mm256_shuffle_epi8(pcnt_lookup, hi_a_0);
			__m256i pc2_a_1 = _mm256_shuffle_epi8(pcnt_lookup, hi_a_1);
			__m256i pc2_b_0 = _mm256_shuffle_epi8(pcnt_lookup, hi_b_0);
			__m256i pc2_b_1 = _mm256_shuffle_epi8(pcnt_lookup, hi_b_1);
			__m256i tot_a_0 = _mm256_add_epi8(pc1_a_0, pc2_a_0);
			__m256i tot_a_1 = _mm256_add_epi8(pc1_a_1, pc2_a_1);
			__m256i tot_b_0 = _mm256_add_epi8(pc1_b_0, pc2_b_0);
			__m256i tot_b_1 = _mm256_add_epi8(pc1_b_1, pc2_b_1);
			__m256i tot_0 = _mm256_add_epi8(tot_a_0, tot_b_0);
			__m256i tot_1 = _mm256_add_epi8(tot_a_1, tot_b_1);
			__m256i total = _mm256_add_epi8(tot_0, tot_1);
			__m256i ii4 = _mm256_sad_epu8(total, _mm256_setzero_si256());
			// four frequencies
			__m256d f = _mm256_i64gather_pd(EXP_LOG_MIN_RARE_FREQ, ii4, 8);
			__m256d f2 = { i2[0].Freq, i2[1].Freq, i2[2].Freq, i2[3].Freq };
			f = ff4 * f2 * f;
			prob += f[0] + f[1] + f[2] + f[3];
		}
	}
	return n;
}


THLAType SIMD_NAME(CAlg_Prediction::_BestGuess)(const CHaplotypeList &Haplo,
	const TGenotype &Geno)
{
	const TGenoStruct GS(Haplo.Num_SNP, Geno);
	THLAType rv;
	rv.Allele1 = rv.Allele2 = NA_INTEGER;
	double max=0, prob;
	THaplotype *I1=Haplo.List, *I2;

	for (int h1=0; h1 < _nHLA; h1++)
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
				m2 = add_geno_freq4(m2, i1, i2, GS, prob);
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
		for (int h2=h1+1; h2 < _nHLA; h2++)
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
					m2 = add_geno_freq4(m2, i1, i2, GS, prob);
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
	const TGenoStruct GS(Haplo.Num_SNP, Geno);
	int H1=HLA.Allele1, H2=HLA.Allele2;
	if (H1 > H2) std::swap(H1, H2);
	int IxHLA = H2 + H1*(2*_nHLA-H1-1)/2;
	int idx = 0;
	double sum=0, hlaProb=0, prob;
	THaplotype *I1=Haplo.List, *I2;

	for (int h1=0; h1 < _nHLA; h1++)
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
				m2 = add_geno_freq4(m2, i1, i2, GS, prob);
				i2 += m1 - 1 - m2;
			}
			for (; m2 > 0; m2--, i2++)
				ADD_FREQ_MUTANT(prob, ff * i2->Freq, hamm_d(GS, *i1, *i2));
		}
		I2 = I1 + n1;
		if (IxHLA == idx) hlaProb = prob;
		idx ++; sum += prob;

		// off-diagonal
		for (int h2=h1+1; h2 < _nHLA; h2++)
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
					m2 = add_geno_freq4(m2, i1, i2, GS, prob);
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


void SIMD_NAME(CAlg_Prediction::_PostProb2)(const CHaplotypeList &Haplo,
	const TGenotype &Geno, double &SumProb)
{
	const TGenoStruct GS(Haplo.Num_SNP, Geno);
	double *pProb = &_PostProb[0], sum;
	THaplotype *I1=Haplo.List, *I2;

	for (int h1=0; h1 < _nHLA; h1++)
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
			const double ff = 2 * i1->Freq;
			THaplotype *i2 = i1 + 1;
			size_t m2 = m1 - 1;
			if (m2 >= 4)
			{
				m2 = add_geno_freq4(m2, i1, i2, GS, sum);
				i2 += m1 - 1 - m2;
			}
			for (; m2 > 0; m2--, i2++)
				ADD_FREQ_MUTANT(sum, ff * i2->Freq, hamm_d(GS, *i1, *i2));
		}
		*pProb++ = sum;
		I2 = I1 + n1;

		// off-diagonal
		for (int h2=h1+1; h2 < _nHLA; h2++)
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
					m2 = add_geno_freq4(m2, i1, i2, GS, sum);
					i2 += n2 - m2;
				}
				for (; m2 > 0; m2--, i2++)
					ADD_FREQ_MUTANT(sum, ff * i2->Freq, hamm_d(GS, *i1, *i2));
			}
			*pProb++ = sum;
			I2 += n2;
		}

		I1 += n1;
	}

	// normalize
	sum = 0;
	double *p = &_PostProb[0];
	for (size_t n = _PostProb.size(); n > 0; n--) sum += *p++;
	SumProb = sum;
	sum = 1 / sum;
	p = &_PostProb[0];
	for (size_t n = _PostProb.size(); n > 0; n--) *p++ *= sum;
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

void SIMD_NAME(CAlg_Prediction::_PostProb2)(const CHaplotypeList &Haplo,
	const TGenotype &Geno, double &SumProb)
{
	THROW_ERROR;
}

#endif
