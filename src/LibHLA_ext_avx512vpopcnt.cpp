// ===============================================================
//
// HIBAG R package (HLA Genotype Imputation with Attribute Bagging)
// Copyright (C) 2021   Xiuwen Zheng (zhengx@u.washington.edu)
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
// Name           : LibHLA_ext_avx512vpopcnt
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


#ifdef HIBAG_CPU_ARCH_X86_AVX512VPOPCNTDQ
extern const bool HIBAG_ALGORITHM_AVX512VPOPCNTDQ = true;
#else
extern const bool HIBAG_ALGORITHM_AVX512VPOPCNTDQ = false;
#endif


#define SIMD_NAME(NAME)  NAME ## _avx512vpopcnt
#define THROW_ERROR      throw ErrHLA("No AVX512VPOPCNT support!")


#ifdef HIBAG_CPU_ARCH_X86_AVX512VPOPCNTDQ

#ifdef __ICC
	#pragma intel optimization_parameter target_arch=ICELAKE
#   define TARGET_AVX512
#else
#   if !defined(__AVX512F__) && !defined(__clang__)
		#pragma GCC target("avx512f")
#   endif
#   if !defined(__AVX512VL__) && !defined(__clang__)
		#pragma GCC target("avx512vl")
#   endif
#   if !defined(__AVX512VPOPCNTDQ__) && !defined(__clang__)
		#pragma GCC target("avx512vpopcntdq")
#   endif
#   define TARGET_AVX512    __attribute__((target("avx512f,avx512vl,avx512vpopcntdq")))
#endif

#include <xmmintrin.h>  // SSE
#include <emmintrin.h>  // SSE2
#include <immintrin.h>  // AVX, AVX2, AVX512F, AVX512VL, AVX512BW, , AVX512VPOPCNT

#undef SIMD_NAME
#define SIMD_NAME(NAME)  TARGET_AVX512 NAME ## _avx512vpopcnt

#if defined(__MINGW32__) || defined(__MINGW64__) || defined(__ICC)
#   define SIMD_ANDNOT_I256(x1, x2)    _mm256_andnot_si256(x1, x2)
#   define SIMD_ANDNOT_I512(x1, x2)    _mm512_andnot_si512(x1, x2)
#else
#   define SIMD_ANDNOT_I256(x1, x2)    (x2) & ~(x1)
#   define SIMD_ANDNOT_I512(x1, x2)    (x2) & ~(x1)
#endif


#define U_POPCOUNT    __builtin_popcountll


/// Prepare the internal genotype structure
struct TGenoStruct_512vpopcnt
{
public:
	__m128i S1, S2;  ///< packed genotypes
	__m512i S1_0, S2_0, S1_1, S2_1;
	bool Low64b;     ///< whether length <= 64 or not
	int64_t *p_H_0, *p_H_1;
	double *p_Freq;

	/// constructor
	inline TGenoStruct_512vpopcnt(const CHaplotypeList &Haplo, const TGenotype &G)
	{
		init(Haplo.Num_SNP, G);
		p_H_0 = Haplo.aux_haplo;
		p_H_1 = Haplo.aux_haplo + Haplo.Num_Haplo;
		p_Freq = Haplo.aux_freq;
	}
	inline TGenoStruct_512vpopcnt(const size_t Num_SNP, const TGenotype &G)
	{
		init(Num_SNP, G);
		p_H_0 = p_H_1 = NULL; p_Freq = NULL;
	}

private:
#ifndef __ICC
	TARGET_AVX512
#endif
	inline void init(const size_t Num_SNP, const TGenotype &G)
	{
		const INT64 *s1 = G.PackedSNP1, *s2 = G.PackedSNP2;
		Low64b = (Num_SNP <= 64);
		if (Low64b)
		{
			__m128i I1 = { s1[0], s2[0] }, I2 = { s2[0], s1[0] };  // genotypes
			S1 = I1; S2 = I2;
			S1_0 = _mm512_set1_epi64(s1[0]);
			S2_0 = _mm512_set1_epi64(s2[0]);
		} else {
			__m128i I1 = { s1[0], s1[1] }, I2 = { s2[0], s2[1] };  // genotypes
			S1 = I1; S2 = I2;
			S1_0 = _mm512_set1_epi64(s1[0]);
			S2_0 = _mm512_set1_epi64(s2[0]);
			S1_1 = _mm512_set1_epi64(s1[1]);
			S2_1 = _mm512_set1_epi64(s2[1]);
		}
	}
};


static ALWAYS_INLINE TARGET_AVX512
	int hamm_d(const TGenoStruct_512vpopcnt &G, const THaplotype &H1, const THaplotype &H2)
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
		__m128i v = _mm_popcnt_epi64(va) + _mm_popcnt_epi64(vb);
		return v[0] + v[1];
	}
}


static inline TARGET_AVX512
	size_t add_geno_freq4(size_t n, const THaplotype *i1, size_t i2,
		const TGenoStruct_512vpopcnt &GS, double &prob)
{
	const double ff = 2 * i1->Freq;
	if (GS.Low64b)
	{
		const __m512i H1_8 = _mm512_set1_epi64(i1[0].PackedHaplo[0]);
		const __m512i S1 = GS.S1_0, S2 = GS.S2_0;
		const __m512i M = SIMD_ANDNOT_I512(S1, S2);  // missing value, 1 is missing
		__m512d sum8 = _mm512_setzero_pd();
		for (; n >= 8; n -= 8, i2 += 8)
		{
			__m512i H2_8 = _mm512_loadu_si512(GS.p_H_0 + i2);
			__m512i MASK = SIMD_ANDNOT_I512(M, (H1_8 ^ S2) | (H2_8 ^ S1));
			__m512i va = (H1_8 ^ S1) & MASK, vb = (H2_8 ^ S2) & MASK;
			// popcount for 64b integers
			__m512i ii4 = _mm512_popcnt_epi64(va) + _mm512_popcnt_epi64(vb);
			// eight frequencies
			__m512d f = _mm512_i64gather_pd(ii4, EXP_LOG_MIN_RARE_FREQ, 8);
			__m512d f2 = _mm512_loadu_pd(GS.p_Freq + i2);
			// maybe different behavior due to rounding error of addition
			sum8 = _mm512_fmadd_pd(_mm512_set1_pd(ff)*f2, f, sum8);
		}
		__m256d sum4 = _mm512_castpd512_pd256(sum8) + _mm512_extractf64x4_pd(sum8,1);
		if (n >= 4)
		{
			__m256i H1 = _mm512_castsi512_si256(H1_8);
			__m256i H2 = _mm256_loadu_si256((__m256i*)(GS.p_H_0 + i2));
			__m256i S1 = _mm512_castsi512_si256(GS.S1_0);
			__m256i S2 = _mm512_castsi512_si256(GS.S2_0);
			__m256i M = SIMD_ANDNOT_I256(S1, S2);  // missing value, 1 is missing
			__m256i MASK = SIMD_ANDNOT_I256(M, (H1 ^ S2) | (H2 ^ S1));
			__m256i va = (H1 ^ S1) & MASK, vb = (H2 ^ S2) & MASK;
			// popcount for 64b integers
			__m256i ii4 = _mm256_popcnt_epi64(va) + _mm256_popcnt_epi64(vb);
			// four frequencies
			__m256d f = _mm256_i64gather_pd(EXP_LOG_MIN_RARE_FREQ, ii4, 8);
			__m256d f2 = _mm256_loadu_pd(GS.p_Freq + i2);
			sum4 += _mm256_set1_pd(ff) * f2 * f;
			n -= 4;
		}
		// maybe different behavior due to rounding error of addition
		__m128d a = _mm256_castpd256_pd128(sum4) + _mm256_extractf128_pd(sum4,1);
		prob += a[0] + a[1];
	} else {
		const __m512i H1_0_8 = _mm512_set1_epi64(i1[0].PackedHaplo[0]);
		const __m512i H1_1_8 = _mm512_set1_epi64(i1[0].PackedHaplo[1]);
		const __m512i S1_0 = GS.S1_0, S2_0 = GS.S2_0;
		const __m512i S1_1 = GS.S1_1, S2_1 = GS.S2_1;
		const __m512i M_0 = SIMD_ANDNOT_I512(S1_0, S2_0);  // missing value, 1 is missing
		const __m512i M_1 = SIMD_ANDNOT_I512(S1_1, S2_1);  // missing value, 1 is missing
		__m512d sum8 = _mm512_setzero_pd();
		for (; n >= 8; n -= 8, i2 += 8)
		{
			__m512i H2_0_8 = _mm512_loadu_si512((__m512i*)(GS.p_H_0 + i2));
			__m512i H2_1_8 = _mm512_loadu_si512((__m512i*)(GS.p_H_1 + i2));
			__m512i MASK_0 = SIMD_ANDNOT_I512(M_0, (H1_0_8 ^ S2_0) | (H2_0_8 ^ S1_0));
			__m512i MASK_1 = SIMD_ANDNOT_I512(M_1, (H1_1_8 ^ S2_1) | (H2_1_8 ^ S1_1));
			__m512i va_0 = (H1_0_8 ^ S1_0) & MASK_0, vb_0 = (H2_0_8 ^ S2_0) & MASK_0;
			__m512i va_1 = (H1_1_8 ^ S1_1) & MASK_1, vb_1 = (H2_1_8 ^ S2_1) & MASK_1;
			// popcount for 64b integers
			__m512i ii4 = _mm512_popcnt_epi64(va_0) + _mm512_popcnt_epi64(vb_0) +
				_mm512_popcnt_epi64(va_1) + _mm512_popcnt_epi64(vb_1);
			// eight frequencies
			__m512d f = _mm512_i64gather_pd(ii4, EXP_LOG_MIN_RARE_FREQ, 8);
			__m512d f2 = _mm512_loadu_pd(GS.p_Freq + i2);
			// maybe different behavior due to rounding error of addition
			sum8 = _mm512_fmadd_pd(_mm512_set1_pd(ff)*f2, f, sum8);
		}
		__m256d sum4 = _mm512_castpd512_pd256(sum8) + _mm512_extractf64x4_pd(sum8,1);
		if (n >= 4)
		{
			__m256i H1_0 = _mm512_castsi512_si256(H1_0_8);
			__m256i H1_1 = _mm512_castsi512_si256(H1_1_8);
			__m256i H2_0 = _mm256_loadu_si256((__m256i*)(GS.p_H_0 + i2));
			__m256i H2_1 = _mm256_loadu_si256((__m256i*)(GS.p_H_1 + i2));
			__m256i S1_0 = _mm512_castsi512_si256(GS.S1_0);
			__m256i S2_0 = _mm512_castsi512_si256(GS.S2_0);
			__m256i S1_1 = _mm512_castsi512_si256(GS.S1_1);
			__m256i S2_1 = _mm512_castsi512_si256(GS.S2_1);
			__m256i M_0 = SIMD_ANDNOT_I256(S1_0, S2_0);  // missing value, 1 is missing
			__m256i M_1 = SIMD_ANDNOT_I256(S1_1, S2_1);  // missing value, 1 is missing
			__m256i MASK_0 = SIMD_ANDNOT_I256(M_0, (H1_0 ^ S2_0) | (H2_0 ^ S1_0));
			__m256i MASK_1 = SIMD_ANDNOT_I256(M_1, (H1_1 ^ S2_1) | (H2_1 ^ S1_1));
			__m256i va_0 = (H1_0 ^ S1_0) & MASK_0, vb_0 = (H2_0 ^ S2_0) & MASK_0;
			__m256i va_1 = (H1_1 ^ S1_1) & MASK_1, vb_1 = (H2_1 ^ S2_1) & MASK_1;
			// popcount for 64b integers
			__m256i ii4 = _mm256_popcnt_epi64(va_0) + _mm256_popcnt_epi64(vb_0) +
				_mm256_popcnt_epi64(va_1) + _mm256_popcnt_epi64(vb_1);
			// four frequencies
			__m256d f = _mm256_i64gather_pd(EXP_LOG_MIN_RARE_FREQ, ii4, 8);
			__m256d f2 = _mm256_loadu_pd(GS.p_Freq + i2);
			sum4 += _mm256_set1_pd(ff) * f2 * f;
			n -= 4;
		}
		// maybe different behavior due to rounding error of addition
		__m128d a = _mm256_castpd256_pd128(sum4) + _mm256_extractf128_pd(sum4,1);
		prob += a[0] + a[1];
	}
	return n;
}


void SIMD_NAME(CAlg_Prediction::_PrepHaploMatch)(const TGenotype &Geno,
	THaplotype *pH1_st, size_t pH1_n, THaplotype *pH2_st, size_t pH2_n,
	size_t Num_SNP, std::vector<CAlg_EM::THaploPair> &HP_PairList, short DiffList[])
{
	const TGenoStruct_512vpopcnt GS(Num_SNP, Geno);
	int MinDiff = Num_SNP * 4;
	short *pD = DiffList;

	if (pH1_st != pH2_st)
	{
		THaplotype *p1 = pH1_st;
		for (size_t n1=pH1_n; n1 > 0; n1--, p1++)
		{
			THaplotype *p2 = pH2_st;
			for (size_t n2=pH2_n; n2 > 0; n2--, p2++)
			{
				int d = hamm_d(GS, *p1, *p2);
				*pD++ = d;
				if (d < MinDiff) MinDiff = d;
				if (d == 0)
					HP_PairList.push_back(CAlg_EM::THaploPair(p1, p2));
			}
		}

		if (MinDiff > 0)
		{
			pD = DiffList;
			THaplotype *p1 = pH1_st;
			for (size_t n1=pH1_n; n1 > 0; n1--, p1++)
			{
				THaplotype *p2 = pH2_st;
				for (size_t n2=pH2_n; n2 > 0; n2--, p2++)
				{
					if (*pD++ == MinDiff)
						HP_PairList.push_back(CAlg_EM::THaploPair(p1, p2));
				}
			}
		}

	} else {
		THaplotype *p1 = pH1_st;
		for (size_t n1=pH1_n; n1 > 0; n1--, p1++)
		{
			THaplotype *p2 = p1;
			for (size_t n2=n1; n2 > 0; n2--, p2++)
			{
				int d = hamm_d(GS, *p1, *p2);
				*pD++ = d;
				if (d < MinDiff) MinDiff = d;
				if (d == 0)
					HP_PairList.push_back(CAlg_EM::THaploPair(p1, p2));
			}
		}

		if (MinDiff > 0)
		{
			pD = DiffList;
			THaplotype *p1 = pH1_st;
			for (size_t n1=pH1_n; n1 > 0; n1--, p1++)
			{
				THaplotype *p2 = p1;
				for (size_t n2=n1; n2 > 0; n2--, p2++)
				{
					if (*pD++ == MinDiff)
						HP_PairList.push_back(CAlg_EM::THaploPair(p1, p2));
				}
			}
		}
	}
}


THLAType SIMD_NAME(CAlg_Prediction::_BestGuess)(const CHaplotypeList &Haplo,
	const TGenotype &Geno)
{
	const TGenoStruct_512vpopcnt GS(Haplo, Geno);
	THLAType rv;
	rv.Allele1 = rv.Allele2 = NA_INTEGER;
	double max=0;
	const int nHLA = Haplo.nHLA();
	THaplotype *base=Haplo.List, *I1=base;

	for (int h1=0; h1 < nHLA; h1++)
	{
		const size_t n1 = Haplo.LenPerHLA[h1];

		// diagonal
		double prob = 0;
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
		THaplotype *I2 = I1 + n1;
		if (max < prob)
		{
			max = prob;
			rv.Allele1 = rv.Allele2 = h1;
		}

		// off-diagonal
		for (int h2=h1+1; h2 < nHLA; h2++)
		{
			const size_t n2 = Haplo.LenPerHLA[h2];
			prob = 0; i1 = I1;
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
	const TGenoStruct_512vpopcnt GS(Haplo, Geno);
	int H1=HLA.Allele1, H2=HLA.Allele2;
	if (H1 > H2) std::swap(H1, H2);
	const int nHLA = Haplo.nHLA();
	int IxHLA = H2 + H1*(2*nHLA-H1-1)/2;
	int idx = 0;
	double sum=0, hlaProb=0;
	THaplotype *base=Haplo.List, *I1=base;

	for (int h1=0; h1 < nHLA; h1++)
	{
		const size_t n1 = Haplo.LenPerHLA[h1];

		// diagonal
		double prob = 0;
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
		THaplotype *I2 = I1 + n1;
		if (IxHLA == idx) hlaProb = prob;
		idx ++; sum += prob;

		// off-diagonal
		for (int h2=h1+1; h2 < nHLA; h2++)
		{
			const size_t n2 = Haplo.LenPerHLA[h2];
			prob = 0; i1 = I1;
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
	const TGenoStruct_512vpopcnt GS(Haplo, Geno);
	double *p = Prob;
	const int nHLA = Haplo.nHLA();
	THaplotype *base=Haplo.List, *I1=base;

	for (int h1=0; h1 < nHLA; h1++)
	{
		const size_t n1 = Haplo.LenPerHLA[h1];

		// diagonal
		double sum = 0;
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
		THaplotype *I2 = I1 + n1;

		// off-diagonal
		for (int h2=h1+1; h2 < nHLA; h2++)
		{
			const size_t n2 = Haplo.LenPerHLA[h2];
			sum = 0; i1 = I1;
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
	double sum = 0;
	for (size_t i=0; i < n; i++) sum += Prob[i];
	const double ff = 1 / sum;
	for (size_t i=0; i < n; i++) Prob[i] *= ff;
	return sum;
}


#else

void SIMD_NAME(CAlg_Prediction::_PrepHaploMatch)(const TGenotype &Geno,
	THaplotype *pH1_st, size_t pH1_n, THaplotype *pH2_st, size_t pH2_n,
	size_t Num_SNP, std::vector<CAlg_EM::THaploPair> &HP_PairList, short DiffList[])
{
	THROW_ERROR;
}

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
