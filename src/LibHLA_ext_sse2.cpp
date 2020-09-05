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
// Name           : LibHLA_ext_sse2
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


#ifdef HIBAG_CPU_ARCH_X86_SSE2
extern const bool HIBAG_ALGORITHM_SSE2 = true;
#else
extern const bool HIBAG_ALGORITHM_SSE2 = false;
#endif


#define SIMD_NAME(NAME)  NAME ## _sse2
#define THROW_ERROR      throw ErrHLA("No SSE2 support!")


#ifdef HIBAG_CPU_ARCH_X86_SSE2
#   include <xmmintrin.h>  // SSE
#   include <emmintrin.h>  // SSE2
#   if !defined(__SSE2__) && !defined(__clang__) && !defined(__ICC)
		#pragma GCC target("sse2")
#   endif


#if defined(__POPCNT__)
	#define U_POPCOUNT    __builtin_popcountll
	extern const bool HIBAG_ALGORITHM_SSE2_POPCNT = true;
#else
	extern const bool HIBAG_ALGORITHM_SSE2_POPCNT = false;
	static const __m128i Z55 = { 0x5555555555555555LL, 0x5555555555555555LL };
	static const __m128i Z33 = { 0x3333333333333333LL, 0x3333333333333333LL };
	static const __m128i Z0F = { 0x0F0F0F0F0F0F0F0FLL, 0x0F0F0F0F0F0F0F0FLL };
	#define SSE2_POPCOUNT_128B(x, cnt)    \
	{ \
		x = _mm_sub_epi64(x, _mm_srli_epi64(x, 1) & Z55); \
		x = _mm_add_epi64(x & Z33, _mm_srli_epi64(x, 2) & Z33); \
		x = _mm_add_epi64(x, _mm_srli_epi64(x, 4)) & Z0F; \
		cnt = ((uint64_t(x[0]) * 0x0101010101010101LLU) >> 56) + \
			((uint64_t(x[1]) * 0x0101010101010101LLU) >> 56); \
	}
#endif


/// Prepare the internal genotype structure
struct TGenoStruct_sse2
{
public:
	__m128i S1, S2;  ///< packed genotypes
	bool Low64b;     ///< whether length <= 64 or not
	/// constructor
	inline TGenoStruct_sse2(size_t Length, const TGenotype &G)
	{
		const INT64 *s1 = G.PackedSNP1, *s2 = G.PackedSNP2;
		Low64b = (Length <= 64);
		if (Low64b)
		{
			__m128i I1 = { s1[0], s2[0] }, I2 = { s2[0], s1[0] };  // genotypes
			S1 = I1; S2 = I2;
		} else {
			__m128i I1 = { s1[0], s1[1] }, I2 = { s2[0], s2[1] };  // genotypes
			S1 = I1; S2 = I2;
		}
	}
};


static ALWAYS_INLINE int hamm_d(const TGenoStruct_sse2 &G, const THaplotype &H1,
	const THaplotype &H2)
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
	#ifndef U_POPCOUNT
		int cnt;
		SSE2_POPCOUNT_128B(v, cnt)
		return cnt;
	#else
		// popcount
		return U_POPCOUNT(v[0]) + U_POPCOUNT(v[1]);
	#endif
	} else {
		// since HIBAG_MAXNUM_SNP_IN_CLASSIFIER = 128
		__m128i H1 = { h1[0], h1[1] }, H2 = { h2[0], h2[1] };  // two haplotypes
		__m128i S1 = G.S1, S2 = G.S2;  // genotypes
		// worry about n < UTYPE_BIT_NUM? unused bits are set to be a missing flag
		__m128i M = _mm_andnot_si128(S1, S2);  // missing value, 1 is missing
		__m128i MASK = _mm_andnot_si128(M, (H1 ^ S2) | (H2 ^ S1));
		__m128i va = (H1 ^ S1) & MASK;  // (H1 ^ S1) & MASK, (H2 ^ S2) & MASK
		__m128i vb = (H2 ^ S2) & MASK;
	#ifndef U_POPCOUNT
		int cnt1, cnt2;
		SSE2_POPCOUNT_128B(va, cnt1);
		SSE2_POPCOUNT_128B(vb, cnt2);
		return cnt1 + cnt2;
	#else
		// popcount
		return U_POPCOUNT(va[0]) + U_POPCOUNT(va[1]) +
			U_POPCOUNT(vb[0]) + U_POPCOUNT(vb[1]);
	#endif
	}
}


THLAType SIMD_NAME(CAlg_Prediction::_BestGuess)(const CHaplotypeList &Haplo,
	const TGenotype &Geno)
{
	const TGenoStruct_sse2 GS(Haplo.Num_SNP, Geno);
	THLAType rv;
	rv.Allele1 = rv.Allele2 = NA_INTEGER;
	double max=0, prob;
	const int nHLA = Haplo.nHLA();
	THaplotype *I1=Haplo.List, *I2;

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
			for (size_t m2=m1-1; m2 > 0; m2--, i2++)
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
				for (size_t m2=n2; m2 > 0; m2--, i2++)
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
	const TGenoStruct_sse2 GS(Haplo.Num_SNP, Geno);
	int H1=HLA.Allele1, H2=HLA.Allele2;
	if (H1 > H2) std::swap(H1, H2);
	const int nHLA = Haplo.nHLA();
	int IxHLA = H2 + H1*(2*nHLA-H1-1)/2;
	int idx = 0;
	double sum=0, hlaProb=0, prob;
	THaplotype *I1=Haplo.List, *I2;

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
			for (size_t m2=m1-1; m2 > 0; m2--, i2++)
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
				for (size_t m2=n2; m2 > 0; m2--, i2++)
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
	const TGenoStruct_sse2 GS(Haplo.Num_SNP, Geno);
	double *p = Prob, sum;
	const int nHLA = Haplo.nHLA();
	THaplotype *I1=Haplo.List, *I2;

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
			const double ff = 2 * i1->Freq;
			THaplotype *i2 = i1 + 1;
			for (size_t m2=m1-1; m2 > 0; m2--, i2++)
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
				for (size_t m2=n2; m2 > 0; m2--, i2++)
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
