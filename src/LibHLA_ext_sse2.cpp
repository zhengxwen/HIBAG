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


#ifdef HIBAG_CPU_ARCH_X86
extern const bool HIBAG_ALGORITHM_SSE2 = true;
#else
extern const bool HIBAG_ALGORITHM_SSE2 = false;
#endif


#ifdef HIBAG_CPU_ARCH_X86
#   include <xmmintrin.h>  // SSE
#   include <emmintrin.h>  // SSE2
#   ifndef __SSE2__
		#pragma GCC target("sse2")
#   endif


typedef int64_t UTYPE;


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


#define GENO_HAMM_DIST_INIT(Length)    \
	TGenoHammDist GenoVar; init_hamm_dist(Length, Geno, GenoVar)
#define GENO_VAR  GenoVar
#define GENO_TYPE TGenoHammDist
#define GENO_HALF_NBIT  64

typedef struct {
	__m128i S1a, S2a;  ///< packed genotypes
	__m128i S1b, S2b;  ///< packed genotypes
} TGenoHammDist;


static ALWAYS_INLINE void init_hamm_dist(size_t Length, const TGenotype &G,
	TGenoHammDist &out)
{
	const UTYPE *s1 = (const UTYPE*)&G.PackedSNP1[0];
	const UTYPE *s2 = (const UTYPE*)&G.PackedSNP2[0];
	if (Length <= GENO_HALF_NBIT)
	{
		__m128i S1 = { *s1, *s2 }, S2 = { *s2, *s1 };  // genotypes
		out.S1a = S1; out.S2a = S2;
	} else {
		__m128i S1a = { s1[0], s1[1] }, S1b = { s2[0], s2[1] };  // genotypes
		__m128i S2a = { s2[0], s2[1] }, S2b = { s1[0], s1[1] };  // genotypes
		out.S1a = S1a; out.S2a = S2a;
		out.S1b = S1b; out.S2b = S2b;
	}
}


static ALWAYS_INLINE int hamm_dist(size_t Length, const GENO_TYPE &G,
	const THaplotype &H1, const THaplotype &H2)
{
	const UTYPE *h1 = (const UTYPE*)&H1.PackedHaplo[0];
	const UTYPE *h2 = (const UTYPE*)&H2.PackedHaplo[0];
	// here, UTYPE = int64_t
	if (Length <= GENO_HALF_NBIT)
	{
		__m128i H  = { *h1, *h2 };  // two haplotypes
		__m128i S1 = G.S1a, S2 = G.S2a;  // genotypes
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
		__m128i Ha = { h1[0], h1[1] }, Hb = { h2[0], h2[1] };  // two haplotypes
		__m128i S1a = G.S1a, S1b = G.S1b;
		__m128i S2a = G.S2a, S2b = G.S2b;
		__m128i m1a = Ha ^ S2a, m1b = Hb ^ S2b;
		__m128i m2a = m1b, m2b = m1a;
		// worry about n < UTYPE_BIT_NUM? unused bits are set to be a missing flag
		__m128i M2 = _mm_andnot_si128(S1a, S2a);  // missing value, 1 is missing
		__m128i MASKa = _mm_andnot_si128(M2, m1a | m2a);
		__m128i MASKb = _mm_andnot_si128(M2, m1b | m2b);
		__m128i va = (Ha ^ S1a) & MASKa;  // (H1 ^ S1) & MASK, (H2 ^ S2) & MASK
		__m128i vb = (Hb ^ S1b) & MASKb;
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


THLAType CAlg_Prediction::_BestGuess_sse2(const CHaplotypeList &Haplo,
	const TGenotype &Geno)
{
	THLAType rv;
	rv.Allele1 = rv.Allele2 = NA_INTEGER;
	double max=0, prob;
	THaplotype *I1, *I2;
	I1 = Haplo.List;
	GENO_HAMM_DIST_INIT(Haplo.Num_SNP);

	for (int h1=0; h1 < _nHLA; h1++)
	{
		size_t n1 = Haplo.LenPerHLA[h1];

		// diagonal
		prob = 0;
		THaplotype *i1 = I1;
		for (size_t m1=n1; m1 > 0; m1--, i1++)
		{
			THaplotype *i2 = i1;
			for (size_t m2=m1; m2 > 0; m2--, i2++)
			{
				ADD_FREQ_MUTANT(prob, (i1 != i2) ?
					(2 * i1->Freq * i2->Freq) : (i1->Freq * i2->Freq),
					hamm_dist(Haplo.Num_SNP, GENO_VAR, *i1, *i2));
			}
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
			size_t n2 = Haplo.LenPerHLA[h2];
			prob = 0;
			THaplotype *i1 = I1;
			for (size_t m1=n1; m1 > 0; m1--, i1++)
			{
				THaplotype *i2 = I2;
				for (size_t m2=n2; m2 > 0; m2--, i2++)
				{
					ADD_FREQ_MUTANT(prob, 2 * i1->Freq * i2->Freq,
						hamm_dist(Haplo.Num_SNP, GENO_VAR, *i1, *i2));
				}
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

double CAlg_Prediction::_PostProb_sse2(const CHaplotypeList &Haplo,
	const TGenotype &Geno, const THLAType &HLA)
{
	int H1=HLA.Allele1, H2=HLA.Allele2;
	if (H1 > H2) std::swap(H1, H2);
	int IxHLA = H2 + H1*(2*_nHLA-H1-1)/2;
	int idx = 0;
	double sum=0, hlaProb=0, prob;
	THaplotype *I1, *I2;
	I1 = Haplo.List;
	GENO_HAMM_DIST_INIT(Haplo.Num_SNP);

	for (int h1=0; h1 < _nHLA; h1++)
	{
		size_t n1 = Haplo.LenPerHLA[h1];

		// diagonal
		prob = 0;
		THaplotype *i1 = I1;
		for (size_t m1=n1; m1 > 0; m1--, i1++)
		{
			THaplotype *i2 = i1;
			for (size_t m2=m1; m2 > 0; m2--, i2++)
			{
				ADD_FREQ_MUTANT(prob, (i1 != i2) ?
					(2 * i1->Freq * i2->Freq) : (i1->Freq * i2->Freq),
					hamm_dist(Haplo.Num_SNP, GENO_VAR, *i1, *i2));
			}
		}
		I2 = I1 + n1;
		if (IxHLA == idx) hlaProb = prob;
		idx ++; sum += prob;

		// off-diagonal
		for (int h2=h1+1; h2 < _nHLA; h2++)
		{
			size_t n2 = Haplo.LenPerHLA[h2];
			prob = 0;
			THaplotype *i1 = I1;
			for (size_t m1=n1; m1 > 0; m1--, i1++)
			{
				THaplotype *i2 = I2;
				for (size_t m2=n2; m2 > 0; m2--, i2++)
				{
					ADD_FREQ_MUTANT(prob, 2 * i1->Freq * i2->Freq,
						hamm_dist(Haplo.Num_SNP, GENO_VAR, *i1, *i2));
				}
			}
			I2 += n2;
			if (IxHLA == idx) hlaProb = prob;
			idx ++; sum += prob;
		}

		I1 += n1;
	}

	return hlaProb / sum;
}

void CAlg_Prediction::_PostProb2_sse2(const CHaplotypeList &Haplo,
	const TGenotype &Geno, double &SumProb)
{
	THaplotype *I1, *I2;
	double *pProb = &_PostProb[0];
	double sum;
	GENO_HAMM_DIST_INIT(Haplo.Num_SNP);

	I1 = Haplo.List;
	for (int h1=0; h1 < _nHLA; h1++)
	{
		size_t n1 = Haplo.LenPerHLA[h1];

		// diagonal
		sum = 0;
		THaplotype *i1 = I1;
		for (size_t m1=n1; m1 > 0; m1--, i1++)
		{
			THaplotype *i2 = i1;
			for (size_t m2=m1; m2 > 0; m2--, i2++)
			{
				ADD_FREQ_MUTANT(sum, (i1 != i2) ?
					(2 * i1->Freq * i2->Freq) : (i1->Freq * i2->Freq),
					hamm_dist(Haplo.Num_SNP, GENO_VAR, *i1, *i2));
			}
		}
		*pProb++ = sum;
		I2 = I1 + n1;

		// off-diagonal
		for (int h2=h1+1; h2 < _nHLA; h2++)
		{
			size_t n2 = Haplo.LenPerHLA[h2];
			sum = 0;
			THaplotype *i1 = I1;
			for (size_t m1=n1; m1 > 0; m1--, i1++)
			{
				THaplotype *i2 = I2;
				for (size_t m2=n2; m2 > 0; m2--, i2++)
				{
					ADD_FREQ_MUTANT(sum, 2 * i1->Freq * i2->Freq,
						hamm_dist(Haplo.Num_SNP, GENO_VAR, *i1, *i2));
				}
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

#endif
