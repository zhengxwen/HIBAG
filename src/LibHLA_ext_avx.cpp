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
#   define HIBAG_CPU_ARCH_X86_AVX
#endif

#if defined(__AVX__) && !defined(HIBAG_CPU_ARCH_X86_AVX)
#   define HIBAG_CPU_ARCH_X86_AVX
#endif

#ifdef HIBAG_CPU_ARCH_X86_AVX
extern const bool HIBAG_ALGORITHM_AVX = true;
#else
extern const bool HIBAG_ALGORITHM_AVX = false;
#endif


#define SIMD_NAME(NAME)  NAME ## _avx
#define THROW_ERROR      throw ErrHLA("No AVX support!")


#ifdef HIBAG_CPU_ARCH_X86_AVX

#   include <xmmintrin.h>  // SSE
#   include <emmintrin.h>  // SSE2
#   include <immintrin.h>  // AVX, AVX2
#   ifndef __AVX__
		#pragma GCC target("avx")
#   endif


typedef int64_t UTYPE;


#define U_POPCOUNT    __builtin_popcountll


#   define GENO_HAMM_DIST_INIT(Length)    \
		TGenoHammDist GenoVar; init_hamm_dist(Length, Geno, GenoVar)
#   define GENO_VAR  GenoVar
#   define GENO_TYPE TGenoHammDist
#   define GENO_HALF_NBIT  64
typedef union {
	__m128i i128;
	__m256i i256;
} t_simd;
typedef struct {
	t_simd S1, S2;  ///< packed genotypes
} TGenoHammDist;


static ALWAYS_INLINE void init_hamm_dist(size_t Length, const TGenotype &G,
	TGenoHammDist &out)
{
	const UTYPE *s1 = (const UTYPE*)&G.PackedSNP1[0];
	const UTYPE *s2 = (const UTYPE*)&G.PackedSNP2[0];
	if (Length <= GENO_HALF_NBIT)
	{
		__m128i S1 = { *s1, *s2 }, S2 = { *s2, *s1 };  // genotypes
		out.S1.i128 = S1; out.S2.i128 = S2;
	} else {
		__m256i S1 = { s1[0], s1[1], s2[0], s2[1] };  // genotypes
		__m256i S2 = { s2[0], s2[1], s1[0], s1[1] };  // genotypes
		out.S1.i256 = S1; out.S2.i256 = S2;
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
		__m128i S1 = G.S1.i128, S2 = G.S2.i128;  // genotypes
		__m128i m1 = H ^ S2, m2 = { m1[1], m1[0] };
		// worry about n < UTYPE_BIT_NUM? unused bits are set to be a missing flag
		__m128i M = S2 & ~S1;  // missing value, 1 is missing
		__m128i M2 = { M[0], M[0] };
		__m128i MASK = (m1 | m2) & ~M2;
		__m128i v = (H ^ S1) & MASK;  // (H1 ^ S1) & MASK, (H2 ^ S2) & MASK
		// popcount
		return U_POPCOUNT(v[0]) + U_POPCOUNT(v[1]);
	} else {
		// since HIBAG_MAXNUM_SNP_IN_CLASSIFIER = 128
		__m256i H  = { h1[0], h1[1], h2[0], h2[1] };  // two haplotypes
		__m256i S1 = G.S1.i256, S2 = G.S2.i256;  // genotypes
		__m256i m1 = H ^ S2, m2 = { m1[2], m1[3], m1[0], m1[1] };
		// worry about n < UTYPE_BIT_NUM? unused bits are set to be a missing flag
		__m256i M = S2 & ~S1;  // missing value, 1 is missing
		__m256i M2 = { M[0], M[1], M[0], M[1] };
		__m256i MASK = (m1 | m2) & ~M2;
		__m256i v = (H ^ S1) & MASK;  // (H1 ^ S1) & MASK, (H2 ^ S2) & MASK
		// popcount
		return U_POPCOUNT(v[0]) + U_POPCOUNT(v[1]) +
			U_POPCOUNT(v[2]) + U_POPCOUNT(v[3]);
	}
}


THLAType SIMD_NAME(CAlg_Prediction::_BestGuess)(const CHaplotypeList &Haplo,
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

double SIMD_NAME(CAlg_Prediction::_PostProb)(const CHaplotypeList &Haplo,
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

void SIMD_NAME(CAlg_Prediction::_PostProb2)(const CHaplotypeList &Haplo,
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
