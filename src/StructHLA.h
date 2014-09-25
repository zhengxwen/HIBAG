// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// Copyright (C) 2011-2014	Xiuwen Zheng (zhengx@u.washington.edu)
//
// This file is part of HIBAG R package.
//
// HIBAG is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// HIBAG is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with HIBAG.
// If not, see <http://www.gnu.org/licenses/>.

// ===========================================================
// Name           : LibHLA
// Author         : Xiuwen Zheng
// Version        : 1.2.4
// Copyright      : Xiuwen Zheng (GPL v3.0)
// Created        : 11/14/2011
// Last modified  : 10/29/2013
// Description    : HLA Genotype Imputation with Attribute Bagging
// ===========================================================

#ifndef _StructHLA_H_
#define _StructHLA_H_

#include <stdint.h>


#ifdef __cplusplus
extern "C" {
#endif


/// Define macro for GPU computing
#define HIBAG_ALLOW_GPU_SUPPORT


/// Define unsigned integers
typedef uint8_t   UINT8;
typedef uint16_t  UINT16;
typedef uint32_t  UINT32;
typedef uint64_t  UINT64;


/// The number of bytes or shorts for packed SNP genotypes
#define HIBAG_PACKED_NUM_IN_UTYPE(x)  (((x)/4) + (((x) % 4)!=0 ? 1 : 0))

/// The max number of SNP markers in an individual classifier
#define HIBAG_MAXNUM_SNP_IN_CLASSIFIER	256

/// The max number of bytes or shorts for packed SNP genotypes
#define HIBAG_PACKED_UTYPE_MAXNUM_SNP		\
	HIBAG_PACKED_NUM_IN_UTYPE(HIBAG_MAXNUM_SNP_IN_CLASSIFIER)




/// A pair of HLA alleles
struct THLAType
{
	int Allele1;  //< the first HLA allele
	int Allele2;  //< the second HLA allele
};



// ************************************************************************* //
// ********                       description                       ********
// SNP allele: 0 (B allele), 1 (A allele)
// SNP genotype: 0 (BB), 1 (AB), 2 (AA), 3 (missing)
// HLA allele: start from 0
//
// Packed SNP storage strategy is used for faster matching
// Packed SNP alleles: s8 s7 s6 s5 s4 s3 s2 s1
//     the 1st allele: (s1), the 2nd allele: (s3)
//     the 3rd allele: (s5), the 4th allele: (s7)
// Packed SNP genotype: s8 s7 s6 s5 s4 s3 s2 s1
//     the 1st genotype: (s2 s1), the 2nd genotype: (s4 s3)
//     the 3rd genotype: (s6 s5), the 4th genotype: (s8 s7)
// ********                                                         ********
// ************************************************************************* //


#ifdef HIBAG_ALLOW_GPU_SUPPORT

/// Packed SNP haplotype structure: 4 SNPs in a byte / short
struct TGPU_Haplotype_F32
{
	UINT8 PackedHaplo[HIBAG_PACKED_UTYPE_MAXNUM_SNP];
	float Frequency;
};

struct TGPU_Haplotype_F64
{
	UINT8 PackedHaplo[HIBAG_PACKED_UTYPE_MAXNUM_SNP];
	double Frequency;
};


/// Packed SNP genotype structure: 4 SNPs in a byte / short
struct TGPU_Genotype
{
	UINT8 PackedSNPs[HIBAG_PACKED_UTYPE_MAXNUM_SNP];
	int BootstrapCount;
	THLAType HLA;
};

#endif


#ifdef __cplusplus
}
#endif

#endif /* _StructHLA_H_ */
