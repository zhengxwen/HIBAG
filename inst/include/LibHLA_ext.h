// ===============================================================
//
// HIBAG R package (HLA Genotype Imputation with Attribute Bagging)
// Copyright (C) 2020-2021   Xiuwen Zheng (zhengx@u.washington.edu)
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
// Name           : LibHLA ext header
// Author         : Xiuwen Zheng
// Kernel Version : 1.5
// Copyright      : Xiuwen Zheng (GPL v3)
// Description    : Define Macro for CPU intrinsics
// ===============================================================

#ifndef LIBHLA_EXT_H_
#define LIBHLA_EXT_H_


// Detect whether x86 microprocessor architecture or not
#if defined(__i386__) || defined(__X86__) || defined(_M_IX86) || defined(__I86__) || defined(__INTEL__) || defined(__amd64__) || defined(__x86_64__) || defined(_M_AMD64)
#   define HIBAG_CPU_ARCH_X86
#endif


// 32-bit or 64-bit registers
#ifdef __LP64__
#   define HIBAG_CPU_LP64
#else
#   ifdef HIBAG_CPU_LP64
#      undef HIBAG_CPU_LP64
#   endif
#endif


// Always inline function
#ifndef ALWAYS_INLINE
#   ifdef _MSC_VER
#       define ALWAYS_INLINE    __forceinline
#   elif defined(__GNUC__)
#       define ALWAYS_INLINE    inline __attribute__((always_inline))
#   else
#       define ALWAYS_INLINE    inline
#   endif
#endif


// remove HIBAG_CPU_ARCH_X86 if define HIBAG_NO_X86_SIMD
#if defined(HIBAG_CPU_ARCH_X86) && defined(HIBAG_NO_X86_SIMD)
#   undef HIBAG_CPU_ARCH_X86
#endif

// whether has __builtin_cpu_supports or not
#if defined(__GNUC__) && ((__GNUC__>4) || (__GNUC__==4 && __GNUC_MINOR__>=8))
#   define HIBAG_BUILTIN_CPU
#elif defined(__clang__) && defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=7))
#   define HIBAG_BUILTIN_CPU
#endif



// SSE2
#ifdef HIBAG_CPU_ARCH_X86
#   if defined(__GNUC__) && ((__GNUC__>3) || (__GNUC__==3 && __GNUC_MINOR__>=1))
#       define HIBAG_CPU_ARCH_X86_SSE2
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=3))
#       define HIBAG_CPU_ARCH_X86_SSE2
#   endif
#   ifdef HIBAG_BUILTIN_CPU
#       define HIBAG_BUILTIN_CPU_SSE2
#   endif
#endif
#if defined(__SSE2__) && !defined(HIBAG_CPU_ARCH_X86_SSE2)
#   define HIBAG_CPU_ARCH_X86_SSE2
#endif


// SSE4.2
#ifdef HIBAG_CPU_ARCH_X86
#   if defined(__GNUC__) && ((__GNUC__>4) || (__GNUC__==4 && __GNUC_MINOR__>=3))
#       define HIBAG_CPU_ARCH_X86_SSE4_2
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=3))
#       define HIBAG_CPU_ARCH_X86_SSE4_2
#   endif
#   ifdef HIBAG_BUILTIN_CPU
#       define HIBAG_BUILTIN_CPU_SSE4_2
#   endif
#endif
#if defined(__SSE4_2__) && !defined(HIBAG_CPU_ARCH_X86_SSE4_2)
#   define HIBAG_CPU_ARCH_X86_SSE4_2
#endif


// AVX
#ifdef HIBAG_CPU_ARCH_X86
#   if defined(__GNUC__) && ((__GNUC__>4) || (__GNUC__==4 && __GNUC_MINOR__>=4))
#       define HIBAG_CPU_ARCH_X86_AVX
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=3))
#       define HIBAG_CPU_ARCH_X86_AVX
#   endif
#   ifdef HIBAG_BUILTIN_CPU
#       define HIBAG_BUILTIN_CPU_AVX
#   endif
#endif
#if defined(__AVX__) && !defined(HIBAG_CPU_ARCH_X86_AVX)
#   define HIBAG_CPU_ARCH_X86_AVX
#endif


// AVX2
#ifdef HIBAG_CPU_ARCH_X86
#   if defined(__GNUC__) && ((__GNUC__>4) || (__GNUC__==4 && __GNUC_MINOR__>=7))
#       define HIBAG_CPU_ARCH_X86_AVX2
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=7))
#       define HIBAG_CPU_ARCH_X86_AVX2
#   endif
#   ifdef HIBAG_BUILTIN_CPU
#       define HIBAG_BUILTIN_CPU_AVX2
#   endif
#endif
#if defined(__AVX2__) && !defined(HIBAG_CPU_ARCH_X86_AVX2)
#   define HIBAG_CPU_ARCH_X86_AVX2
#endif


// AVX512F
#ifdef HIBAG_CPU_ARCH_X86
#   if defined(__GNUC__) && (__GNUC__>=5)
#       define HIBAG_CPU_ARCH_X86_AVX512F
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=9))
#       define HIBAG_CPU_ARCH_X86_AVX512F
#   endif
#   if defined(__GNUC__) && (__GNUC__>=6)
#       define HIBAG_BUILTIN_CPU_AVX512F
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=9))
#       define HIBAG_BUILTIN_CPU_AVX512F
#   elif defined(__ICC)
#       define HIBAG_BUILTIN_CPU_AVX512F
#   endif
#endif
#if defined(__AVX512F__) && !defined(HIBAG_CPU_ARCH_X86_AVX512F)
#   define HIBAG_CPU_ARCH_X86_AVX512F
#endif


// AVX512BW
#ifdef HIBAG_CPU_ARCH_X86
#   if defined(__GNUC__) && (__GNUC__>=5)
#       define HIBAG_CPU_ARCH_X86_AVX512BW
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=9))
#       define HIBAG_CPU_ARCH_X86_AVX512BW
#   endif
#   if defined(__GNUC__) && (__GNUC__>=6)
#       define HIBAG_BUILTIN_CPU_AVX512BW
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=9))
#       define HIBAG_BUILTIN_CPU_AVX512BW
#   elif defined(__ICC)
#       define HIBAG_BUILTIN_CPU_AVX512BW
#   endif
#endif
#if defined(__AVX512BW__) && !defined(HIBAG_CPU_ARCH_X86_AVX512BW)
#   define HIBAG_CPU_ARCH_X86_AVX512BW
#endif


#include <stdint.h>
#include <stdlib.h>

#ifndef HIBAG_STRUCTURE_HEAD_ONLY
#   include <string>
#endif

#ifdef __cplusplus
namespace HLA_LIB
{
	/// Kernel Version, Major Number (0x01) / Minor Number (0x05)
	#define HIBAG_KERNEL_VERSION    0x0105

	/// Define 8-bit and 64-bit integers
	typedef uint8_t     UINT8;
	typedef int64_t     INT64;

	/// The max number of SNP markers in an individual classifier.
	//  Don't modify this value since the code is optimized for this value!!!
	const size_t HIBAG_MAXNUM_SNP_IN_CLASSIFIER = 128;

	/// The max number of INT64 for packed SNP genotypes.
	const size_t HIBAG_PACKED_INT64_MAXNUM =
		HIBAG_MAXNUM_SNP_IN_CLASSIFIER / (8*sizeof(INT64));


	// ===================================================================== //
	// ========                     Description                     ========
	//
	// Packed SNP storage strategy is used for faster matching
	//
	// HLA allele: start from 0
	//
	// THaplotype: packed SNP alleles (little endianness):
	//     (s8 s7 s6 s5 s4 s3 s2 s1)
	//     the 1st allele: (s1), the 2nd allele: (s2), ...
	//     SNP allele: 0 (B allele), 1 (A allele)
	//
	// TGenotype: packed SNP genotype (little endianness):
	//     array_1 = (s1_8 s1_7 s1_6 s1_5 s1_4 s1_3 s1_2 s1_1),
	//     array_2 = (s2_8 s2_7 s2_6 s2_5 s2_4 s2_3 s2_2 s2_1),
	//     array_3 = (s3_8 s3_7 s3_6 s3_5 s3_4 s3_3 s3_2 s3_1)
	//     the 1st genotype: (s1_1 s2_1 s3_1),
	//     the 2nd genotype: (s1_1 s2_1 s3_1), ...
	//     SNP genotype: 0 (BB) -- (s1_1=0 s2_1=0 s3_1=1),
	//                   1 (AB) -- (s1_1=1 s2_1=0 s3_1=1),
	//                   2 (AA) -- (s1_1=1 s2_1=1 s3_1=1),
	//                   -1 or other value (missing)
	//                          -- (s1_1=0 s2_1=1 s3_1=0)
	//
	// ========                                                     ========
	// ===================================================================== //

	/// Packed bi-allelic SNP haplotype structure: 8 alleles in a byte
	struct THaplotype
	{
	public:
		friend class CHaplotypeList;

		/// packed SNP alleles
		INT64 PackedHaplo[HIBAG_PACKED_INT64_MAXNUM];
		/// haplotype frequency
		double Freq;
		/// auxiliary variables, sizeof(THaplotype)=32
		union type_aux
		{
			double OldFreq;  /// old haplotype frequency
			struct type_aux2 {
				float Freq_f32;  /// 32-bit haplotype frequency, used in GPU
				int HLA_allele;  /// the associated HLA allele, used in GPU
			} a2;
		} aux;

    #ifndef HIBAG_STRUCTURE_HEAD_ONLY
		/// constructor
		THaplotype();
		THaplotype(double _freq);
		THaplotype(const char *str, double _freq);

		/// get SNP allele, idx starts from ZERO
		UINT8 GetAllele(size_t idx) const;
		/// set SNP allele, idx starts from ZERO
		void SetAllele(size_t idx, UINT8 val);
		/// get a string of "0" and "1" from packed SNP alleles
		std::string HaploToStr(size_t Length) const;
		/// set packed SNP alleles from a string of "0" and "1"
		void StrToHaplo(const std::string &str);

	private:
		/// set SNP allele, idx starts from ZERO, without checking
		inline void _SetAllele(size_t idx, UINT8 val);
    #endif
	};


	/// An unordered pair of HLA alleles
	struct THLAType
	{
		int Allele1;  //< the first HLA allele
		int Allele2;  //< the second HLA allele
	};


	/// Packed bi-allelic SNP genotype structure: 8 SNPs in a byte, sizeof(TGenotype)=48
	struct TGenotype
	{
	public:
		friend class CGenotypeList;
		friend class CAlg_EM;
		friend class CAlg_Prediction;

		/// packed SNP genotypes, allele 1
		INT64 PackedSNP1[HIBAG_PACKED_INT64_MAXNUM];
		/// packed SNP genotypes, allele 2
		INT64 PackedSNP2[HIBAG_PACKED_INT64_MAXNUM];
		/// the count in the bootstrapped data
		int BootstrapCount;

		/// auxiliary correct HLA type
		THLAType aux_hla_type;
		/// auxiliary integer to make sizeof(TGenotype)=48
		int aux_temp;

    #ifndef HIBAG_STRUCTURE_HEAD_ONLY
		/// constructor
		TGenotype();
		/// get SNP genotype (0, 1, 2) at the specified locus, idx starts from ZERO
		int GetSNP(size_t idx) const;
		/// set SNP genotype (0, 1, 2) at the specified locus, idx starts from ZERO
		void SetSNP(size_t idx, int val);
		/// get a string of SNP genotypes, consisting of "0" or "1"
		std::string SNPToString(size_t Length) const;
		/// set SNP genotypes by a string of "0" and "1"
		void StringToSNP(const std::string &str);
		/// export SNPs to a vector of integers
		void SNPToInt(size_t Length, int OutArray[]) const;
		/// import SNPs from an integer vector 'GenoBase' with 'Index'
		void IntToSNP(size_t Length, const int GenoBase[], const int Index[]);
		/// compute the Hamming distance between SNPs and H1+H2
		int HammingDistance(size_t Length, const THaplotype &H1, const THaplotype &H2) const;

	protected:
		/// set SNP genotype (0, 1, 2) without checking
		void _SetSNP(size_t idx, int val);
	#endif
	};


	// ===================================================================== //

	/// the data structure of functions using GPU
	struct TypeGPUExtProc
	{
		/// initialize the internal structure for building a model
		void (*build_init)(int nHLA, int nSample);
		/// finalize the structure for building a model
		void (*build_done)();
		/// initialize bottstrapping
		void (*build_set_bootstrap)(const int oob_cnt[]);
		/// initialize haplotypes and SNPs genotypes
		void (*build_set_haplo_geno)(const THaplotype haplo[], int n_haplo,
			const TGenotype geno[], int n_snp);
		/// calculate the out-of-bag accuracy (the number of correct alleles)
		int (*build_acc_oob)();
		/// calculate the in-bag log likelihood
		double (*build_acc_ib)();

		/// initialize the internal structure for predicting
		//    nHaplo[2*nClassifier]:
		//    nHaplo[0] -- total # of haplotypes
		//    nHaplo[1] -- # of SNPs
		void (*predict_init)(int nHLA, int nClassifier,
			const THaplotype *const pHaplo[], const int nHaplo[], const int nSNP[]);
		/// finalize the structure for predicting
		void (*predict_done)();
		/// average the posterior probabilities among classifiers for predicting
		void (*predict_avg_prob)(const TGenotype geno[], const double weight[],
			double out_prob[], double out_match[]);
	};

}
#endif


#endif /* LIBHLA_EXT_H_ */
