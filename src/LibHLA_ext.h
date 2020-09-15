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


#endif /* LIBHLA_EXT_H_ */
