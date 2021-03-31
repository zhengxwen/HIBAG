// ===============================================================
//
// HIBAG R package (HLA Genotype Imputation with Attribute Bagging)
// Copyright (C) 2011-2021   Xiuwen Zheng (zhengx@u.washington.edu)
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
#include <Rmath.h>
#include <Rinternals.h>
#ifdef length
#   undef length
#endif
#include <RcppParallel.h>
#include <tbb/parallel_for.h>
#include <mutex>
#include <algorithm>

#ifdef HIBAG_CPU_ARCH_X86
#   include <xmmintrin.h>  // SSE
#   include <emmintrin.h>  // SSE2
#endif


using namespace std;
using namespace HLA_LIB;


// whether compile the algorithm with specified targets or not
extern const bool HIBAG_ALGORITHM_SSE2;
extern const bool HIBAG_ALGORITHM_SSE2_POPCNT;
extern const bool HIBAG_ALGORITHM_SSE4_2;
extern const bool HIBAG_ALGORITHM_AVX;
extern const bool HIBAG_ALGORITHM_AVX2;
extern const bool HIBAG_ALGORITHM_AVX512F;
extern const bool HIBAG_ALGORITHM_AVX512BW;

// target-specific functions
static CAlg_Prediction::F_PrepHaploMatch fc_PrepHaploMatch = NULL;
static CAlg_Prediction::F_BestGuess fc_BestGuess = NULL;
static CAlg_Prediction::F_PostProb  fc_PostProb  = NULL;
static CAlg_Prediction::F_PostProb2 fc_PostProb2 = NULL;


// indicate CPU flags and selection of target-specific functions
static string HIBAG_CPU_Info;

/// Get CPU information
const char *HLA_LIB::CPU_Info()
{
	return HIBAG_CPU_Info.c_str();
}


// ========================================================================= //

// Parameters -- EM algorithm

/// the max number of iterations
int HLA_LIB::EM_MaxNum_Iterations = 500;
/// the initial value of EM algorithm
static const double EM_INIT_VAL_FRAC = 0.001;
/// the reltol convergence tolerance, sqrt(machine.epsilon) by default, used in EM algorithm
double HLA_LIB::EM_FuncRelTol = sqrt(DBL_EPSILON);


// Parameters -- reduce the number of possible haplotypes

/// The fraction of one haplotype that can be ignored
static const double FRACTION_HAPLO = 1.0/10;


// Parameters -- search SNP markers

/// the reltol for the stopping rule of adding a new SNP marker
static const double STOP_RELTOL_LOGLIK_ADDSNP = 0.001;
/// the reltol for erasing the SNP marker if prune = TRUE
static const double PRUNE_RELTOL_LOGLIK = 0.1;


/// Random number: return an integer from 0 to n-1 with equal probability
static inline int RandomNum(int n)
{
	// 'unif_rand()' returns [0 .. 1]
	int v = (int)(n * unif_rand());
	if (v >= n) v = n - 1;
	return v;
}


// ========================================================================= //
// Define Intel TBB Macro, requiring C++11

#if RCPP_PARALLEL_USE_TBB

#define PARALLEL_FOR(i, SIZE)    \
	tbb::parallel_for(tbb::blocked_range<size_t>(0, SIZE),  \
		[&](const tbb::blocked_range<size_t> &rng)  \
	{  \
		for (size_t i=rng.begin(); i < rng.end(); i++)

#define PARALLEL_FOR_CHUNK(i, SIZE, CHUNK)    \
	tbb::parallel_for(tbb::blocked_range<size_t>(0, SIZE, CHUNK),  \
		[&](const tbb::blocked_range<size_t> &rng)  \
	{  \
		for (size_t i=rng.begin(); i < rng.end(); i++)

#define PARALLEL_END    });

static inline int thread_num() { return tbb::this_task_arena::max_concurrency(); }
static inline int thread_idx() { return tbb::this_task_arena::current_thread_index(); }

#else

#define PARALLEL_FOR(i, SIZE)                 {  for (size_t i=0; i < SIZE; i++)
#define PARALLEL_FOR_CHUNK(i, SIZE, CHUNK)    {  for (size_t i=0; i < SIZE; i++)
#define PARALLEL_END    }

static inline int thread_num() { return 1; }
static inline int thread_idx() { return 0; }

#endif



// ========================================================================= //

/// exp(cnt * log(MIN_RARE_FREQ)), cnt is the hamming distance
double HLA_LIB::EXP_LOG_MIN_RARE_FREQ[HIBAG_MAXNUM_SNP_IN_CLASSIFIER*2 + 1];

class CInit
{
public:
	CInit()
	{
		// initialize internal variables
		const int n = 2 * HIBAG_MAXNUM_SNP_IN_CLASSIFIER;
		for (int i=0; i <= n; i++)
			EXP_LOG_MIN_RARE_FREQ[i] = exp(i * log(MIN_RARE_FREQ));
		EXP_LOG_MIN_RARE_FREQ[0] = 1;
		for (int i=0; i <= n; i++)
		{
			if (!R_FINITE(EXP_LOG_MIN_RARE_FREQ[i]))
				EXP_LOG_MIN_RARE_FREQ[i] = 0;
		}
		// select CPU target
		CAlg_Prediction::Init_Target_IFunc("max");
	}
};

static CInit _Init;


// GPU extensible component
TypeGPUExtProc *HLA_LIB::GPUExtProcPtr = NULL;;



// ========================================================================= //
// CdProgression

static void check_interrupt_fc(void *p) { R_CheckUserInterrupt(); }

// this will call the above in a top-level context so it won't longjmp-out of your context
inline void CheckInterrupt()
{
	if (R_ToplevelExec(check_interrupt_fc, NULL) == FALSE)
		throw ErrHLA("User interrupts the progress.");
}

static char date_buffer[256];

const char *HLA_LIB::date_text()
{
	time_t rawtime;
	time(&rawtime);
	struct tm *p = localtime(&rawtime);
	sprintf(date_buffer, "%04d-%02d-%02d %02d:%02d:%02d", p->tm_year+1900,
		p->tm_mon+1, p->tm_mday, p->tm_hour, p->tm_min, p->tm_sec);
	return date_buffer;
}

static const clock_t TimeInterval = 15*CLOCKS_PER_SEC;

static mutex progress_add;

CdProgression::CdProgression()
{
	Init(0, false);
}

void CdProgression::Init(INT64 TotalCnt, bool ShowInit)
{
	if (TotalCnt < 0) TotalCnt = 0;
	fTotal = TotalCnt;
	fCurrent = fPercent = 0;
	OldTime = clock();
	if (ShowInit) ShowProgress();
}

bool CdProgression::Forward(INT64 step, bool Show)
{
	progress_add.lock();
	fCurrent += step;
	int p = int(double(TotalPercent)*fCurrent / fTotal);
	if ((p != fPercent) || (p == TotalPercent))
	{
		clock_t Now = clock();
		if (((Now - OldTime) >= TimeInterval) || (p == TotalPercent))
		{
			fPercent = p;
			if (Show) ShowProgress();
			OldTime = Now;
			progress_add.unlock();
			return true;
		}
	}
	progress_add.unlock();
	return false;
}

void CdProgression::ShowProgress()
{
	Rprintf("%s (%s)\t%d%%\n", Info.c_str(), date_text(),
		int(fPercent*StepPercent));
}


/// The progression information
CdProgression HLA_LIB::Progress;



// ========================================================================= //
// Packed bi-allelic SNP haplotype structure: 8 alleles in a byte

THaplotype::THaplotype()
{
	Freq = aux.OldFreq = 0;
}

THaplotype::THaplotype(double _freq)
{
	Freq = _freq;
	aux.OldFreq = 0;
}

THaplotype::THaplotype(const char *str, double _freq)
{
	Freq = _freq;
	aux.OldFreq = 0;
	StrToHaplo(str);
}

UINT8 THaplotype::GetAllele(size_t idx) const
{
	HIBAG_CHECKING(idx >= HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"THaplotype::GetAllele, invalid index.");
	const UINT8 *p = (const UINT8*)PackedHaplo;
	return (p[idx >> 3] >> (idx & 0x07)) & 0x01;
}

void THaplotype::SetAllele(size_t idx, UINT8 val)
{
	HIBAG_CHECKING(idx >= HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"THaplotype::SetAllele, invalid index.");
	HIBAG_CHECKING(val!=0 && val!=1,
		"THaplotype::SetAllele, the value should be 0 or 1.");
	_SetAllele(idx, val);
}

string THaplotype::HaploToStr(size_t Length) const
{
	HIBAG_CHECKING(Length > HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"THaplotype::HaploToStr, the length is invalid.");
	string rv;
	if (Length > 0)
	{
		rv.resize(Length);
		const UINT8 *p = (const UINT8*)PackedHaplo;
		for (size_t i=0; i < Length; i++)
			rv[i] = ((p[i >> 3] >> (i & 0x07)) & 0x01) ? '1' : '0';
	}
	return rv;
}

void THaplotype::StrToHaplo(const string &str)
{
	HIBAG_CHECKING(str.size() > HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"THaplotype::StrToHaplo, the input string is too long.");
	for (size_t i=0; i < str.size(); i++)
	{
		const char ch = str[i];
		HIBAG_CHECKING(ch!='0' && ch!='1',
			"THaplotype::StrToHaplo, the input string should be '0' or '1'");
		_SetAllele(i, ch-'0');
	}
}

inline void THaplotype::_SetAllele(size_t idx, UINT8 val)
{
	size_t r = idx & 0x07;
	UINT8 mask = ~(0x01 << r);
	UINT8 &ch = ((UINT8*)PackedHaplo)[idx >> 3];
	ch = (ch & mask) | (val << r);
}



// ========================================================================= //
// Haplotype list with an HLA gene and SNP alleles

CHaplotypeList::CHaplotypeList()
{
	Num_Haplo = Num_SNP = 0;
	reserve_size = 0; base_ptr = NULL;
	List = NULL;
	aux_haplo = NULL; aux_freq = NULL;
}

CHaplotypeList::CHaplotypeList(const CHaplotypeList &src)
{
	Num_Haplo = Num_SNP = 0;
	reserve_size = 0; base_ptr = NULL;
	List = NULL;
	aux_haplo = NULL; aux_freq = NULL;
	*this = src;
}

CHaplotypeList::CHaplotypeList(size_t reserve_num)
{
	Num_Haplo = Num_SNP = 0;
	reserve_size = 0; base_ptr = NULL;
	List = NULL;
	aux_haplo = NULL; aux_freq = NULL;
	if (reserve_num > 0)
		alloc_mem(reserve_size = reserve_num);
}

CHaplotypeList::~CHaplotypeList()
{
	if (base_ptr) free(base_ptr);
	base_ptr = NULL;
}

void CHaplotypeList::alloc_mem(size_t num)
{
	const size_t size = sizeof(THaplotype) * num + 32;
	base_ptr = realloc(base_ptr, size);
	if (base_ptr == NULL)
		throw ErrHLA("Fails to allocate memory.");
	UINT8 *p = (UINT8 *)base_ptr;
	size_t r = (size_t)p & 0x1F;
	if (r > 0) p += 32 - r;
	List = (THaplotype *)p;
}

CHaplotypeList& CHaplotypeList::operator= (const CHaplotypeList &src)
{
	Num_SNP = src.Num_SNP;
	LenPerHLA = src.LenPerHLA;
	ResizeHaplo(src.Num_Haplo);
	memmove(List, src.List, sizeof(THaplotype)*src.Num_Haplo);
	return *this;
}

void CHaplotypeList::ResizeHaplo(size_t num)
{
	if (Num_Haplo != num)
	{
		Num_Haplo = num;
		if (num > reserve_size)
			alloc_mem(reserve_size = num);
	}
}

void CHaplotypeList::DoubleHaplos(CHaplotypeList &OutHaplos) const
{
	HIBAG_CHECKING(Num_SNP >= HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"CHaplotypeList::DoubleHaplos, there are too many SNP markers.");

	OutHaplos.Num_SNP = Num_SNP + 1;
	OutHaplos.ResizeHaplo(Num_Haplo*2);
	const THaplotype *pSrc = List;
	THaplotype *pDst = OutHaplos.List;

	// double haplotypes
	for (size_t i=0; i < Num_Haplo; i++)
	{
		*pDst = *pSrc;
		pDst->_SetAllele(Num_SNP, 0); pDst ++;
		*pDst = *pSrc;
		pDst->_SetAllele(Num_SNP, 1); pDst ++;
		pSrc ++;
	}

	// double the count of haplotypes for each HLA allele
	size_t nhla = nHLA();
	OutHaplos.LenPerHLA.resize(nhla);
	const size_t *s = &LenPerHLA[0];
	size_t *p = &OutHaplos.LenPerHLA[0];
	for (; nhla > 0; nhla--) *p++ = (*s++) * 2;
}

void CHaplotypeList::DoubleHaplosInitFreq(CHaplotypeList &OutHaplos,
	double AFreq) const
{
	HIBAG_CHECKING(Num_Haplo*2 != OutHaplos.Num_Haplo,
		"CHaplotypeList::DoubleHaplosInitFreq, the total number of haplotypes is not correct.");

	const double p0 = 1-AFreq, p1 = AFreq;
	const THaplotype *s = List;
	THaplotype *p = OutHaplos.List;
	for (size_t n=Num_Haplo; n > 0; n--)
	{
		p[0].Freq = p0 * s->Freq + EM_INIT_VAL_FRAC;
		p[1].Freq = p1 * s->Freq + EM_INIT_VAL_FRAC;
		p += 2; s ++;
	}
}

void CHaplotypeList::EraseDoubleHaplos(double RareProb, CHaplotypeList &OutHaplos) const
{
	// count the number of haplotypes after filtering rare haplotypes
	size_t num = 0;
	const THaplotype *p = List;
	for (size_t i=0; i < Num_Haplo; i+=2, p+=2)
	{
		if ((p[0].Freq < RareProb) || (p[1].Freq < RareProb))
		{
			if ((p[0].Freq + p[1].Freq) >= MIN_RARE_FREQ)
				num ++;
		} else {
			num += 2;
		}
	}

	// initialize the output
	OutHaplos.Num_SNP = Num_SNP;
	OutHaplos.ResizeHaplo(num);
	OutHaplos.LenPerHLA.resize(LenPerHLA.size());

	// assign haplotypes
	p = List;
	THaplotype *pOut = OutHaplos.List;
	double sum = 0;
	for (size_t i=0; i < LenPerHLA.size(); i++)
	{
		num = 0;
		for (size_t n=LenPerHLA[i]; n > 0; n-=2, p+=2)
		{
			double sumfreq = p[0].Freq + p[1].Freq;
			if ((p[0].Freq < RareProb) || (p[1].Freq < RareProb))
			{
				if (sumfreq >= MIN_RARE_FREQ)
				{
					if (p[0].Freq >= p[1].Freq)
						*pOut = p[0];
					else
						*pOut = p[1];
					pOut->Freq = sumfreq; pOut ++;
					sum += sumfreq;
					num ++;
				}
			} else {
				*pOut++ = p[0];
				*pOut++ = p[1];
				sum += sumfreq;
				num += 2;
			}
		}
		OutHaplos.LenPerHLA[i] = num;
	}

	OutHaplos.ScaleFrequency(1/sum);
}

void CHaplotypeList::SaveClearFrequency()
{
	THaplotype *p = List;
	for (size_t n=Num_Haplo; n > 0; n--, p++)
	{
		p->aux.OldFreq = p->Freq;
		p->Freq = 0;
	}
}

void CHaplotypeList::ScaleFrequency(double scale)
{
	THaplotype *p = List;
	for (size_t n=Num_Haplo; n > 0; n--, p++)
		p->Freq *= scale;
}

size_t CHaplotypeList::StartHaploHLA(int hla) const
{
	HIBAG_CHECKING(hla < 0 || hla >= (int)LenPerHLA.size(),
		"CHaplotypeList::StartHLA, invalid HLA allele.");
	size_t rv = 0;
	for (int i=0; i < hla; i++) rv += LenPerHLA[i];
	return rv;
}

void CHaplotypeList::SetHaploAux(int64_t buf_haplo[], double buf_freq[])
{
	aux_haplo = buf_haplo;
	aux_freq = buf_freq;
	if (Num_SNP <= 64)
	{
		for (size_t i=0; i < Num_Haplo; i++)
		{
			buf_haplo[i] = List[i].PackedHaplo[0];
			buf_freq[i] = List[i].Freq;
		}
	} else {
		int64_t *buf_haplo2 = buf_haplo + Num_Haplo;
		for (size_t i=0; i < Num_Haplo; i++)
		{
			buf_haplo[i] = List[i].PackedHaplo[0];
			buf_haplo2[i] = List[i].PackedHaplo[1];
			buf_freq[i] = List[i].Freq;
		}
	}
}

void CHaplotypeList::SetHaploAux_GPU()
{
	THaplotype *p = List;
	size_t *s = &LenPerHLA[0];
	size_t n = LenPerHLA.size();
	for (size_t i=0; i < n; i++)
	{
		for (size_t m=*s++; m > 0; m--, p++)
		{
			p->aux.a2.Freq_f32 = p->Freq;
			p->aux.a2.HLA_allele = i;
		}
	}
}



// ========================================================================= //
// Packed bi-allelic SNP genotype structure: 8 SNPs in a byte

TGenotype::TGenotype()
{
	BootstrapCount = 0;
}

int TGenotype::GetSNP(size_t idx) const
{
	HIBAG_CHECKING(idx >= HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"TGenotype::GetSNP, invalid index.");
	const size_t i = idx >> 3, r = idx & 0x07;
	const UINT8 *S1 = (const UINT8*)PackedSNP1;
	const UINT8 *S2 = (const UINT8*)PackedSNP2;
	const int b1 = (S1[i] >> r) & 0x01;
	const int b2 = (S2[i] >> r) & 0x01;
	return (b1==0 && b2==1) ? -1 : b1+b2;
}

void TGenotype::SetSNP(size_t idx, int val)
{
	HIBAG_CHECKING(idx >= HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"TGenotype::SetSNP, invalid index.");
	_SetSNP(idx, val);
}

void TGenotype::_SetSNP(size_t idx, int val)
{
	const size_t i = idx >> 3, r = idx & 0x07;
	const UINT8 SET = (UINT8(0x01) << r), CLEAR = ~SET;
	UINT8 &S1 = ((UINT8*)PackedSNP1)[i];
	UINT8 &S2 = ((UINT8*)PackedSNP2)[i];
	switch (val)
	{
		case 0:  S1 &= CLEAR; S2 &= CLEAR; break;
		case 1:  S1 |= SET;   S2 &= CLEAR; break;
		case 2:  S1 |= SET;   S2 |= SET;   break;
		default: S1 &= CLEAR; S2 |= SET;
	}
}

string TGenotype::SNPToString(size_t Length) const
{
	HIBAG_CHECKING(Length > HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"TGenotype::SNPToString, the length is too large.");
	string rv;
	if (Length > 0)
	{
		rv.resize(Length);
		for (size_t i=0; i < Length; i++)
		{
			UINT8 ch = GetSNP(i);
			rv[i] = (ch < 3) ? (ch + '0') : '?';
		}
	}
	return rv;
}

void TGenotype::StringToSNP(const string &str)
{
	HIBAG_CHECKING(str.size() > HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"TGenotype::StringToSNP, the input string is too long.");
	for (size_t i=0; i < str.size(); i++)
	{
		const char ch = str[i];
		HIBAG_CHECKING(ch!='0' && ch!='1' && ch!='2' && ch!='?',
			"TGenotype::StringToSNP, the input string should be '0', '1', '2' or '?'.");
		_SetSNP(i, ch-'0');
	}
}

void TGenotype::SNPToInt(size_t Length, int OutArray[]) const
{
	HIBAG_CHECKING(Length > HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"TGenotype::SNPToInt, the length is invalid.");
	for (size_t i=0; i < Length; i++)
		OutArray[i] = GetSNP(i);
}

static ALWAYS_INLINE size_t geno4(int g)
{
	return (0<=g && g<=2) ? g : 3;
}

void TGenotype::IntToSNP(size_t Length, const int GenoBase[], const int Index[])
{
	HIBAG_CHECKING(Length > HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"TGenotype::IntToSNP, the length is invalid.");
	// fill all bytes with missing value
	memset(PackedSNP1,  0, sizeof(PackedSNP1));
	memset(PackedSNP2, -1, sizeof(PackedSNP2));
	// set genotypes
	const static UINT8 P1[4] = { 0, 1, 1, 0 };
	const static UINT8 P2[4] = { 0, 0, 1, 1 };
	UINT8 *p1 = (UINT8*)PackedSNP1;
	UINT8 *p2 = (UINT8*)PackedSNP2;
	// packed fill
	for (; Length >= 8; Length -= 8, Index += 8)
	{
		const size_t i0 = geno4(GenoBase[Index[0]]);
		const size_t i1 = geno4(GenoBase[Index[1]]);
		const size_t i2 = geno4(GenoBase[Index[2]]);
		const size_t i3 = geno4(GenoBase[Index[3]]);
		const size_t i4 = geno4(GenoBase[Index[4]]);
		const size_t i5 = geno4(GenoBase[Index[5]]);
		const size_t i6 = geno4(GenoBase[Index[6]]);
		const size_t i7 = geno4(GenoBase[Index[7]]);
		*p1++ = P1[i0] | (P1[i1] << 1) | (P1[i2] << 2) | (P1[i3] << 3) |
			(P1[i4] << 4) | (P1[i5] << 5) | (P1[i6] << 6) | (P1[i7] << 7);
		*p2++ = P2[i0] | (P2[i1] << 1) | (P2[i2] << 2) | (P2[i3] << 3) |
			(P2[i4] << 4) | (P2[i5] << 5) | (P2[i6] << 6) | (P2[i7] << 7);
	}
	if (Length > 0)
	{
		*p1 = 0; *p2 = 0xFF;
		for (size_t i=0; i < Length; i++)
		{
			const size_t ii = geno4(GenoBase[*Index++]);
			*p1 = (*p1) | (P1[ii] << i);
			*p2 = ((*p2) & ~(1 << i)) | (P2[ii] << i);
		}
		p1++; p2++;
	}
}


// --------------------------------
// Compute the Hamming distance between SNPs and H1+H2 without checking

#ifdef HIBAG_CPU_ARCH_X86
	typedef int64_t UTYPE;
	#define U_POPCOUNT     __builtin_popcountll
#else
	#ifdef HIBAG_CPU_LP64
		typedef uint64_t UTYPE;
	#else
		typedef uint32_t UTYPE;
	#endif
#endif

#if !(defined(HIBAG_CPU_ARCH_X86) && defined(__SSE2__))
static const ssize_t UTYPE_BIT_NUM = sizeof(UTYPE)*8;
#endif

#ifndef U_POPCOUNT
#   define U_POPCOUNT  u_popcount
#   ifdef HIBAG_CPU_LP64
	static ALWAYS_INLINE int u_popcount(uint64_t v)
	{
		v -= ((v >> 1) & 0x5555555555555555LLU);
		v = (v & 0x3333333333333333LLU) + ((v >> 2) & 0x3333333333333333LLU);
		return (((v + (v >> 4)) & 0x0F0F0F0F0F0F0F0FLLU) *
			0x0101010101010101LLU) >> 56;
	}
#   else
	static ALWAYS_INLINE int u_popcount(uint32_t v)
	{
		v -= ((v >> 1) & 0x55555555);
		v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
		return (((v + (v >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
	}
#   endif
#endif

static ALWAYS_INLINE int hamm_d(size_t Length, const TGenotype &G,
	const THaplotype &H1, const THaplotype &H2)
{
	const UTYPE *h1 = (const UTYPE*)H1.PackedHaplo;
	const UTYPE *h2 = (const UTYPE*)H2.PackedHaplo;
	const UTYPE *s1 = (const UTYPE*)G.PackedSNP1;
	const UTYPE *s2 = (const UTYPE*)G.PackedSNP2;
#if defined(HIBAG_CPU_ARCH_X86) && defined(__SSE2__)
	// here, UTYPE = int64_t
	if (Length <= 64)
	{
		__m128i H  = { h1[0], h2[0] };  // two haplotypes
		__m128i S1 = { s1[0], s2[0] }, S2 = { s2[0], s1[0] };  // genotypes
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
		__m128i S1 = { s1[0], s1[1] }, S2 = { s2[0], s2[1] };  // genotypes
		// worry about n < UTYPE_BIT_NUM? unused bits are set to be a missing flag
		__m128i M = _mm_andnot_si128(S1, S2);  // missing value, 1 is missing
		__m128i MASK = _mm_andnot_si128(M, (H1 ^ S2) | (H2 ^ S1));
		__m128i va = (H1 ^ S1) & MASK;  // (H1 ^ S1) & MASK, (H2 ^ S2) & MASK
		__m128i vb = (H2 ^ S2) & MASK;
		// popcount
		return U_POPCOUNT(va[0]) + U_POPCOUNT(va[1]) +
			U_POPCOUNT(vb[0]) + U_POPCOUNT(vb[1]);
	}
#else
	size_t ans = 0;
	for (ssize_t n=Length; n > 0; n -= UTYPE_BIT_NUM)
	{
		UTYPE H1 = *h1++, H2 = *h2++;  // two haplotypes
		UTYPE S1 = *s1++, S2 = *s2++;  // genotypes
		// worry about n < UTYPE_BIT_NUM? unused bits are set to be a missing flag
		UTYPE M = S2 & ~S1;  // missing value, 1 is missing
		UTYPE MASK = ((H1 ^ S2) | (H2 ^ S1)) & ~M;
		// popcount
		ans += U_POPCOUNT((H1 ^ S1) & MASK) + U_POPCOUNT((H2 ^ S2) & MASK);
	}
	return ans;
#endif
}

// --------------------------------

int TGenotype::HammingDistance(size_t Length,
	const THaplotype &H1, const THaplotype &H2) const
{
	HIBAG_CHECKING(Length > HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"THaplotype::HammingDistance, the length is too large.");
	return hamm_d(Length, *this, H1, H2);
}


// -------------------------------------------------------------------------
// The class of SNP genotype list

CSNPGenoMatrix::CSNPGenoMatrix()
{
	Num_Total_SNP = Num_Total_Samp = 0;
	pGeno = NULL;
}

int CSNPGenoMatrix::Get(const int IdxSamp, const int IdxSNP) const
{
	return pGeno[IdxSamp*Num_Total_SNP + IdxSNP];
}

int *CSNPGenoMatrix::Get(const int IdxSamp)
{
	return pGeno + IdxSamp * Num_Total_SNP;
}


// -------------------------------------------------------------------------
// The class of SNP genotype list

CGenotypeList::CGenotypeList()
{
	Num_SNP = 0;
}

void CGenotypeList::AddSNP(int IdxSNP, const CSNPGenoMatrix &SNPMat)
{
	HIBAG_CHECKING(nSamp() != SNPMat.Num_Total_Samp,
		"CGenotypeList::AddSNP, SNPMat should have the same number of samples.");
	HIBAG_CHECKING(Num_SNP >= HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"CGenotypeList::AddSNP, there are too many SNP markers.");
	const int *pG = SNPMat.pGeno + IdxSNP;
	for (int i=0; i < SNPMat.Num_Total_Samp; i++)
	{
		int g = *pG;
		pG += SNPMat.Num_Total_SNP;
		List[i]._SetSNP(Num_SNP, g);
	}
	Num_SNP ++;
}

void CGenotypeList::ReduceSNP()
{
	HIBAG_CHECKING(Num_SNP <= 0,
		"CGenotypeList::ReduceSNP, there is no SNP marker.");
	Num_SNP --;
}

void CGenotypeList::SetAllMissing()
{
	size_t n = List.size();
	for (TGenotype *p = &List[0]; n > 0; n--, p++)
	{
		memset(p->PackedSNP1,  0, sizeof(p->PackedSNP1));
		memset(p->PackedSNP2, -1, sizeof(p->PackedSNP2));
	}
}

void CGenotypeList::SetMissing(int idx)
{
	const size_t i = idx >> 3, r = idx & 0x07;
	const UINT8 SET = (UINT8(0x01) << r), CLEAR = ~SET;
	size_t n = List.size();
	for (TGenotype *p = &List[0]; n > 0; n--, p++)
	{
		((UINT8*)p->PackedSNP1)[i] &= CLEAR;
		((UINT8*)p->PackedSNP2)[i] |= SET;
	}
}



// -------------------------------------------------------------------------
// A list of HLA types

CHLATypeList::CHLATypeList() { }

inline int CHLATypeList::Compare(const THLAType &H1, const THLAType &H2)
{
	int P1=H1.Allele1, P2=H1.Allele2;
	int T1=H2.Allele1, T2=H2.Allele2;
	int cnt = 0;
	if ((P1==T1) || (P1==T2))
	{
		cnt = 1;
		if (P1==T1) T1 = -1; else T2 = -1;
	}
	if ((P2==T1) || (P2==T2)) cnt++;
	return cnt;
}


// -------------------------------------------------------------------------
// CSamplingWithoutReplace

CSamplingWithoutReplace::CSamplingWithoutReplace()
{
	_m_try = 0;
}

CBaseSampling *CSamplingWithoutReplace::Init(int m_total)
{
	_m_try = 0;
	_IdxArray.resize(m_total);
	for (int i=0; i < m_total; i++)
		_IdxArray[i] = i;
	return this;
}

int CSamplingWithoutReplace::TotalNum() const
{
	return _IdxArray.size();
}

void CSamplingWithoutReplace::RandomSelect(int m_try)
{
	const int n_tmp = _IdxArray.size();
	if (m_try > n_tmp) m_try = n_tmp;
	if (m_try < n_tmp)
	{
		for (int i=0; i < m_try; i++)
		{
			int I = RandomNum(n_tmp - i);
			std::swap(_IdxArray[I], _IdxArray[n_tmp-i-1]);
		}
	}
	_m_try = m_try;
}

int CSamplingWithoutReplace::NumOfSelection() const
{
	return _m_try;
}

void CSamplingWithoutReplace::Remove(int idx)
{
	idx = _IdxArray.size() - _m_try + idx;
	_IdxArray.erase(_IdxArray.begin() + idx);
}

void CSamplingWithoutReplace::RemoveSelection()
{
	_IdxArray.resize(_IdxArray.size() - _m_try);
}

void CSamplingWithoutReplace::RemoveFlag()
{
	const int n_tmp = _IdxArray.size();
	for (int i=n_tmp-1; i >= n_tmp - _m_try; i--)
	{
		vector<int>::iterator p = _IdxArray.begin() + i;
		if (*p < 0) _IdxArray.erase(p);
	}
}

int &CSamplingWithoutReplace::operator[] (int idx)
{
	return _IdxArray[_IdxArray.size() - _m_try + idx];
}



// -------------------------------------------------------------------------
// The class of SNP genotype list

CAlg_EM::CAlg_EM(CVariableSelection &v): vs(v) { }

void CAlg_EM::PrepareHaplotypes(const CHaplotypeList &CurHaplo,
	const CGenotypeList &GenoList, CHaplotypeList &NextHaplo)
{
	// create a next set of haplotypes
	CurHaplo.DoubleHaplos(NextHaplo);
	const int nhla = NextHaplo.nHLA();  // # of unique HLA alleles, not too large

	// the number of in-bag samples
	const size_t n_samp_ib = vs.idx_inbag.size();

	if (GPUExtProcPtr && *GPUExtProcPtr->build_haplomatch)
	{
		// prepare data and check
		int StartIdx[nhla], st=0;
		for (int i=0; i < nhla; i++)
		{
			if (CurHaplo.LenPerHLA[i] > 65535)
				throw "There are too many HLA allele-specific haplotypes (# > 65535).";
			StartIdx[i] = st;
			st += CurHaplo.LenPerHLA[i];
		}

		// initialize pair lists
		_SampHaploPair.clear();
		_SampHaploPair.resize(n_samp_ib);
		for (size_t i=0; i < n_samp_ib; i++)
		{
			size_t k = vs.idx_inbag[i];
			THaploPairList &HP = _SampHaploPair[i];
			HP.BootstrapCount = GenoList.List[k].BootstrapCount;
			HP.SampIndex = k;
		}

		// call GPU (have called 'build_set_bootstrap')
		size_t n_buf=0;
		UINT32 *buf = (*GPUExtProcPtr->build_haplomatch)(CurHaplo.List,
			&CurHaplo.LenPerHLA[0], CurHaplo.Num_SNP, &GenoList.List[0], n_buf);
		// set pair lists from the GPU output
		if (buf)
		{
			UINT32 n = *buf >> 1;
			for (UINT32 *p=buf+1; n > 0; n--, p+=2)
			{
				const size_t k = p[0];  // in-bag sample index
				std::vector<CAlg_EM::THaploPair> &PL = _SampHaploPair[k].PairList;
				const TGenotype &pG = GenoList.List[vs.idx_inbag[k]];
				const size_t pH1_st = StartIdx[pG.aux_hla_type.Allele1];
				const size_t pH2_st = StartIdx[pG.aux_hla_type.Allele2];
				const UINT32 v = p[1];
				THaplotype *p1 = &NextHaplo.List[2*(pH1_st + (v & 0xFFFF))];
				THaplotype *p2 = &NextHaplo.List[2*(pH2_st + (v >> 16))];
				// always p1 <= p2
				PL.push_back(CAlg_EM::THaploPair(p1, p2));
				PL.push_back(CAlg_EM::THaploPair(p1, p2+1));
				if (p1+1 <= p2)
					PL.push_back(CAlg_EM::THaploPair(p1+1, p2));
				PL.push_back(CAlg_EM::THaploPair(p1+1, p2+1));
			}
			// free the buffer
			free(buf);
		}

		// check in case GPU does not work appropriately
		for (size_t i=0; i < n_samp_ib; i++)
		{
			vector<THaploPair> &PL = _SampHaploPair[i].PairList;
			// std::sort(PL.begin(), PL.end(), THaploPair::Less);  // enable when debugging
			if (PL.empty())
				throw "PairList should not be empty in PrepareHaplotypes().";
		}

	} else {

		// prepare data
		int StartIdx[nhla], st=0;
		for (int i=0; i < nhla; i++)
		{
			StartIdx[i] = st;
			st += NextHaplo.LenPerHLA[i];
		}

		size_t n_max_diff = 0;
		vector<int>::const_iterator pEnd = vs.idx_inbag.end();
		for (vector<int>::const_iterator p=vs.idx_inbag.begin(); p != pEnd; p++)
		{
			const THLAType &pHLA = GenoList.List[*p].aux_hla_type;
			size_t m, pH1_n = NextHaplo.LenPerHLA[pHLA.Allele1];
			if (pHLA.Allele1 != pHLA.Allele2)
			{
				size_t pH2_n = NextHaplo.LenPerHLA[pHLA.Allele2];
				m = pH1_n * pH2_n;
			} else {
				m = pH1_n * (pH1_n + 1) / 2;
			}
			if (m > n_max_diff) n_max_diff = m;
		}

		_SampHaploPair.clear();
		_SampHaploPair.resize(n_samp_ib);
		vector<short> DiffList(n_max_diff*thread_num());

		// get haplotype pairs for each in-bag sample
		PARALLEL_FOR(i, n_samp_ib)
		{
			size_t k = vs.idx_inbag[i];
			const TGenotype &pG = GenoList.List[k];
			THaploPairList &HP = _SampHaploPair[i];
			HP.BootstrapCount = pG.BootstrapCount;
			HP.SampIndex = k;

			const THLAType &pHLA = pG.aux_hla_type;
			const size_t pH1_st  = StartIdx[pHLA.Allele1];
			const size_t pH1_n   = NextHaplo.LenPerHLA[pHLA.Allele1];
			const size_t pH2_st  = StartIdx[pHLA.Allele2];
			const size_t pH2_n   = NextHaplo.LenPerHLA[pHLA.Allele2];
			const size_t Num_SNP = NextHaplo.Num_SNP - 1;

			(*fc_PrepHaploMatch)(pG, &NextHaplo.List[pH1_st], pH1_n,
				&NextHaplo.List[pH2_st], pH2_n, Num_SNP, HP.PairList,
				&DiffList[n_max_diff*thread_idx()]);
		}
		PARALLEL_END
	}
}

bool CAlg_EM::PrepareNewSNP(const int NewSNP, const CHaplotypeList &CurHaplo,
	const CSNPGenoMatrix &SNPMat, CGenotypeList &GenoList,
	CHaplotypeList &NextHaplo)
{
	HIBAG_CHECKING((NewSNP<0) || (NewSNP>=SNPMat.Num_Total_SNP),
		"CAlg_EM::PrepareNewSNP, invalid NewSNP.");
	HIBAG_CHECKING(SNPMat.Num_Total_Samp != GenoList.nSamp(),
		"CAlg_EM::PrepareNewSNP, SNPMat and GenoList should have the same number of SNPs.");

	// compute the allele frequency of NewSNP
	int allele_cnt=0, valid_cnt=0;
	vector<int>::const_iterator pEnd = vs.idx_inbag.end();
	for (vector<int>::const_iterator p=vs.idx_inbag.begin(); p != pEnd; p++)
	{
		const size_t i = *p;
		int dup = GenoList.List[i].BootstrapCount;
		int g = SNPMat.Get(i, NewSNP);
		if ((0<=g) && (g<=2))
		{
			allele_cnt += g*dup;
			valid_cnt  += 2*dup;
		}
	}
	if ((allele_cnt==0) || (allele_cnt==valid_cnt))
		return false;

	// initialize the haplotype frequencies in NextHaplo
	CurHaplo.DoubleHaplosInitFreq(NextHaplo, double(allele_cnt)/valid_cnt);

	// update haplotype pair
	const size_t IdxNewSNP = NextHaplo.Num_SNP - 1;
	const size_t num = _SampHaploPair.size();
	PARALLEL_FOR_CHUNK(i, num, 16)
	{
		THaploPairList &PL = _SampHaploPair[i];
		// SNP genotype
		int geno = SNPMat.Get(PL.SampIndex, NewSNP);
		if ((0<=geno) && (geno<=2))
		{
			// check for each pair
			vector<THaploPair>::iterator p;
			for (p = PL.PairList.begin(); p != PL.PairList.end(); p++)
			{
				p->Flag = ((p->H1->GetAllele(IdxNewSNP) +
					p->H2->GetAllele(IdxNewSNP)) == geno);
			}
		} else {
			// missing genotype
			vector<THaploPair>::iterator p;
			for (p = PL.PairList.begin(); p != PL.PairList.end(); p++)
				p->Flag = true;
		}
	}
	PARALLEL_END

	return true;
}

void CAlg_EM::ExpectationMaximization(CHaplotypeList &NextHaplo)
{
	const int TotalNumSamp = vs.nSamp();  // sum of BootstrapCount = total # of samples
	double ConvTol = 0, LogLik = -1e+30;  // the converage tolerance

	if (log_buf.size() < (size_t)TotalNumSamp)
		log_buf.resize(TotalNumSamp);

	// iterate ...
	for (int iter=0; iter <= EM_MaxNum_Iterations; iter++)
	{
		// save old values
		// old log likelihood
		double Old_LogLik = LogLik;
		// save old haplotype frequencies
		NextHaplo.SaveClearFrequency();

		// for-loop each in-bag sample
		const size_t num = _SampHaploPair.size();
		PARALLEL_FOR_CHUNK(i, num, 16)
		{
			THaploPairList &PL = _SampHaploPair[i];
			double psum = 0;
			vector<THaploPair>::iterator p;
			for (p = PL.PairList.begin(); p != PL.PairList.end(); p++)
			{
				if (p->Flag)
				{
					p->GenoFreq = (p->H1 != p->H2) ?
						(2 * p->H1->aux.OldFreq * p->H2->aux.OldFreq) :
						(p->H1->aux.OldFreq * p->H2->aux.OldFreq);
					psum += p->GenoFreq;
				}
			}
			// always "PL.BootstrapCount > 0"
			log_buf[i] = PL.BootstrapCount * log(psum);
			psum = PL.BootstrapCount / psum;
			// update
			for (p = PL.PairList.begin(); p != PL.PairList.end(); p++)
				if (p->Flag) p->GenoFreq *= psum;
		}
		PARALLEL_END

		// update
		LogLik = 0;
		for (size_t i=0; i < num; i++)
		{
			LogLik += log_buf[i];
			vector<THaploPair> &PL = _SampHaploPair[i].PairList;
			for (vector<THaploPair>::iterator p=PL.begin(); p != PL.end(); p++)
			{
				if (p->Flag)
				{
					double r = p->GenoFreq;
					p->H1->Freq += r; p->H2->Freq += r;
				}
			}
		}

		// finally
		NextHaplo.ScaleFrequency(0.5/TotalNumSamp);
		if (iter > 0)
		{
			if (fabs(LogLik - Old_LogLik) <= ConvTol)
				break;
		} else {
			ConvTol = EM_FuncRelTol * (fabs(LogLik) + EM_FuncRelTol);
			if (ConvTol < 0) ConvTol = 0;
		}
	}
}



// -------------------------------------------------------------------------
// The algorithm of prediction

static bool need_auxiliary_haplo = false;

CAlg_Prediction::CAlg_Prediction() { }

void CAlg_Prediction::Init_Target_IFunc(const char *cpu)
{
	bool cpu_auto, cpu_max;
	if (!cpu) cpu = "";
	if (strcmp(cpu, "auto.avx2")==0) cpu = "";
	cpu_auto = (strlen(cpu)==0);
	cpu_max = (strcmp(cpu, "max")==0);
	if (cpu_max) cpu_auto = true;
	string cpu_info;
	bool need_aux_haplo = false;

#ifdef HIBAG_CPU_LP64
	cpu_info = "64-bit";
#else
	cpu_info = "32-bit";
#endif

#ifdef HIBAG_CPU_ARCH_X86

#ifdef HIBAG_BUILTIN_CPU
	__builtin_cpu_init();
#endif

#if defined(__AVX512F__) && defined(__AVX512BW__)
	const bool has_avx512bw = true;
#elif defined(HIBAG_BUILTIN_CPU_AVX512BW)
#if defined(__ICC) && defined(_FEATURE_AVX512F) && defined(_FEATURE_AVX512BW)
	// since __builtin_cpu_supports("avx512bw") always return 0
	const bool has_avx512bw = _may_i_use_cpu_feature(
		_FEATURE_AVX512F | _FEATURE_AVX512BW) && HIBAG_ALGORITHM_AVX512BW;
#else
	const bool has_avx512bw = __builtin_cpu_supports("avx512f") &&
		__builtin_cpu_supports("avx512bw") && HIBAG_ALGORITHM_AVX512BW;
#endif
#else
	const bool has_avx512bw = false;
#endif

#if defined(__AVX512F__)
	const bool has_avx512f = true;
#elif defined(HIBAG_BUILTIN_CPU_AVX512F)
	const bool has_avx512f = __builtin_cpu_supports("avx512f") &&
		HIBAG_ALGORITHM_AVX512F;
#else
	const bool has_avx512f = false;
#endif

#if defined(__AVX2__)
	const bool has_avx2 = true;
#elif defined(HIBAG_BUILTIN_CPU_AVX2)
	const bool has_avx2 = __builtin_cpu_supports("avx2") && HIBAG_ALGORITHM_AVX2;
#else
	const bool has_avx2 = false;
#endif

#if defined(__AVX__)
	const bool has_avx = true;
#elif defined(HIBAG_BUILTIN_CPU_AVX)
	const bool has_avx  = __builtin_cpu_supports("avx") && HIBAG_ALGORITHM_AVX;
#else
	const bool has_avx = false;
#endif

#if defined(__SSE4_2__)
	const bool has_sse4 = true;
#elif defined(HIBAG_BUILTIN_CPU) && defined(HIBAG_CPU_ARCH_X86_SSE4_2)
	const bool has_sse4 = __builtin_cpu_supports("sse4.2") &&
		__builtin_cpu_supports("popcnt") && HIBAG_ALGORITHM_SSE4_2;
#else
	const bool has_sse4 = false;
#endif

#if defined(__SSE2__)
	const bool has_sse2 = true;
#elif defined(HIBAG_BUILTIN_CPU_SSE2)
	const bool has_sse2 = __builtin_cpu_supports("sse2") && HIBAG_ALGORITHM_SSE2;
#else
	const bool has_sse2 = false;
#endif

	bool flag_popcnt = false;
	if (cpu_max && has_avx512bw)
	{
		if (!has_avx512bw)
			error("Not support AVX512F+AVX512BW.");
		fc_PrepHaploMatch = &CAlg_Prediction::_PrepHaploMatch_avx512bw;
		fc_BestGuess = &CAlg_Prediction::_BestGuess_avx512bw;
		fc_PostProb  = &CAlg_Prediction::_PostProb_avx512bw;
		fc_PostProb2 = &CAlg_Prediction::_PostProb2_avx512bw;
		cpu_info.append(", AVX512F+AVX512BW");
		need_aux_haplo = true;
	} else if (strcmp(cpu, "base")==0)
	{
		fc_PrepHaploMatch = &CAlg_Prediction::_PrepHaploMatch_def;
		fc_BestGuess = &CAlg_Prediction::_BestGuess_def;
		fc_PostProb  = &CAlg_Prediction::_PostProb_def;
		fc_PostProb2 = &CAlg_Prediction::_PostProb2_def;
	#ifdef __POPCNT__
		flag_popcnt = true;
	#endif
	} else if (strcmp(cpu, "avx2")==0 || (cpu_auto && has_avx2))
	{
		if (!has_avx2)
			error("Not support AVX2.");
		fc_PrepHaploMatch = &CAlg_Prediction::_PrepHaploMatch_avx2;
		fc_BestGuess = &CAlg_Prediction::_BestGuess_avx2;
		fc_PostProb  = &CAlg_Prediction::_PostProb_avx2;
		fc_PostProb2 = &CAlg_Prediction::_PostProb2_avx2;
		cpu_info.append(", AVX2");
		need_aux_haplo = true;
	} else if (strcmp(cpu, "avx")==0 || (cpu_auto && has_avx))
	{
		if (!has_avx)
			error("Not support AVX.");
		fc_PrepHaploMatch = &CAlg_Prediction::_PrepHaploMatch_avx;
		fc_BestGuess = &CAlg_Prediction::_BestGuess_avx;
		fc_PostProb  = &CAlg_Prediction::_PostProb_avx;
		fc_PostProb2 = &CAlg_Prediction::_PostProb2_avx;
		cpu_info.append(", AVX");
		need_aux_haplo = true;
	} else if (strcmp(cpu, "sse4")==0 || (cpu_auto && has_sse4))
	{
		if (!has_sse4)
			error("Not support SSE4.2.");
		fc_PrepHaploMatch = &CAlg_Prediction::_PrepHaploMatch_sse4_2;
		fc_BestGuess = &CAlg_Prediction::_BestGuess_sse4_2;
		fc_PostProb  = &CAlg_Prediction::_PostProb_sse4_2;
		fc_PostProb2 = &CAlg_Prediction::_PostProb2_sse4_2;
		cpu_info.append(", SSE4.2");
		flag_popcnt = true;
	} else if (strcmp(cpu, "sse2")==0 || (cpu_auto && has_sse2))
	{
		if (!has_sse2)
			error("Not support SSE2.");
		fc_PrepHaploMatch = &CAlg_Prediction::_PrepHaploMatch_sse2;
		fc_BestGuess = &CAlg_Prediction::_BestGuess_sse2;
		fc_PostProb  = &CAlg_Prediction::_PostProb_sse2;
		fc_PostProb2 = &CAlg_Prediction::_PostProb2_sse2;
		cpu_info.append(", SSE2");
		flag_popcnt = HIBAG_ALGORITHM_SSE2_POPCNT;
	} else if (strcmp(cpu, "avx512bw")==0 || (cpu_auto && has_avx512bw))
	{
		if (!has_avx512bw)
			error("Not support AVX512F+AVX512BW.");
		fc_PrepHaploMatch = &CAlg_Prediction::_PrepHaploMatch_avx512bw;
		fc_BestGuess = &CAlg_Prediction::_BestGuess_avx512bw;
		fc_PostProb  = &CAlg_Prediction::_PostProb_avx512bw;
		fc_PostProb2 = &CAlg_Prediction::_PostProb2_avx512bw;
		cpu_info.append(", AVX512F+AVX512BW");
		need_aux_haplo = true;
	} else if (strcmp(cpu, "avx512f")==0 || (cpu_auto && has_avx512f))
	{
		if (!has_avx512f)
			error("Not support AVX512F.");
		fc_PrepHaploMatch = &CAlg_Prediction::_PrepHaploMatch_avx512f;
		fc_BestGuess = &CAlg_Prediction::_BestGuess_avx512f;
		fc_PostProb  = &CAlg_Prediction::_PostProb_avx512f;
		fc_PostProb2 = &CAlg_Prediction::_PostProb2_avx512f;
		cpu_info.append(", AVX512F");
		need_aux_haplo = true;
	} else {
		fc_PrepHaploMatch = &CAlg_Prediction::_PrepHaploMatch_def;
		fc_BestGuess = &CAlg_Prediction::_BestGuess_def;
		fc_PostProb  = &CAlg_Prediction::_PostProb_def;
		fc_PostProb2 = &CAlg_Prediction::_PostProb2_def;
	}
	if (flag_popcnt)
		cpu_info.append(", POPCNT");
#ifdef __FMA__
	cpu_info.append(", FMA");
#endif
#endif
	HIBAG_CPU_Info = cpu_info;
	need_auxiliary_haplo = need_aux_haplo;
}

void CAlg_Prediction::InitPrediction(int n_hla)
{
	HIBAG_CHECKING(n_hla<=0, "CAlg_Prediction::Init, n_hla error.");
	_nHLA = n_hla;
	const size_t size = n_hla*(n_hla+1)/2;
	_PostProb.resize(size);
	_SumPostProb.resize(size);
}

void CAlg_Prediction::InitPostProb()
{
	memset(&_PostProb[0], 0, _PostProb.size()*sizeof(double));
}

void CAlg_Prediction::InitSumPostProb()
{
	memset(&_SumPostProb[0], 0, _SumPostProb.size()*sizeof(double));
	_Sum_Weight = 0;
}

void CAlg_Prediction::AddProbToSum(double weight)
{
	if (weight > 0)
	{
		const size_t n = _SumPostProb.size();
		double *p = &_PostProb[0], *s = &_SumPostProb[0];
		for (size_t i=0; i < n; i++)
			s[i] += p[i] * weight;
		_Sum_Weight += weight;
	}
}

void CAlg_Prediction::NormalizeSumPostProb()
{
	if (_Sum_Weight > 0)
	{
		const double ff = 1.0 / _Sum_Weight;
		const size_t n = _SumPostProb.size();
		double *s = &_SumPostProb[0];
		for (size_t i=0; i < n; i++) s[i] *= ff;
	}
}

double &CAlg_Prediction::IndexPostProb(int H1, int H2)
{
	if (H1 > H2) std::swap(H1, H2);
	return _PostProb[H2 + H1*(2*_nHLA-H1-1)/2];
}

double &CAlg_Prediction::IndexSumPostProb(int H1, int H2)
{
	if (H1 > H2) std::swap(H1, H2);
	return _SumPostProb[H2 + H1*(2*_nHLA-H1-1)/2];
}

inline void CAlg_Prediction::PredictPostProb(const CHaplotypeList &Haplo,
	const TGenotype &Geno, double &SumProb)
{
	SumProb = (*fc_PostProb2)(Haplo, Geno, &_PostProb[0]);
}

THLAType CAlg_Prediction::BestGuess()
{
	return _BestGuess(&_PostProb[0]);
}

THLAType CAlg_Prediction::BestGuessEnsemble()
{
	return _BestGuess(&_SumPostProb[0]);
}

/// the best-guess HLA type from prob
inline THLAType CAlg_Prediction::_BestGuess(double hla_prob[])
{
	THLAType rv;
	rv.Allele1 = rv.Allele2 = NA_INTEGER;
	double *p=hla_prob, max = 0;
	for (int h1=0; h1 < _nHLA; h1++)
	{
		for (int h2=h1; h2 < _nHLA; h2++, p++)
		{
			if (max < *p)
			{
				max = *p;
				rv.Allele1 = h1; rv.Allele2 = h2;
			}
		}
	}
	return rv;
}


void CAlg_Prediction::_PrepHaploMatch_def(const TGenotype &Geno,
	THaplotype *pH1_st, size_t pH1_n, THaplotype *pH2_st, size_t pH2_n,
	size_t Num_SNP, std::vector<CAlg_EM::THaploPair> &HP_PairList, short DiffList[])
{
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
				int d = hamm_d(Num_SNP, Geno, *p1, *p2);
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
				int d = hamm_d(Num_SNP, Geno, *p1, *p2);
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

THLAType CAlg_Prediction::_BestGuess_def(const CHaplotypeList &Haplo,
	const TGenotype &Geno)
{
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
			ADD_FREQ_MUTANT(prob, i1->Freq * i1->Freq,
				hamm_d(Haplo.Num_SNP, Geno, *i1, *i1));
			// i2 > i1
			const double ff = 2 * i1->Freq;
			THaplotype *i2 = i1 + 1;
			for (size_t m2=m1-1; m2 > 0; m2--, i2++)
			{
				ADD_FREQ_MUTANT(prob, ff * i2->Freq,
					hamm_d(Haplo.Num_SNP, Geno, *i1, *i2));
			}
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
				{
					ADD_FREQ_MUTANT(prob, ff * i2->Freq,
						hamm_d(Haplo.Num_SNP, Geno, *i1, *i2));
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

double CAlg_Prediction::_PostProb_def(const CHaplotypeList &Haplo,
	const TGenotype &Geno, const THLAType &HLA)
{
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
			ADD_FREQ_MUTANT(prob, i1->Freq * i1->Freq,
				hamm_d(Haplo.Num_SNP, Geno, *i1, *i1));
			// i2 > i1
			const double ff = 2 * i1->Freq;
			THaplotype *i2 = i1 + 1;
			for (size_t m2=m1-1; m2 > 0; m2--, i2++)
			{
				ADD_FREQ_MUTANT(prob, ff * i2->Freq,
					hamm_d(Haplo.Num_SNP, Geno, *i1, *i2));
			}
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
				{
					ADD_FREQ_MUTANT(prob, ff * i2->Freq,
						hamm_d(Haplo.Num_SNP, Geno, *i1, *i2));
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

double CAlg_Prediction::_PostProb2_def(const CHaplotypeList &Haplo,
	const TGenotype &Geno, double Prob[])
{
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
			ADD_FREQ_MUTANT(sum, i1->Freq * i1->Freq,
				hamm_d(Haplo.Num_SNP, Geno, *i1, *i1));
			// i2 > i1
			const double ff = 2 * i1->Freq;
			THaplotype *i2 = i1 + 1;
			for (size_t m2=m1-1; m2 > 0; m2--, i2++)
			{
				ADD_FREQ_MUTANT(sum, ff * i2->Freq,
					hamm_d(Haplo.Num_SNP, Geno, *i1, *i2));
			}
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
				{
					ADD_FREQ_MUTANT(sum, ff * i2->Freq,
						hamm_d(Haplo.Num_SNP, Geno, *i1, *i2));
				}
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



// -------------------------------------------------------------------------
// The algorithm of variable selection

CVariableSelection::CVariableSelection(): _EM(*this)
{
	_SNPMat = NULL;
	_HLAList = NULL;
}

void CVariableSelection::InitSelection(CSNPGenoMatrix &snpMat,
	CHLATypeList &hlaList, const int _BootstrapCnt[])
{
	HIBAG_CHECKING(snpMat.Num_Total_Samp != hlaList.nSamp(),
		"CVariableSelection::InitSelection, snpMat and hlaList should have the same number of samples.");
	_SNPMat = &snpMat;
	_HLAList = &hlaList;
	// initialize genotype list
	_GenoList.List.resize(snpMat.Num_Total_Samp);
	log_inbag.resize(snpMat.Num_Total_Samp);
	idx_inbag.clear();
	idx_outbag.clear();
	idx_inbag.reserve(snpMat.Num_Total_Samp);
	idx_outbag.reserve(snpMat.Num_Total_Samp);
	// for-loop
	for (int i=0; i < snpMat.Num_Total_Samp; i++)
	{
		TGenotype &g = _GenoList.List[i];
		g.BootstrapCount = _BootstrapCnt[i];
		g.aux_hla_type = hlaList.List[i];
		if (g.aux_hla_type.Allele2 < g.aux_hla_type.Allele1)
		{
			int w = g.aux_hla_type.Allele2;
			g.aux_hla_type.Allele2 = g.aux_hla_type.Allele1;
			g.aux_hla_type.Allele1 = w;
		}
		// push back indices
		if (_BootstrapCnt[i] > 0)
			idx_inbag.push_back(i);
		else
			idx_outbag.push_back(i);
	}
	_GenoList.Num_SNP = 0;
	_GenoList.SetAllMissing();
	_Predict.InitPrediction(nHLA());
}

void CVariableSelection::_InitHaplotype(CHaplotypeList &Haplo)
{
	const size_t n_hla = _HLAList->Num_HLA_Allele();
	vector<int> tmp(n_hla, 0);
	int SumCnt = 0;
	for (vector<int>::const_iterator p=idx_inbag.begin(); p != idx_inbag.end(); p++)
	{
		const TGenotype &G = _GenoList.List[*p];
		const int cnt = G.BootstrapCount;
		tmp[G.aux_hla_type.Allele1] += cnt;
		tmp[G.aux_hla_type.Allele2] += cnt;
		SumCnt += cnt;
	}

	Haplo.LenPerHLA.resize(n_hla);
	int n_valid = 0;
	for (size_t i=0; i < n_hla; i++)
	{
		if (tmp[i] > 0) n_valid ++;
		Haplo.LenPerHLA[i] = (tmp[i] > 0) ? 1 : 0;
	}

	Haplo.Num_SNP = 0;
	Haplo.ResizeHaplo(n_valid);
	const double scale = 0.5 / SumCnt;
	n_valid = 0;
	for (size_t i=0; i < n_hla; i++)
	{
		if (tmp[i] > 0)
			Haplo.List[n_valid++].Freq = tmp[i] * scale;
	}
}

void CVariableSelection::_Init_EvalAcc(CHaplotypeList &Haplo,
	CGenotypeList &Geno)
{
	if (GPUExtProcPtr && *GPUExtProcPtr->build_set_haplo_geno)
	{
		Haplo.SetHaploAux_GPU();
		(*GPUExtProcPtr->build_set_haplo_geno)(Haplo.List, Haplo.Num_Haplo,
			&Geno.List[0], Haplo.Num_SNP);
	} else {
		if (need_auxiliary_haplo)
		{
			aux_haplo.resize(Haplo.Num_Haplo*2);
			aux_freq.resize(Haplo.Num_Haplo);
			Haplo.SetHaploAux(&aux_haplo[0], &aux_freq[0]);
		}
	}
}

void CVariableSelection::_Done_EvalAcc()
{ }

int CVariableSelection::_OutOfBagAccuracy(CHaplotypeList &Haplo)
{
	HIBAG_CHECKING(Haplo.Num_SNP != _GenoList.Num_SNP,
		"CVariableSelection::_OutOfBagAccuracy, Haplo and GenoList should have the same number of SNP markers.");
	if (GPUExtProcPtr && *GPUExtProcPtr->build_acc_oob)
	{
		return (*GPUExtProcPtr->build_acc_oob)();
	} else {
		const int nthread = thread_num();
		vector<int> cnt(nthread, 0);
		PARALLEL_FOR(i, idx_outbag.size())
		{
			TGenotype &p = _GenoList.List[idx_outbag[i]];
			THLAType g = (*fc_BestGuess)(Haplo, p);
			cnt[thread_idx()] += CHLATypeList::Compare(g, p.aux_hla_type);
		}
		PARALLEL_END
		int CorrectCnt = 0;
		for (int i=0; i < nthread; i++) CorrectCnt += cnt[i];
		return CorrectCnt;
	}
}

double CVariableSelection::_InBagLogLik(CHaplotypeList &Haplo)
{
	HIBAG_CHECKING(Haplo.Num_SNP != _GenoList.Num_SNP,
		"CVariableSelection::_InBagLogLik, Haplo and GenoList should have the same number of SNP markers.");
	if (GPUExtProcPtr && *GPUExtProcPtr->build_acc_ib)
	{
		return (*GPUExtProcPtr->build_acc_ib)();
	} else {
		const size_t n = idx_inbag.size();
		PARALLEL_FOR(i, n)
		{
			TGenotype &p = _GenoList.List[idx_inbag[i]];
			log_inbag[i] = p.BootstrapCount *
				log((*fc_PostProb)(Haplo, p, p.aux_hla_type));
		}
		PARALLEL_END
		// reduce sum
		double LogLik = 0;
		for (size_t i=0; i < n; i++) LogLik += log_inbag[i];
		LogLik *= -2;
		return LogLik;
	}
}

void CVariableSelection::Search(CBaseSampling &VarSampling,
	CHaplotypeList &OutHaplo, vector<int> &OutSNPIndex,
	double &Out_Global_Max_OutOfBagAcc, int mtry, bool prune,
	bool verbose, bool verbose_detail)
{
	// rare probability
	const double RARE_PROB = std::max(FRACTION_HAPLO/(2*nSamp()), MIN_RARE_FREQ);

	// initialize output
	_InitHaplotype(OutHaplo);
	OutSNPIndex.clear();

	// initialize internal variables
	const int NumOOB = idx_outbag.size();
	int Global_Max_OutOfBagAcc = 0;  // # of correct alleles
	double Global_Min_Loss = 1e+30;

	// reserve memory for haplotype lists
	const size_t reserve_num_haplo = nSamp() * 2;
	CHaplotypeList NextHaplo(reserve_num_haplo);
	CHaplotypeList NextReducedHaplo(reserve_num_haplo);
	CHaplotypeList MinHaplo(reserve_num_haplo);

	while (VarSampling.TotalNum() > 0 &&
		OutSNPIndex.size() < HIBAG_MAXNUM_SNP_IN_CLASSIFIER)
	{
		// prepare for growing the individual classifier
		_EM.PrepareHaplotypes(OutHaplo, _GenoList, NextHaplo);

		int max_OutOfBagAcc = Global_Max_OutOfBagAcc;
		double min_loss = Global_Min_Loss;
		int min_i = -1;

		// sample mtry from all candidate SNP markers
		VarSampling.RandomSelect(mtry);

		// for-loop
		for (int i=0; i < VarSampling.NumOfSelection(); i++)
		{
			if (_EM.PrepareNewSNP(VarSampling[i], OutHaplo, *_SNPMat, _GenoList, NextHaplo))
			{
				// run EM algorithm
				_EM.ExpectationMaximization(NextHaplo);
				// remove rare haplotypes
				NextHaplo.EraseDoubleHaplos(RARE_PROB, NextReducedHaplo);
				// add a SNP to the SNP genotype list
				_GenoList.AddSNP(VarSampling[i], *_SNPMat);

				// evaluate losses
				_Init_EvalAcc(NextReducedHaplo, _GenoList);
				double loss = 0;
				int acc = _OutOfBagAccuracy(NextReducedHaplo);
				if (acc >= max_OutOfBagAcc)
					loss = _InBagLogLik(NextReducedHaplo);
				_Done_EvalAcc();

				// remove the last SNP in the SNP genotype list
				_GenoList.ReduceSNP();

				// compare
				if (acc > max_OutOfBagAcc)
				{
					min_i = i;
					min_loss = loss;
					max_OutOfBagAcc = acc;
					MinHaplo = NextReducedHaplo;
				} else if (acc == max_OutOfBagAcc)
				{
					if (loss < min_loss)
					{
						min_i = i;
						min_loss = loss;
						MinHaplo = NextReducedHaplo;
					}
				}
				// check and delete
				if (prune)
				{
					if (acc < Global_Max_OutOfBagAcc)
					{
						VarSampling[i] = -1;
					} else if (acc == Global_Max_OutOfBagAcc)
					{
						if ((loss > Global_Min_Loss*(1+PRUNE_RELTOL_LOGLIK)) && (min_i != i))
							VarSampling[i] = -1;
					}
				}
			}
		}

		// compare ...
		bool sign = false;
		if (max_OutOfBagAcc > Global_Max_OutOfBagAcc)
		{
			sign = true;
		} else if (max_OutOfBagAcc == Global_Max_OutOfBagAcc)
		{
			if (min_i >= 0)
			{
				sign = ((min_loss >= STOP_RELTOL_LOGLIK_ADDSNP) &&
					(min_loss < Global_Min_Loss*(1-STOP_RELTOL_LOGLIK_ADDSNP)));
			} else
				sign = false;
		} else
			sign = false;

		// handle ...
		if (sign)
		{
			// add a new SNP predictor
			Global_Max_OutOfBagAcc = max_OutOfBagAcc;
			Global_Min_Loss = min_loss;
			OutHaplo = MinHaplo;
			OutSNPIndex.push_back(VarSampling[min_i]);
			_GenoList.AddSNP(VarSampling[min_i], *_SNPMat);
			if (prune)
			{
				VarSampling[min_i] = -1;
				VarSampling.RemoveFlag();
			} else {
				VarSampling.Remove(min_i);
			}
			// show ...
			if (verbose_detail)
			{
				Rprintf("    %2d, SNP: %d, loss: %g, oob acc: %0.2f%%, # of haplo: %d\n",
					OutSNPIndex.size(), OutSNPIndex.back()+1,
					Global_Min_Loss,
					double(Global_Max_OutOfBagAcc) / NumOOB * 50,
					OutHaplo.Num_Haplo);
			}
			CheckInterrupt();
		} else {
			// only keep "n_tmp - m" predictors
			VarSampling.RemoveSelection();
			// remove the last SNP in the SNP genotype list (set it to missing)
			_GenoList.SetMissing(_GenoList.Num_SNP);
		}
	}

	Out_Global_Max_OutOfBagAcc = 0.5 * Global_Max_OutOfBagAcc / NumOOB;
}



// -------------------------------------------------------------------------
// The individual classifier

CAttrBag_Classifier::CAttrBag_Classifier(CAttrBag_Model &_owner)
{
	_Owner = &_owner;
	_OutOfBag_Accuracy = 0;
}

void CAttrBag_Classifier::InitBootstrapCount(const int SampCnt[])
{
	_BootstrapCount.assign(&SampCnt[0], &SampCnt[_Owner->nSamp()]);
	_SNPIndex.clear();
	_OutOfBag_Accuracy = 0;
}

void CAttrBag_Classifier::Assign(int n_snp, const int snpidx[],
	const int samp_num[], int n_haplo, const double *freq, const int *hla,
	const char *haplo[], double *_acc)
{
	// SNP markers
	_SNPIndex.assign(&snpidx[0], &snpidx[n_snp]);
	// The number of samples
	if (samp_num)
	{
		const int n = _Owner->nSamp();
		_BootstrapCount.assign(&samp_num[0], &samp_num[n]);
	}
	// The haplotypes
	_Haplo.Num_SNP = n_snp;
	_Haplo.ResizeHaplo(n_haplo);
	_Haplo.LenPerHLA.resize(_Owner->nHLA());
	for (int i=0; i < n_haplo; i++)
	{
		_Haplo.List[i] = THaplotype(haplo[i], freq[i]);
		_Haplo.LenPerHLA[hla[i]] ++;
	}
	// Accuracies
	_OutOfBag_Accuracy = (_acc) ? (*_acc) : 0;
}

void CAttrBag_Classifier::Grow(CBaseSampling &VarSampling, int mtry,
	bool prune, bool verbose, bool verbose_detail)
{
	_Owner->_VarSelect.InitSelection(_Owner->_SNPMat,
		_Owner->_HLAList, &_BootstrapCount[0]);
	_Owner->_VarSelect.Search(VarSampling, _Haplo, _SNPIndex,
		_OutOfBag_Accuracy, mtry, prune, verbose, verbose_detail);
}


// -------------------------------------------------------------------------
// the attribute bagging model

CAttrBag_Model::CAttrBag_Model() { }

void CAttrBag_Model::InitTraining(int n_snp, int n_samp, int n_hla)
{
	HIBAG_CHECKING(n_snp < 0, "CAttrBag_Model::InitTraining, n_snp error.")
	HIBAG_CHECKING(n_samp < 0, "CAttrBag_Model::InitTraining, n_samp error.")
	HIBAG_CHECKING(n_hla < 0, "CAttrBag_Model::InitTraining, n_hla error.")

	_SNPMat.Num_Total_Samp = n_samp;
	_SNPMat.Num_Total_SNP = n_snp;
	_SNPMat.pGeno = NULL;

	_HLAList.List.resize(n_samp);
	_HLAList.Str_HLA_Allele.resize(n_hla);
}

void CAttrBag_Model::InitTraining(int n_snp, int n_samp, int *snp_geno,
	int n_hla, int *H1, int *H2)
{
	HIBAG_CHECKING(n_snp < 0, "CAttrBag_Model::InitTraining, n_snp error.")
	HIBAG_CHECKING(n_samp < 0, "CAttrBag_Model::InitTraining, n_samp error.")
	HIBAG_CHECKING(n_hla < 0, "CAttrBag_Model::InitTraining, n_hla error.")

	_SNPMat.Num_Total_Samp = n_samp;
	_SNPMat.Num_Total_SNP = n_snp;
	_SNPMat.pGeno = snp_geno;

	_HLAList.List.resize(n_samp);
	_HLAList.Str_HLA_Allele.resize(n_hla);
	for (int i=0; i < n_samp; i++)
	{
		HIBAG_CHECKING(H1[i]<0 || H1[i]>=n_hla,
			"CAttrBag_Model::InitTraining, H1 error.");
		HIBAG_CHECKING(H2[i]<0 || H2[i]>=n_hla,
			"CAttrBag_Model::InitTraining, H2 error.");
		_HLAList.List[i].Allele1 = H1[i];
		_HLAList.List[i].Allele2 = H2[i];
	}
}

CAttrBag_Classifier *CAttrBag_Model::NewClassifierBootstrap()
{
	_ClassifierList.push_back(CAttrBag_Classifier(*this));
	CAttrBag_Classifier *I = &_ClassifierList.back();

	const int n = nSamp();
	vector<int> S(n);
	int n_unique;

	do {
		// initialize S
		for (int i=0; i < n; i++) S[i] = 0;
		n_unique = 0;

		for (int i=0; i < n; i++)
		{
			int k = RandomNum(n);
			if (S[k] == 0) n_unique ++;
			S[k] ++;
		}
	} while (n_unique >= n); // to avoid the case of no out-of-bag individuals

	I->InitBootstrapCount(&S[0]);

	return I;
}

CAttrBag_Classifier *CAttrBag_Model::NewClassifierAllSamp()
{
	_ClassifierList.push_back(CAttrBag_Classifier(*this));
	CAttrBag_Classifier *I = &_ClassifierList.back();
	vector<int> S(nSamp(), 1);
	I->InitBootstrapCount(&S[0]);
	return I;
}

CAttrBag_Model::try_final_train_gpu::try_final_train_gpu(CAttrBag_Model *m)
{
	if (GPUExtProcPtr && *GPUExtProcPtr->build_init)
		(*GPUExtProcPtr->build_init)(m->nHLA(), m->nSamp());
}

CAttrBag_Model::try_final_train_gpu::~try_final_train_gpu()
{
	if (GPUExtProcPtr && *GPUExtProcPtr->build_done)
		(*GPUExtProcPtr->build_done)();
}

void CAttrBag_Model::BuildClassifiers(int nclassifier, int mtry, bool prune,
	bool verbose, bool verbose_detail)
{
	try_final_train_gpu gpu(this);
	CSamplingWithoutReplace VarSampling;

	for (int k=0; k < nclassifier; k++)
	{
		VarSampling.Init(nSNP());

		CAttrBag_Classifier *I = NewClassifierBootstrap();
		if (verbose)
		{
			const vector<int> &bc = I->BootstrapCount();
			int nOOB = 0;
			for (size_t i=0; i < bc.size(); i++)
				if (bc[i] == 0) nOOB++;
			Rprintf("=== building individual classifier %d, out-of-bag (%d/%.1f%%) ===\n",
				(int)_ClassifierList.size(), nOOB, 100.0*nOOB/bc.size());
		}

		// initialize bootstrap samples in GPU implementation
		if (GPUExtProcPtr && *GPUExtProcPtr->build_set_bootstrap)
		{
			(*GPUExtProcPtr->build_set_bootstrap)(&(I->BootstrapCount()[0]));
		}

		I->Grow(VarSampling, mtry, prune, verbose, verbose_detail);
		if (verbose)
		{
			Rprintf(
				"[%d] %s, oob acc: %0.2f%%, # of SNPs: %d, # of haplo: %d\n",
				(int)_ClassifierList.size(), date_text(),
				I->OutOfBag_Accuracy()*100, I->nSNP(), I->nHaplo());
			CheckInterrupt();
		}
	}
}

CAttrBag_Model::try_final_pred_gpu::try_final_pred_gpu(CAttrBag_Model *m)
{
	(owner = m)->_Init_GPU_PredHLA();
}

CAttrBag_Model::try_final_pred_gpu::~try_final_pred_gpu()
{
	owner->_Done_GPU_PredHLA();
}

void CAttrBag_Model::PredictHLA(const int *genomat, int n_samp, int vote_method,
	int OutH1[], int OutH2[], double OutMaxProb[], double OutMatching[],
	double OutDosage[], double OutProbArray[], bool verbose)
{
	if ((vote_method < 1) || (vote_method > 2))
		throw ErrHLA("Invalid 'vote_method'.");

	// prediction object
	const int nthread = thread_num();
	vector<CAlg_Prediction> pred_lst(nthread);
	for (int i=0; i < nthread; i++)
		pred_lst[i].InitPrediction(nHLA());

	// set auxiliary variables if need_auxiliary_haplo=true
	vector<int64_t> aux_haplo;
	vector<double> aux_freq;
	if (need_auxiliary_haplo)
	{
		const int ncl = _ClassifierList.size();
		size_t sum_num_haplo = 0;
		for (int i=0; i < ncl; i++)
			sum_num_haplo += _ClassifierList[i]._Haplo.Num_Haplo;
		aux_haplo.resize(2*sum_num_haplo);
		aux_freq.resize(sum_num_haplo);
		// call SetHaploAux
		size_t p = 0;
		for (int i=0; i < ncl; i++)
		{
			CHaplotypeList &haplo = _ClassifierList[i]._Haplo;
			haplo.SetHaploAux(&aux_haplo[2*p], &aux_freq[p]);
			p += haplo.Num_Haplo;
		}
	}

	// get SNP weights according to the number of classifiers using that SNP
	vector<double> c_weight(_ClassifierList.size()*nthread);
	vector<int> snp_weight(nSNP());
	_GetSNPWeights(&snp_weight[0]);
	const size_t nn = nHLA()*(nHLA()+1)/2;

	// progress information
	Progress.Info = "Predicting";
	Progress.Init(n_samp, verbose);

	try_final_pred_gpu gpu(this);
	PARALLEL_FOR(i, n_samp)
	{
		const int th_idx = thread_idx();
		CAlg_Prediction &pred = pred_lst[th_idx];
		double match_prob;
		_PredictHLA(pred, genomat + i*nSNP(), &snp_weight[0], vote_method,
			match_prob, &c_weight[th_idx*_ClassifierList.size()]);

		THLAType HLA = pred.BestGuessEnsemble();
		if (OutH1 && OutH2)
		{
			OutH1[i] = HLA.Allele1;
			OutH2[i] = HLA.Allele2;
		}
		if (OutMaxProb)
		{
			if ((HLA.Allele1 != NA_INTEGER) && (HLA.Allele2 != NA_INTEGER))
				OutMaxProb[i] = pred.IndexSumPostProb(HLA.Allele1, HLA.Allele2);
			else
				OutMaxProb[i] = 0;
		}
		if (OutMatching)  // matching proportion
		{
			OutMatching[i] = match_prob;
		}
		if (OutDosage)  // dosages
		{
			const size_t n = nHLA();
			const double *s = &pred._SumPostProb[0];
			double *p = OutDosage + i*n;
			memset(p, 0, sizeof(double)*n);
			for (size_t h1=0; h1 < n; h1++)
			{
				p[h1] += 2 * (*s++);  // expected hom. dosage
				for (size_t h2=h1+1; h2 < n; h2++)
				{
					const double v = *s++;
					p[h1] += v; p[h2] += v;
				}
			}
		}
		if (OutProbArray)  // probabilities for all allele pairs
		{
			memcpy(OutProbArray+i*nn, &pred._SumPostProb[0], sizeof(double)*nn);
		}

		Progress.Forward(1, verbose);
		if (th_idx == 0) CheckInterrupt();  // run on the baseline thread
	}
	PARALLEL_END
}

void CAttrBag_Model::_PredictHLA(CAlg_Prediction &pred, const int geno[],
	const int snp_weight[], int vote_method, double &OutMatching, double c_weight[])
{
	// weight for each classifier, based on missing proportion
	vector<CAttrBag_Classifier>::iterator p = _ClassifierList.begin();
	for (size_t w_i=0; p != _ClassifierList.end(); p++)
	{
		const int n = p->nSNP();
		int nw=0, sum=0;
		for (int i=0; i < n; i++)
		{
			int k = p->_SNPIndex[i];
			sum += snp_weight[k];
			if ((0 <= geno[k]) && (geno[k] <= 2))
				nw += snp_weight[k];
		}
		c_weight[w_i++] = (sum > 0) ? (double(nw) / sum) : 0;
	}

	if (GPUExtProcPtr && *GPUExtProcPtr->predict_avg_prob)
	{
		p = _ClassifierList.begin();
		for (size_t w_i=0; p != _ClassifierList.end(); p++, w_i++)
		{
			gpu_geno_buf[w_i].IntToSNP(p->nSNP(), geno, &(p->_SNPIndex[0]));
		}
		(*GPUExtProcPtr->predict_avg_prob)(&gpu_geno_buf[0], c_weight,
			&pred._SumPostProb[0], &OutMatching);
	} else {
		// initialize probability
		pred.InitSumPostProb();
		TGenotype Geno;
		double sum_matching=0, num_matching=0;

		p = _ClassifierList.begin();
		for (size_t w_i=0; p != _ClassifierList.end(); p++, w_i++)
		{
			if (c_weight[w_i] <= 0) continue;
			// get the packed SNP genotypes
			Geno.IntToSNP(p->nSNP(), geno, &(p->_SNPIndex[0]));
			// predict
			double pm;
			pred.PredictPostProb(p->_Haplo, Geno, pm);
			// add matching probability
			sum_matching += pm * c_weight[w_i];
			num_matching += c_weight[w_i];
			// add prob to the ensemble
			if (vote_method == 1)
			{
				// predicting based on the averaged posterior probabilities
				pred.AddProbToSum(c_weight[w_i]);
			} else if (vote_method == 2)
			{
				// predicting by class majority voting
				THLAType pd = pred.BestGuess();
				if ((pd.Allele1 != NA_INTEGER) && (pd.Allele2 != NA_INTEGER))
				{
					pred.InitPostProb();  // fill by ZERO
					pred.IndexPostProb(pd.Allele1, pd.Allele2) = 1.0;
					pred.AddProbToSum(1.0);
				}
			}
		}

		// normalize the sum of posterior prob
		pred.NormalizeSumPostProb();
		OutMatching = sum_matching / num_matching;
	}
}

void CAttrBag_Model::_GetSNPWeights(int OutSNPWeight[])
{
	// initialize
	memset(OutSNPWeight, 0, sizeof(int)*nSNP());
	// for each classifier
	vector<CAttrBag_Classifier>::const_iterator p;
	for (p = _ClassifierList.begin(); p != _ClassifierList.end(); p++)
	{
		const size_t n = p->nSNP();
		for (size_t i=0; i < n; i++)
			OutSNPWeight[ p->_SNPIndex[i] ] ++;
	}
}

void CAttrBag_Model::_Init_GPU_PredHLA()
{
	if (GPUExtProcPtr && *GPUExtProcPtr->predict_init)
	{
		// prepare data structure for GPU
		const size_t n_classifier = _ClassifierList.size();
		THaplotype* haplo[n_classifier];
		int p_n_haplo[n_classifier];
		int p_n_snp[n_classifier];

		gpu_geno_buf.resize(n_classifier);
		vector<CAttrBag_Classifier>::iterator p = _ClassifierList.begin();
		for (size_t c_i=0; p != _ClassifierList.end(); p++, c_i++)
		{
			CHaplotypeList &hl = p->_Haplo;
			hl.SetHaploAux_GPU();
			haplo[c_i] = hl.List;
			p_n_haplo[c_i] = p->nHaplo();
			p_n_snp[c_i] = p->nSNP();
		}

		// call the external function
		(*GPUExtProcPtr->predict_init)(nHLA(), n_classifier, haplo,
			p_n_haplo, p_n_snp);
	}
}

void CAttrBag_Model::_Done_GPU_PredHLA()
{
	if (GPUExtProcPtr && *GPUExtProcPtr->predict_done)
	{
		(*GPUExtProcPtr->predict_done)();
	}
}
