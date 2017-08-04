// ===============================================================
//
// HIBAG R package (HLA Genotype Imputation with Attribute Bagging)
// Copyright (C) 2011-2017   Xiuwen Zheng (zhengx@u.washington.edu)
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
// Kernel Version : 1.4
// Copyright      : Xiuwen Zheng (GPL v3)
// Description    : HLA imputation C++ library
// ===============================================================


#include "LibHLA.h"

#define HIBAG_ENABLE_TIMING

#ifdef HIBAG_ENABLE_TIMING
#   include <time.h>
#endif


using namespace std;
using namespace HLA_LIB;



// ========================================================================= //

#ifdef HIBAG_ENABLE_TIMING

static clock_t timing_array[5];

template<size_t I> struct TTiming
{
	clock_t start;
	inline TTiming() { start = clock(); }
	inline ~TTiming() { Stop(); }
	inline void Stop() { timing_array[I] += clock() - start; }
};

#define HIBAG_TIMING(i)    TTiming<i> tm;
#define TM_TOTAL      0
#define TM_ACC_OOB    1
#define TM_ACC_IB     2
#define TM_PRE_HAPLO  3
#define TM_EM_ALG     4

// _OutOfBagAccuracy: 41.00%
// _InBagLogLik: 54.58%
//  PrepareHaplotypes: 0.34%
//  ExpectationMaximization: 3.30%

#else
#   define HIBAG_TIMING(i)
#endif



// ========================================================================= //

// Parameters -- EM algorithm

/// the max number of iterations
int HLA_LIB::EM_MaxNum_Iterations = 500;
/// the initial value of EM algorithm
static const double EM_INIT_VAL_FRAC = 0.001;
/// the reltol convergence tolerance, sqrt(machine.epsilon) by default, used in EM algorithm
double HLA_LIB::EM_FuncRelTol = sqrt(DBL_EPSILON);


// Parameters -- reduce the number of possible haplotypes

/// The minimum rare frequency to store haplotypes
static const double MIN_RARE_FREQ = 1e-5;
/// The fraction of one haplotype that can be ignored
static const double FRACTION_HAPLO = 1.0/10;


// Parameters -- search SNP markers

/// the reltol for the stopping rule of adding a new SNP marker
static const double STOP_RELTOL_LOGLIK_ADDSNP = 0.001;
/// the reltol for erasing the SNP marker is prune = TRUE
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

/// Frequency Calculation
#define FREQ_MUTANT(p, cnt)    ((p) * EXP_LOG_MIN_RARE_FREQ[cnt])

/// exp(cnt * log(MIN_RARE_FREQ)), cnt is the hamming distance
static double EXP_LOG_MIN_RARE_FREQ[HIBAG_MAXNUM_SNP_IN_CLASSIFIER*2];

class CInit
{
public:
	CInit()
	{
		const int n = 2 * HIBAG_MAXNUM_SNP_IN_CLASSIFIER;
		for (int i=0; i < n; i++)
			EXP_LOG_MIN_RARE_FREQ[i] = exp(i * log(MIN_RARE_FREQ));
		EXP_LOG_MIN_RARE_FREQ[0] = 1;
		for (int i=0; i < n; i++)
		{
			if (!R_finite(EXP_LOG_MIN_RARE_FREQ[i]))
				EXP_LOG_MIN_RARE_FREQ[i] = 0;
		}
	}
};

static CInit _Init;


// GPU extensible component
TypeGPUExtProc *HLA_LIB::GPUExtProcPtr = NULL;;



// ========================================================================= //
// CdProgression

static char date_buffer[256];

inline static const char *date_text()
{
	time_t rawtime;
	time(&rawtime);
	struct tm *p = localtime(&rawtime);
	sprintf(date_buffer, "%04d-%02d-%02d %02d:%02d:%02d", p->tm_year+1900,
		p->tm_mon+1, p->tm_mday, p->tm_hour, p->tm_min, p->tm_sec);
	return date_buffer;
}

static const clock_t TimeInterval = 15*CLOCKS_PER_SEC;

CdProgression::CdProgression()
{
	Init(0, false);
}

void CdProgression::Init(long TotalCnt, bool ShowInit)
{
	if (TotalCnt < 0) TotalCnt = 0;
	fTotal = TotalCnt;
	fCurrent = fPercent = 0;
	OldTime = clock();
	if (ShowInit) ShowProgress();
}

bool CdProgression::Forward(long step, bool Show)
{
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
			return true;
		}
	}
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
	return (PackedHaplo[idx >> 3] >> (idx & 0x07)) & 0x01;
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
		for (size_t i=0; i < Length; i++)
			rv[i] = (GetAllele(i)==0) ? '0' : '1';
	}
	return rv;
}

void THaplotype::StrToHaplo(const string &str)
{
	HIBAG_CHECKING(str.size() > HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"THaplotype::StrToHaplo, the input string is too long.");
	for (size_t i=0; i < str.size(); i++)
	{
		char ch = str[i];
		HIBAG_CHECKING(ch!='0' && ch!='1',
			"THaplotype::StrToHaplo, the input string should be '0' or '1'");
		_SetAllele(i, ch-'0');
	}
}

inline void THaplotype::_SetAllele(size_t idx, UINT8 val)
{
	size_t r = idx & 0x07;
	UINT8 mask = ~(0x01 << r);
	UINT8 &ch = PackedHaplo[idx >> 3];
	ch = (ch & mask) | (val << r);
}



// ========================================================================= //
// Haplotype list with an HLA gene and SNP alleles

CHaplotypeList::CHaplotypeList()
{
	Num_Haplo = Num_SNP = 0;
	reserve_size = 0;
	base_ptr = NULL;
	List = NULL;
}

CHaplotypeList::CHaplotypeList(const CHaplotypeList &src)
{
	Num_Haplo = Num_SNP = 0;
	reserve_size = 0;
	base_ptr = NULL;
	List = NULL;
	*this = src;
}

CHaplotypeList::CHaplotypeList(size_t reserve_num)
{
	Num_Haplo = Num_SNP = 0;
	reserve_size = 0;
	base_ptr = NULL;
	List = NULL;
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

	// TODO
	int ss = 0;
	for (size_t i=0; i < OutHaplos.LenPerHLA.size(); i++)
	{
		ss += OutHaplos.LenPerHLA[i];
	}
	if (ss != OutHaplos.Num_Haplo)
		throw "Erroroeewrweq";

	OutHaplos.ScaleFrequency(1/sum);
}

void CHaplotypeList::SaveClearFrequency()
{
	THaplotype *p = List;
	for (size_t n=Num_Haplo; n > 0; n--)
	{
		p->aux.OldFreq = p->Freq;
		p->Freq = 0;
		p ++;
	}
}

void CHaplotypeList::ScaleFrequency(double scale)
{
	THaplotype *p = List;
	for (size_t n=Num_Haplo; n > 0; n--)
	{
		p->Freq *= scale;
		p ++;
	}
}

size_t CHaplotypeList::StartHaploHLA(int hla) const
{
	HIBAG_CHECKING(hla < 0 || hla >= (int)LenPerHLA.size(),
		"CHaplotypeList::StartHLA, invalid HLA allele.");
	size_t rv = 0;
	for (int i=0; i < hla; i++) rv += LenPerHLA[i];
	return rv;
}

void CHaplotypeList::SetHaploAux()
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
	size_t i = idx >> 3, r = idx & 0x07;
	if ((PackedMissing[i] >> r) & 0x01)
		return ((PackedSNP1[i] >> r) & 0x01) + ((PackedSNP2[i] >> r) & 0x01);
	else
		return -1;
}

void TGenotype::SetSNP(size_t idx, int val)
{
	HIBAG_CHECKING(idx >= HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"TGenotype::SetSNP, invalid index.");
	_SetSNP(idx, val);
}

void TGenotype::_SetSNP(size_t idx, int val)
{
	size_t i = idx >> 3, r = idx & 0x07;
	UINT8 &S1 = PackedSNP1[i];
	UINT8 &S2 = PackedSNP2[i];
	UINT8 &M  = PackedMissing[i];
	UINT8 SET = (UINT8(0x01) << r);
	UINT8 CLEAR = ~SET;

	switch (val)
	{
		case 0:
			S1 &= CLEAR; S2 &= CLEAR; M |= SET; break;
		case 1:
			S1 |= SET; S2 &= CLEAR; M |= SET; break;
		case 2:
			S1 |= SET; S2 |= SET; M |= SET; break;
		default:
			S1 &= CLEAR; S2 &= CLEAR; M &= CLEAR; break;
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
		char ch = str[i];
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

void TGenotype::IntToSNP(size_t Length, const int InBase[], const int Index[])
{
	HIBAG_CHECKING(Length > HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"TGenotype::IntToSNP, the length is invalid.");

	const static UINT8 P1[4] = { 0, 1, 1, 0 };
	const static UINT8 P2[4] = { 0, 0, 1, 0 };
	const static UINT8 PM[4] = { 1, 1, 1, 0 };

	UINT8 *p1 = PackedSNP1;     // --> P1
	UINT8 *p2 = PackedSNP2;     // --> P2
	UINT8 *pM = PackedMissing;  // --> PM

	for (; Length >= 8; Length -= 8, Index += 8)
	{
		unsigned g1 = InBase[Index[0]];
		size_t i1 = (g1 < 3) ? g1 : 3;
		unsigned g2 = InBase[Index[1]];
		size_t i2 = (g2 < 3) ? g2 : 3;
		unsigned g3 = InBase[Index[2]];
		size_t i3 = (g3 < 3) ? g3 : 3;
		unsigned g4 = InBase[Index[3]];
		size_t i4 = (g4 < 3) ? g4 : 3;
		unsigned g5 = InBase[Index[4]];
		size_t i5 = (g5 < 3) ? g5 : 3;
		unsigned g6 = InBase[Index[5]];
		size_t i6 = (g6 < 3) ? g6 : 3;
		unsigned g7 = InBase[Index[6]];
		size_t i7 = (g7 < 3) ? g7 : 3;
		unsigned g8 = InBase[Index[7]];
		size_t i8 = (g8 < 3) ? g8 : 3;

		*p1++ = P1[i1] | (P1[i2] << 1) | (P1[i3] << 2) | (P1[i4] << 3) |
			(P1[i5] << 4) | (P1[i6] << 5) | (P1[i7] << 6) | (P1[i8] << 7);
		*p2++ = P2[i1] | (P2[i2] << 1) | (P2[i3] << 2) | (P2[i4] << 3) |
			(P2[i5] << 4) | (P2[i6] << 5) | (P2[i7] << 6) | (P2[i8] << 7);
		*pM++ = PM[i1] | (PM[i2] << 1) | (PM[i3] << 2) | (PM[i4] << 3) |
			(PM[i5] << 4) | (PM[i6] << 5) | (PM[i7] << 6) | (PM[i8] << 7);
	}

	if (Length > 0)
	{
		*p1 = *p2 = *pM = 0;
		for (size_t i=0; i < Length; i++)
		{
			unsigned g1 = InBase[*Index++];
			size_t i1 = (g1 < 3) ? g1 : 3;
			*p1 |= (P1[i1] << i);
			*p2 |= (P2[i1] << i);
			*pM |= (PM[i1] << i);
		}
		pM ++;
	}

	for (UINT8 *pEnd=PackedMissing+sizeof(PackedMissing); pM < pEnd; )
		*pM++ = 0;
}

int TGenotype::HammingDistance(size_t Length,
	const THaplotype &H1, const THaplotype &H2) const
{
	HIBAG_CHECKING(Length > HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"THaplotype::HammingDistance, the length is too large.");
	return _HamDist(Length, H1, H2);
}


// compute the Hamming distance between SNPs and H1+H2 without checking

#ifdef HIBAG_SIMD_OPTIMIZE_HAMMING_DISTANCE

	// signed integer for initializing XMM
	typedef int64_t UTYPE;
	#ifndef HIBAG_HARDWARE_POPCNT
	static const __m128i Z05 = _mm_set1_epi8(0x55);
	static const __m128i Z03 = _mm_set1_epi8(0x33);
	static const __m128i Z0F = _mm_set1_epi8(0x0F);
	#endif

#else
	#ifdef HIBAG_REG_BIT64
		typedef uint64_t UTYPE;
	#else
		typedef uint32_t UTYPE;
	#endif
#endif

static const ssize_t UTYPE_BIT_NUM = sizeof(UTYPE)*8;

inline int TGenotype::_HamDist(size_t Length,
	const THaplotype &H1, const THaplotype &H2) const
{
	size_t ans = 0;

	const UTYPE *h1 = (const UTYPE*)&H1.PackedHaplo[0];
	const UTYPE *h2 = (const UTYPE*)&H2.PackedHaplo[0];
	const UTYPE *s1 = (const UTYPE*)&PackedSNP1[0];
	const UTYPE *s2 = (const UTYPE*)&PackedSNP2[0];
	const UTYPE *sM = (const UTYPE*)&PackedMissing[0];

#ifdef HIBAG_SIMD_OPTIMIZE_HAMMING_DISTANCE

	// for-loop
	for (ssize_t n=Length; n > 0; n -= UTYPE_BIT_NUM)
	{
		__m128i H  = _mm_set_epi64x(*h2++, *h1++);  // *h1, *h2
		__m128i S1 = _mm_set_epi64x(*s2++, *s1++);  // *s1, *s2
		__m128i S2 = _mm_shuffle_epi32(S1, _MM_SHUFFLE(1,0,3,2)); // *s2, *s1

		__m128i mask1 = _mm_xor_si128(H, S2);
		__m128i mask2 = _mm_shuffle_epi32(mask1, _MM_SHUFFLE(1,0,3,2));

		// worry about n < UTYPE_BIT_NUM? unused bits have been set to zero
		__m128i M = _mm_set1_epi64x(*sM++);
		__m128i MASK = _mm_and_si128(_mm_or_si128(mask1, mask2), M);

		// val = '(H1 ^ S1) & MASK' / '(H2 ^ S2) & MASK'
		__m128i val = _mm_and_si128(_mm_xor_si128(H, S1), MASK);

		// popcount for val

	#ifdef HIBAG_HARDWARE_POPCNT

    #   ifdef HIBAG_REG_BIT64
			ans += _mm_popcnt_u64(M128_I64_0(val)) + _mm_popcnt_u64(M128_I64_1(val));
	#   else
			ans += _mm_popcnt_u32(M128_I32_0(val)) + _mm_popcnt_u32(M128_I32_1(val)) +
				_mm_popcnt_u32(M128_I32_2(val)) + _mm_popcnt_u32(M128_I32_3(val));
	#   endif

	#else

		// two 64-bit integers
		// suggested by
		// http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel

		// val -= ((val >> 1) & 0x5555555555555555);
		val = _mm_sub_epi64(val, _mm_and_si128(_mm_srli_epi64(val, 1), Z05));
		// val = (val & 0x3333333333333333) + ((val >> 2) & 0x3333333333333333);
		val = _mm_add_epi64(_mm_and_si128(val, Z03),
			_mm_and_si128(_mm_srli_epi64(val, 2), Z03));
		// val = (val + (val >> 4)) & 0x0F0F0F0F0F0F0F0F
		val = _mm_and_si128(_mm_add_epi64(val, _mm_srli_epi64(val, 4)), Z0F);

		// ans += (val * 0x0101010101010101LLU) >> 56;
		uint64_t r0 = _mm_cvtsi128_si64(val);
		uint64_t r1 = _mm_cvtsi128_si64(_mm_unpackhi_epi64(val, val));
		ans += ((r0 * 0x0101010101010101LLU) >> 56) +
			((r1 * 0x0101010101010101LLU) >> 56);

	#endif
	}

#else

	// for-loop
	for (ssize_t n=Length; n > 0; n -= UTYPE_BIT_NUM)
	{
		UTYPE H1 = *h1++;
		UTYPE H2 = *h2++;
		UTYPE S1 = *s1++;
		UTYPE S2 = *s2++;

		// worry about n < UTYPE_BIT_NUM? unused bits have been set to zero
		UTYPE M  = *sM++;  // missing value
		UTYPE MASK = ((H1 ^ S2) | (H2 ^ S1)) & M;

		// popcount for '(H1 ^ S1) & MASK'
		// suggested by
		// http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel

		UTYPE v1 = (H1 ^ S1) & MASK;
	#ifdef HIBAG_REG_BIT64
		// 64-bit integers
		v1 -= ((v1 >> 1) & 0x5555555555555555LLU);
		v1 = (v1 & 0x3333333333333333LLU) + ((v1 >> 2) & 0x3333333333333333LLU);
		ans += (((v1 + (v1 >> 4)) & 0x0F0F0F0F0F0F0F0FLLU) *
			0x0101010101010101LLU) >> 56;
	#else
		// 32-bit integers
		v1 -= ((v1 >> 1) & 0x55555555);
		v1 = (v1 & 0x33333333) + ((v1 >> 2) & 0x33333333);
		ans += (((v1 + (v1 >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
	#endif

		// popcount for '(H2 ^ S2) & MASK'
		UTYPE v2 = (H2 ^ S2) & MASK;
	#ifdef HIBAG_REG_BIT64
		// 64-bit integers
		v2 -= ((v2 >> 1) & 0x5555555555555555LLU);
		v2 = (v2 & 0x3333333333333333LLU) + ((v2 >> 2) & 0x3333333333333333LLU);
		ans += (((v2 + (v2 >> 4)) & 0x0F0F0F0F0F0F0F0FLLU) *
			0x0101010101010101LLU) >> 56;
	#else
		// 32-bit integers
		v2 -= ((v2 >> 1) & 0x55555555);
		v2 = (v2 & 0x33333333) + ((v2 >> 2) & 0x33333333);
		ans += (((v2 + (v2 >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
	#endif
	}

#endif

	return ans;
}



// -------------------------------------------------------------------------
// The class of SNP genotype list

CSNPGenoMatrix::CSNPGenoMatrix()
{
	Num_Total_SNP = Num_Total_Samp = 0;
	pGeno = NULL;
}

const int CSNPGenoMatrix::Get(const int IdxSamp, const int IdxSNP) const
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
#ifdef HIBAG_SIMD_OPTIMIZE_HAMMING_DISTANCE
	// since HIBAG_MAXNUM_SNP_IN_CLASSIFIER = 128
	__m128i zero = _mm_setzero_si128();
#endif
	for (TGenotype *p = &List[0]; n > 0; n--)
	{
	#ifdef HIBAG_SIMD_OPTIMIZE_HAMMING_DISTANCE
		_mm_storeu_si128((__m128i*)p->PackedMissing, zero);
	#else
		memset(p->PackedMissing, 0, sizeof(p->PackedMissing));
	#endif
		p ++;
	}
}

void CGenotypeList::SetMissing(int idx)
{
	size_t i = idx >> 3, r = idx & 0x07;
	UINT8 CLEAR = ~(UINT8(0x01) << r);
	size_t n = List.size();
	for (TGenotype *p = &List[0]; n > 0; n--)
	{
		p->PackedMissing[i] &= CLEAR;
		p ++;
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
	if ((P2==T1) || (P2==T2)) cnt ++;
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

CAlg_EM::CAlg_EM() {}

void CAlg_EM::PrepareHaplotypes(const CHaplotypeList &CurHaplo,
	const CGenotypeList &GenoList, const CHLATypeList &HLAList,
	CHaplotypeList &NextHaplo)
{
	HIBAG_TIMING(TM_PRE_HAPLO)
	HIBAG_CHECKING(GenoList.nSamp() != HLAList.nSamp(),
		"CAlg_EM::PrepareHaplotypes, GenoList and HLAList should have the same number of samples.");

	_SampHaploPair.clear();
	_SampHaploPair.reserve(GenoList.nSamp());
	CurHaplo.DoubleHaplos(NextHaplo);

	vector<int> DiffList(GenoList.nSamp()*(2*GenoList.nSamp() + 1));

	// get haplotype pairs for each sample
	for (int iSamp=0; iSamp < GenoList.nSamp(); iSamp++)
	{
		const TGenotype &pG   = GenoList.List[iSamp];
		const THLAType  &pHLA = HLAList.List[iSamp];

		if (pG.BootstrapCount > 0)
		{
			_SampHaploPair.push_back(THaploPairList());
			THaploPairList &HP = _SampHaploPair.back();
			HP.BootstrapCount = pG.BootstrapCount;
			HP.SampIndex = iSamp;

			size_t pH1_st = NextHaplo.StartHaploHLA(pHLA.Allele1);
			size_t pH1_n  = NextHaplo.LenPerHLA[pHLA.Allele1];
			size_t pH2_st = NextHaplo.StartHaploHLA(pHLA.Allele2);
			size_t pH2_n  = NextHaplo.LenPerHLA[pHLA.Allele2];
			THaplotype *p1, *p2;
			int MinDiff = GenoList.Num_SNP * 4;

			if (pHLA.Allele1 != pHLA.Allele2)
			{
				const size_t m = pH1_n * pH2_n;
				if (m > DiffList.size()) DiffList.resize(m);
				int *pD = &DiffList[0];

				p1 = &NextHaplo.List[pH1_st];
				for (size_t n1=pH1_n; n1 > 0; n1--, p1++)
				{
					p2 = &NextHaplo.List[pH2_st];
					for (size_t n2=pH2_n; n2 > 0; n2--, p2++)
					{
						int d = *pD++ = pG._HamDist(CurHaplo.Num_SNP, *p1, *p2);
						if (d < MinDiff) MinDiff = d;
						if (d == 0)
							HP.PairList.push_back(THaploPair(p1, p2));
					}
				}

				if (MinDiff > 0)
				{
					int *pD = &DiffList[0];
					p1 = &NextHaplo.List[pH1_st];
					for (size_t n1=pH1_n; n1 > 0; n1--, p1++)
					{
						p2 = &NextHaplo.List[pH2_st];
						for (size_t n2=pH2_n; n2 > 0; n2--, p2++)
						{
							if (*pD++ == MinDiff)
								HP.PairList.push_back(THaploPair(p1, p2));
						}
					}
				}

			} else {
				const size_t m = pH1_n * (pH1_n + 1) / 2;
				if (m > DiffList.size()) DiffList.resize(m);
				int *pD = &DiffList[0];

				p1 = &NextHaplo.List[pH1_st];
				for (size_t n1=pH1_n; n1 > 0; n1--, p1++)
				{
					p2 = p1;
					for (size_t n2=n1; n2 > 0; n2--, p2++)
					{
						int d = *pD++ = pG._HamDist(CurHaplo.Num_SNP, *p1, *p2);
						if (d < MinDiff) MinDiff = d;
						if (d == 0)
							HP.PairList.push_back(THaploPair(p1, p2));
					}
				}

				if (MinDiff > 0)
				{
					int *pD = &DiffList[0];
					p1 = &NextHaplo.List[pH1_st];
					for (size_t n1=pH1_n; n1 > 0; n1--, p1++)
					{
						p2 = p1;
						for (size_t n2=n1; n2 > 0; n2--, p2++)
						{
							if (*pD++ == MinDiff)
								HP.PairList.push_back(THaploPair(p1, p2));
						}
					}
				}
			}
		}
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
	int allele_cnt = 0, valid_cnt = 0;
	for (int iSamp=0; iSamp < SNPMat.Num_Total_Samp; iSamp++)
	{
		int dup = GenoList.List[iSamp].BootstrapCount;
		if (dup > 0)
		{
			int g = SNPMat.Get(iSamp, NewSNP);
			if ((0<=g) && (g<=2))
				{ allele_cnt += g*dup; valid_cnt += 2*dup; }
		}
	}
	if ((allele_cnt==0) || (allele_cnt==valid_cnt)) return false;

	// initialize the haplotype frequencies in NextHaplo
	CurHaplo.DoubleHaplosInitFreq(NextHaplo, double(allele_cnt)/valid_cnt);

	// update haplotype pair
	const int IdxNewSNP = NextHaplo.Num_SNP - 1;
	vector<THaploPairList>::iterator it;

	for (it = _SampHaploPair.begin(); it != _SampHaploPair.end(); it++)
	{
		vector<THaploPair>::iterator p;
		// SNP genotype
		int geno = SNPMat.Get(it->SampIndex, NewSNP);
		if ((0<=geno) && (geno<=2))
		{
			// for -- loop
			for (p = it->PairList.begin(); p != it->PairList.end(); p++)
			{
				p->Flag = ((p->H1->GetAllele(IdxNewSNP) +
					p->H2->GetAllele(IdxNewSNP)) == geno);
			}
		} else {
			// for -- loop
			for (p = it->PairList.begin(); p != it->PairList.end(); p++)
				p->Flag = true;
		}
	}

	return true;
}

void CAlg_EM::ExpectationMaximization(CHaplotypeList &NextHaplo)
{
	HIBAG_TIMING(TM_EM_ALG)

	// the converage tolerance
	double ConvTol = 0, LogLik = -1e+30;

	// iterate ...
	for (int iter=0; iter <= EM_MaxNum_Iterations; iter++)
	{
		// save old values
		// old log likelihood
		double Old_LogLik = LogLik;
		// save old haplotype frequencies
		NextHaplo.SaveClearFrequency();

		// for-loop each sample
		vector<THaploPairList>::iterator s;
		vector<THaploPair>::iterator p;
		int TotalNumSamp = 0;
		LogLik = 0;

		for (s = _SampHaploPair.begin(); s != _SampHaploPair.end(); s++)
		{
			// always "s->BootstrapCount > 0"
			TotalNumSamp += s->BootstrapCount;

			double psum = 0;
			for (p = s->PairList.begin(); p != s->PairList.end(); p++)
			{
				if (p->Flag)
				{
					p->GenoFreq = (p->H1 != p->H2) ?
						(2 * p->H1->aux.OldFreq * p->H2->aux.OldFreq) :
						(p->H1->aux.OldFreq * p->H2->aux.OldFreq);
					psum += p->GenoFreq;
				}
			}
			LogLik += s->BootstrapCount * log(psum);
			psum = s->BootstrapCount / psum;

			// update
			for (p = s->PairList.begin(); p != s->PairList.end(); p++)
			{
				if (p->Flag)
				{
					double r = p->GenoFreq * psum;
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

CAlg_Prediction::CAlg_Prediction() { }

void CAlg_Prediction::InitPrediction(int n_hla)
{
	HIBAG_CHECKING(n_hla<=0, "CAlg_Prediction::Init, n_hla error.");

	_nHLA = n_hla;
	const int size = n_hla*(n_hla+1)/2;
	_PostProb.resize(size);
	_SumPostProb.resize(size);
}

void CAlg_Prediction::InitPostProbBuffer()
{
	memset(&_PostProb[0], 0, _PostProb.size()*sizeof(double));
}

void CAlg_Prediction::InitSumPostProbBuffer()
{
	memset(&_SumPostProb[0], 0, _SumPostProb.size()*sizeof(double));
	_Sum_Weight = 0;
}

void CAlg_Prediction::AddProbToSum(double weight)
{
	if (weight > 0)
	{
		double *p = &_PostProb[0];
		double *s = &_SumPostProb[0];
		for (size_t n = _SumPostProb.size(); n > 0; n--)
			(*s++) += (*p++) * weight;
		_Sum_Weight += weight;
	}
}

void CAlg_Prediction::NormalizeSumPostProb()
{
	if (_Sum_Weight > 0)
	{
		const double scale = 1.0 / _Sum_Weight;
		double *s = &_SumPostProb[0];
		for (size_t n = _SumPostProb.size(); n > 0; n--)
			*s++ *= scale;
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

void CAlg_Prediction::PredictPostProb(const CHaplotypeList &Haplo,
	const TGenotype &Geno, double &SumProb)
{
	THaplotype *I1, *I2;
	double *pProb = &_PostProb[0];
	double sum;

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
				sum += FREQ_MUTANT((i1 != i2) ?
					(2 * i1->Freq * i2->Freq) : (i1->Freq * i2->Freq),
					Geno._HamDist(Haplo.Num_SNP, *i1, *i2));
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
					sum += FREQ_MUTANT(2 * i1->Freq * i2->Freq,
						Geno._HamDist(Haplo.Num_SNP, *i1, *i2));
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

THLAType CAlg_Prediction::_PredBestGuess(const CHaplotypeList &Haplo,
	const TGenotype &Geno)
{
	THLAType rv;
	rv.Allele1 = rv.Allele2 = NA_INTEGER;
	double max=0, prob;

	THaplotype *I1, *I2;
	I1 = Haplo.List;

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
				prob += FREQ_MUTANT((i1 != i2) ?
					(2 * i1->Freq * i2->Freq) : (i1->Freq * i2->Freq),
					Geno._HamDist(Haplo.Num_SNP, *i1, *i2));
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
					prob += FREQ_MUTANT(2 * i1->Freq * i2->Freq,
						Geno._HamDist(Haplo.Num_SNP, *i1, *i2));
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

double CAlg_Prediction::_PredPostProb(const CHaplotypeList &Haplo,
	const TGenotype &Geno, const THLAType &HLA)
{
	int H1=HLA.Allele1, H2=HLA.Allele2;
	if (H1 > H2) std::swap(H1, H2);
	int IxHLA = H2 + H1*(2*_nHLA-H1-1)/2;
	int idx = 0;

	double sum=0, hlaProb=0, prob;
	THaplotype *I1, *I2;
	I1 = Haplo.List;

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
				prob += FREQ_MUTANT((i1 != i2) ?
					(2 * i1->Freq * i2->Freq) : (i1->Freq * i2->Freq),
					Geno._HamDist(Haplo.Num_SNP, *i1, *i2));
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
					prob += FREQ_MUTANT(2 * i1->Freq * i2->Freq,
						Geno._HamDist(Haplo.Num_SNP, *i1, *i2));
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

THLAType CAlg_Prediction::BestGuess()
{
	THLAType rv;
	rv.Allele1 = rv.Allele2 = NA_INTEGER;

	double *p = &_PostProb[0];
	double max = 0;
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

THLAType CAlg_Prediction::BestGuessEnsemble()
{
	THLAType rv;
	rv.Allele1 = rv.Allele2 = NA_INTEGER;

	double *p = &_SumPostProb[0];
	double max = 0;
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



// -------------------------------------------------------------------------
// The algorithm of variable selection

CVariableSelection::CVariableSelection()
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
	for (int i=0; i < nSamp(); i++)
	{
		int cnt = _GenoList.List[i].BootstrapCount;
		tmp[_HLAList->List[i].Allele1] += cnt;
		tmp[_HLAList->List[i].Allele2] += cnt;
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
	if (GPUExtProcPtr)
	{
		Haplo.SetHaploAux();
		(*GPUExtProcPtr->build_set_haplo_geno)(Haplo.List, Haplo.Num_Haplo,
			&Geno.List[0], Haplo.Num_SNP);
	}
}

void CVariableSelection::_Done_EvalAcc()
{

}

int CVariableSelection::_OutOfBagAccuracy(CHaplotypeList &Haplo)
{
	HIBAG_TIMING(TM_ACC_OOB)
	HIBAG_CHECKING(Haplo.Num_SNP != _GenoList.Num_SNP,
		"CVariableSelection::_OutOfBagAccuracy, Haplo and GenoList should have the same number of SNP markers.");

	int CorrectCnt=0;
	if (GPUExtProcPtr)
	{
		CorrectCnt = (*GPUExtProcPtr->build_acc_oob)();
	} else {
		vector<TGenotype>::const_iterator p = _GenoList.List.begin();
		for (; p != _GenoList.List.end(); p++)
		{
			if (p->BootstrapCount <= 0)
			{
				THLAType g = _Predict._PredBestGuess(Haplo, *p);
				CorrectCnt += CHLATypeList::Compare(g, p->aux_hla_type);
			}
		}
	}

	return CorrectCnt;
}

double CVariableSelection::_InBagLogLik(CHaplotypeList &Haplo)
{
	HIBAG_TIMING(TM_ACC_IB)
	HIBAG_CHECKING(Haplo.Num_SNP != _GenoList.Num_SNP,
		"CVariableSelection::_InBagLogLik, Haplo and GenoList should have the same number of SNP markers.");

	double LogLik = 0;
	if (GPUExtProcPtr)
	{
		LogLik = (*GPUExtProcPtr->build_acc_ib)();
	} else {
		vector<TGenotype>::const_iterator p = _GenoList.List.begin();
		for (; p != _GenoList.List.end(); p++)
		{
			if (p->BootstrapCount > 0)
			{
				LogLik += p->BootstrapCount *
					log(_Predict._PredPostProb(Haplo, *p, p->aux_hla_type));
			}
		}
		LogLik *= -2;
	}
	return LogLik;
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
	int Global_Max_OutOfBagAcc = 0;  // # of correct alleles
	double Global_Min_Loss = 1e+30;
	int NumOOB = 0;
	{
		vector<TGenotype>::const_iterator p = _GenoList.List.begin();
		for (; p != _GenoList.List.end(); p++)
		{
			if (p->BootstrapCount <= 0) NumOOB ++;
		}
		if (NumOOB <= 0) NumOOB = 1;
	}

	// reserve memory for haplotype lists
	const size_t reserve_num_haplo = nSamp() * 2;
	CHaplotypeList NextHaplo(reserve_num_haplo);
	CHaplotypeList NextReducedHaplo(reserve_num_haplo);
	CHaplotypeList MinHaplo(reserve_num_haplo);

	while (VarSampling.TotalNum() > 0 &&
		OutSNPIndex.size() < HIBAG_MAXNUM_SNP_IN_CLASSIFIER)
	{
		// prepare for growing the individual classifier
		_EM.PrepareHaplotypes(OutHaplo, _GenoList, *_HLAList, NextHaplo);

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
				 Rprintf("    %2d, SNP: %d, Loss: %g, OOB Acc: %0.2f%%, # of Haplo: %d\n",
					OutSNPIndex.size(), OutSNPIndex.back()+1,
					Global_Min_Loss,
					double(Global_Max_OutOfBagAcc) / NumOOB * 50,
					OutHaplo.Num_Haplo);
			}
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

void CAttrBag_Classifier::InitBootstrapCount(int SampCnt[])
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

void CAttrBag_Model::BuildClassifiers(int nclassifier, int mtry, bool prune,
	bool verbose, bool verbose_detail)
{
#ifdef HIBAG_ENABLE_TIMING
	memset(timing_array, 0, sizeof(timing_array));
	HIBAG_TIMING(TM_TOTAL)
#endif

	if (verbose)
		Rprintf("[-] %s\n", date_text());

	if (GPUExtProcPtr)
		(*GPUExtProcPtr->build_init)(nHLA(), nSamp());

	CSamplingWithoutReplace VarSampling;

	for (int k=0; k < nclassifier; k++)
	{
		VarSampling.Init(nSNP());

		CAttrBag_Classifier *I = NewClassifierBootstrap();
		if (GPUExtProcPtr)
			(*GPUExtProcPtr->build_set_bootstrap)(&(I->BootstrapCount()[0]));

		I->Grow(VarSampling, mtry, prune, verbose, verbose_detail);
		if (verbose)
		{
			Rprintf(
				"[%d] %s, OOB Acc: %0.2f%%, # of SNPs: %d, # of Haplo: %d\n",
				k+1, date_text(), I->OutOfBag_Accuracy()*100, I->nSNP(),
				I->nHaplo());
		}
	}

	if (GPUExtProcPtr)
		(*GPUExtProcPtr->build_done)();

#ifdef HIBAG_ENABLE_TIMING
	tm.Stop();
	Rprintf("It took %0.2f seconds in total:\n"
			"    _OutOfBagAccuracy(): %0.2f%%, %0.2fs\n"
			"    _InBagLogLik(): %0.2f%%, %0.2fs\n"
			"    PrepareHaplotypes(): %0.2f%%, %0.2fs\n"
			"    ExpectationMaximization(): %0.2f%%, %0.2fs\n",
		((double)timing_array[TM_TOTAL]) / CLOCKS_PER_SEC,
		100.0 * timing_array[TM_ACC_OOB] / timing_array[TM_TOTAL],
		((double)timing_array[TM_ACC_OOB]) / CLOCKS_PER_SEC,
		100.0 * timing_array[TM_ACC_IB] / timing_array[TM_TOTAL],
		((double)timing_array[TM_ACC_IB]) / CLOCKS_PER_SEC,
		100.0 * timing_array[TM_PRE_HAPLO] / timing_array[TM_TOTAL],
		((double)timing_array[TM_PRE_HAPLO]) / CLOCKS_PER_SEC,
		100.0 * timing_array[TM_EM_ALG] / timing_array[TM_TOTAL],
		((double)timing_array[TM_EM_ALG]) / CLOCKS_PER_SEC
	);
#endif
}

void CAttrBag_Model::PredictHLA(const int *genomat, int n_samp, int vote_method,
	int OutH1[], int OutH2[], double OutMaxProb[], double OutMatching[],
	double OutProbArray[], bool ShowInfo)
{
	if ((vote_method < 1) || (vote_method > 2))
		throw ErrHLA("Invalid 'vote_method'.");

	_Predict.InitPrediction(nHLA());
	Progress.Info = "Predicting";
	Progress.Init(n_samp, ShowInfo);

	vector<int> snp_weight(nSNP());
	_GetSNPWeights(&snp_weight[0]);
	const size_t nn = nHLA()*(nHLA()+1)/2;

	_Init_PredictHLA();
	for (int i=0; i < n_samp; i++, genomat+=nSNP())
	{
		double pb;
		_PredictHLA(genomat, &snp_weight[0], vote_method, pb);

		THLAType HLA = _Predict.BestGuessEnsemble();
		OutH1[i] = HLA.Allele1; OutH2[i] = HLA.Allele2;

		if ((HLA.Allele1 != NA_INTEGER) && (HLA.Allele2 != NA_INTEGER))
			OutMaxProb[i] = _Predict.IndexSumPostProb(HLA.Allele1, HLA.Allele2);
		else
			OutMaxProb[i] = 0;

		if (OutProbArray)
		{
			memcpy(OutProbArray, &_Predict.SumPostProb()[0], sizeof(double)*nn);
			OutProbArray += nn;
		}
		if (OutMatching) OutMatching[i] = pb;

		Progress.Forward(1, ShowInfo);
	}
	_Done_PredictHLA();
}

void CAttrBag_Model::_PredictHLA(const int geno[], const int snp_weight[],
	int vote_method, double &OutMatching)
{
	// weight for each classifier, based on missing proportion
	double weight[_ClassifierList.size()];
	vector<CAttrBag_Classifier>::const_iterator p = _ClassifierList.begin();
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
		weight[w_i++] = double(nw) / sum;
	}

	if (GPUExtProcPtr)
	{
		p = _ClassifierList.begin();
		for (size_t w_i=0; p != _ClassifierList.end(); p++, w_i++)
		{
			gpu_geno_buf[w_i].IntToSNP(p->nSNP(), geno, &(p->_SNPIndex[0]));
		}
		(*GPUExtProcPtr->predict_avg_prob)(&gpu_geno_buf[0], weight,
			&_Predict._SumPostProb[0], &OutMatching);

	} else {
		// initialize probability
		_Predict.InitSumPostProbBuffer();
		TGenotype Geno;
		double sum_pb=0, pb;

		p = _ClassifierList.begin();
		for (size_t w_i=0; p != _ClassifierList.end(); p++, w_i++)
		{
			if (weight[w_i] <= 0) continue;

			Geno.IntToSNP(p->nSNP(), geno, &(p->_SNPIndex[0]));
			_Predict.PredictPostProb(p->_Haplo, Geno, pb);
			sum_pb += pb;

			if (vote_method == 1)
			{
				// predicting based on the averaged posterior probabilities
				_Predict.AddProbToSum(weight[w_i]);
			} else if (vote_method == 2)
			{
				// predicting by class majority voting
				THLAType pd = _Predict.BestGuess();
				if ((pd.Allele1 != NA_INTEGER) && (pd.Allele2 != NA_INTEGER))
				{
					_Predict.InitPostProbBuffer();  // fill by ZERO
					_Predict.IndexPostProb(pd.Allele1, pd.Allele2) = 1.0;
					_Predict.AddProbToSum(1.0);
				}
			}
		}

		// normalize the sum of posterior prob
		_Predict.NormalizeSumPostProb();
		OutMatching = sum_pb / _ClassifierList.size();
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

void CAttrBag_Model::_Init_PredictHLA()
{
	if (GPUExtProcPtr)
	{
		// prepare data structure for GPU
		const size_t n_classifier = _ClassifierList.size();
		THaplotype* haplo[n_classifier];

		gpu_geno_buf.resize(n_classifier);
		gpu_num_haplo.resize(n_classifier*2);
		int *pg = &gpu_num_haplo[0];

		vector<CAttrBag_Classifier>::iterator p;
		p = _ClassifierList.begin();
		for (size_t c_i=0; p != _ClassifierList.end(); p++, c_i++)
		{
			CHaplotypeList &hl = p->_Haplo;
			hl.SetHaploAux();
			haplo[c_i] = hl.List;
			*pg++ = p->nHaplo();
			*pg++ = p->nSNP();
		}

		(*GPUExtProcPtr->predict_init)(nHLA(), n_classifier, haplo,
			&gpu_num_haplo[0]);
	}
}

void CAttrBag_Model::_Done_PredictHLA()
{
	if (GPUExtProcPtr)
	{
		(*GPUExtProcPtr->predict_done)();
	}
}
