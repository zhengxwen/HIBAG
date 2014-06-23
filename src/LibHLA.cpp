// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// Copyright (C) 2013	Xiuwen Zheng (zhengx@u.washington.edu)
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

// ===============================================================
// Name           : LibHLA
// Author         : Xiuwen Zheng
// Version        : 1.2.1
// Copyright      : Xiuwen Zheng (GPL v3.0)
// Created        : 11/14/2011
// Last modified  : 10/29/2013
// Description    : HLA Genotype Imputation with Attribute Bagging
// ===============================================================


#if defined(unix) || defined(__unix) || defined(__unix__) || defined(macintosh) || defined(__APPLE__) || defined(__APPLE_CC__)
#  define HIBAG_SYS_UNIX
#elif defined(_WIN32) || defined(__WIN32__) || defined(WIN32) || defined(_WIN32_WCE)
#  define HIBAG_SYS_WIN
#endif


#include <StructHLA.h>

#ifdef HIBAG_ALLOW_GPU_SUPPORT
#  if defined(HIBAG_SYS_UNIX)
#    include <dlfcn.h>
#  elif defined(HIBAG_SYS_WIN)
#    include <windows.h>
#  else
#    error "not supported"
#  endif
#endif

#include <LibHLA.h>

#define HIBAG_TIMING	0
// 0: No timing
// 1: Spends ~90% of time on 'CVariableSelection::_OutOfBagAccuracy'
//    and 'CVariableSelection::_InBagLogLik'
//    Could be improved on speed by GPU programming
// 2: ~9.2% of time on CAlg_EM::ExpectationMaximization
// 3: ~0.5% of time on CAlg_EM::PrepareHaplotypes

#if (HIBAG_TIMING > 0)
#  include <time.h>
#endif


using namespace std;
using namespace HLA_LIB;


// ************************************************************************* //
// ************************************************************************* //

// parameters -- EM algorithm
/// The max number of iterations
int HLA_LIB::EM_MaxNum_Iterations = 500;
/// The initial value of EM algorithm
static const TFLOAT EM_INIT_VAL_FRAC = 0.001;
/// The reltol convergence tolerance, sqrt(machine.epsilon) by default, used in EM algorithm
TFLOAT HLA_LIB::EM_FuncRelTol = sqrt(FLOAT_EPSILON);

// parameters -- reduce the number of possible haplotypes
/// The minimum rare frequency to store haplotypes
static const TFLOAT MIN_RARE_FREQ = 1e-5;
/// The fraction of one haplotype that can be ignored
static const TFLOAT FRACTION_HAPLO = 1.0/10;



// the parameter of searching SNP markers

/// The reltol for the stopping rule of adding a new SNP marker
static const TFLOAT STOP_RELTOL_LOGLIK_ADDSNP = 0.001;
/// The reltol for erasing the SNP marker is prune = TRUE
static const TFLOAT PRUNE_RELTOL_LOGLIK = 0.1;



#ifdef HIBAG_ALLOW_GPU_SUPPORT

// whether run in double precision, depending on the device
static bool _gpuRunInDoublePrecision = false;

typedef void (*TGPUFunc)();
typedef void (*TGPUInitFunc)(UINT8 *hamdist);
typedef int (*TGPU_OutOfBagAcc_F32)(int nHLA, int *_HLA_HapIdx,
	int nHaplo, TGPU_Haplotype_F32 *_HapList,
	int nGeno, TGPU_Genotype *_GList, int nSNP);
typedef int (*TGPU_OutOfBagAcc_F64)(int nHLA, int *_HLA_HapIdx,
	int nHaplo, TGPU_Haplotype_F64 *_HapList,
	int nGeno, TGPU_Genotype *_GList, int nSNP);

static TGPUInitFunc _gpuInitialize = NULL;
static TGPUFunc     _gpuFinalize = NULL;
static TGPUFunc     _gpuDeviceQuery = NULL;
static TGPU_OutOfBagAcc_F32  _gpuOutOfBagAcc_F32 = NULL;
static TGPU_OutOfBagAcc_F64  _gpuOutOfBagAcc_F64 = NULL;

#endif




// ************************************************************************* //
// ************************************************************************* //

// return 0 .. (Range-1)
static inline int RandomNum(int Range)
{
	int rv = (int)(unif_rand()*(Range-1) + 0.5);
	if (rv >= Range) rv = Range -1;
	return rv;
}



// ************************************************************************* //

/// Frequency Calculation
#define FREQ_MUTANT(p, cnt)    ((p) * EXP_LOG_MIN_RARE_FREQ[cnt]);

/// exp(cnt * log(MIN_RARE_FREQ)), cnt is the hamming distance
static TFLOAT EXP_LOG_MIN_RARE_FREQ[HIBAG_MAXNUM_SNP_IN_CLASSIFIER*2];

/// the hamming distance
/// SNP genotypes (first dimension) and a pair of haplotypes (second dimension)
static UINT8 _Packed_HammingDistance[256][256] HIBAG_SSE_VAR_ALIGN;

#ifdef HIBAG_SSE_OPTIMIZE_HAMMING_DISTANCE
// 32 -- the length of genotypes
static __m128 _Packed_HamLenMask[32] HIBAG_SSE_VAR_ALIGN;
#endif

class CInit
{
public:
	CInit()
	{
		for (int i=0; i < HIBAG_MAXNUM_SNP_IN_CLASSIFIER*2; i++)
			EXP_LOG_MIN_RARE_FREQ[i] = FLOAT_EXP(i * FLOAT_LOG(MIN_RARE_FREQ));
		EXP_LOG_MIN_RARE_FREQ[0] = 1;
		for (int i=0; i < HIBAG_MAXNUM_SNP_IN_CLASSIFIER*2; i++)
		{
			if (!R_finite(EXP_LOG_MIN_RARE_FREQ[i]))
				EXP_LOG_MIN_RARE_FREQ[i] = 0;
		}

		// the hamming distance
		for (int I1=0; I1 < 256; I1++)
		{
			int G[4] = { I1 & 0x03, (I1 >> 2) & 0x03,
				(I1 >> 4) & 0x03, (I1 >> 6) & 0x03 };
			for (int I2=0; I2 < 256; I2++)
			{
				int H1[4] = { (I2 >> 0) & 0x01, (I2 >> 2) & 0x01,
					(I2 >> 4) & 0x01, (I2 >> 6) & 0x01 };
				int H2[4] = { (I2 >> 1) & 0x01, (I2 >> 3) & 0x01,
					(I2 >> 5) & 0x01, (I2 >> 7) & 0x01 };
				int sum = 0;
				for (int k=0; k < 4; k++)
				{
					if (G[k] <= 2)
						sum += abs(G[k] - H1[k] - H2[k]);
				}
				_Packed_HammingDistance[I1][I2] = sum;
			}
		}

	#ifdef HIBAG_SSE_OPTIMIZE_HAMMING_DISTANCE

		// the hamming mask for length
		memset(_Packed_HamLenMask, 0, sizeof(_Packed_HamLenMask));

		for (int Length=1; Length < 32; Length++)
		{
			UINT8 *p = (UINT8*)(&_Packed_HamLenMask[Length]);
			int L = Length;
			for (; L >= 4; L-=4, p+=2)
				p[0] = p[1] = 0xFF;
			if (L > 0)
				p[0] = p[1] = ~(0xFF << (L*2));
		}

	#endif
	}
};

static CInit _Init;


// ************************************************************************* //

#if (HIBAG_TIMING > 0)

static clock_t _timing_ = 0;
static clock_t _timing_last_point;

static inline void _put_timing()
{
	_timing_last_point = clock();
}
static inline void _inc_timing()
{
	clock_t t = clock();
	_timing_ += t - _timing_last_point;
	_timing_last_point = t;
}

#endif




// ************************************************************************* //
// ************************************************************************* //

// CdProgression

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
	int p = int(TFLOAT(TotalPercent)*fCurrent / fTotal);
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
	time_t tm; time(&tm);
	string s(ctime(&tm));
	s.erase(s.size()-1, 1);
	if (Info.empty())
	{
		Rprintf("%s\t%d%%\n", s.c_str(), int(fPercent*StepPercent));
	} else {
		Rprintf("%s\t%s\t%d%%\n", Info.c_str(), s.c_str(),
			int(fPercent*StepPercent));
	}
}


/// The progression information
CdProgression HLA_LIB::Progress;



// ************************************************************************* //
// ************************************************************************* //

// -------------------------------------------------------------------------
// The class of haplotype structure

THaplotype::THaplotype()
{
	memset(PackedHaplo, 0, sizeof(PackedHaplo));
}

THaplotype::THaplotype(const TFLOAT _freq)
{
	memset(PackedHaplo, 0, sizeof(PackedHaplo));
	Frequency = _freq;
	OldFreq = 0;
}

THaplotype::THaplotype(const char *str, const TFLOAT _freq)
{
	memset(PackedHaplo, 0, sizeof(PackedHaplo));
	Frequency = _freq;
	OldFreq = 0;
	StrToHaplo(str);
}

UINT8 THaplotype::GetAllele(const int idx) const
{
	HIBAG_CHECKING((idx<0) || (idx>=HIBAG_MAXNUM_SNP_IN_CLASSIFIER),
		"THaplotype::GetAllele, invalid index.");
	UINT8 ch = PackedHaplo[idx >> 2];
	return (ch >> (2*(idx & 0x03))) & 0x01;
}

void THaplotype::SetAllele(const int idx, UINT8 val)
{
	HIBAG_CHECKING((idx<0) || (idx>=HIBAG_MAXNUM_SNP_IN_CLASSIFIER),
		"THaplotype::SetAllele, invalid index.");
	HIBAG_CHECKING(val!=0 && val!=1,
		"THaplotype::SetAllele, the value should be 0 or 1.");
	UTYPE &ch = PackedHaplo[idx >> 2];
	UINT8 shift = (idx & 0x03)*2;
	UTYPE mask = ~(0x03 << shift);
	ch = (ch & mask) | (val << shift);
}

string THaplotype::HaploToStr(const int Length) const
{
	HIBAG_CHECKING((Length<0) || (Length>HIBAG_MAXNUM_SNP_IN_CLASSIFIER),
		"THaplotype::HaploToStr, the length is invalid.");
	string rv;
	if (Length > 0)
	{
		rv.resize(Length);
		for (int i=0; i < Length; i++)
			rv[i] = (GetAllele(i)==0) ? '0' : '1';
	}
	return rv;
}

void THaplotype::StrToHaplo(const string &str)
{
	HIBAG_CHECKING((int)str.size() > HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"THaplotype::StrToHaplo, the input string is too long.");
	for (int i=0; i < (int)str.size(); i++)
	{
		char ch = str[i];
		HIBAG_CHECKING(ch!='0' && ch!='1',
			"THaplotype::StrToHaplo, the input string should be '0' or '1'");
		SetAllele(i, ch-'0');
	}
}



// -------------------------------------------------------------------------
// The class of genotype structure

TGenotype::TGenotype()
{
	memset(PackedSNPs, 0, sizeof(PackedSNPs));
}

UINT8 TGenotype::GetSNP(const int idx) const
{
	HIBAG_CHECKING((idx<0) || (idx>=HIBAG_MAXNUM_SNP_IN_CLASSIFIER),
		"TGenotype::GetSNP, invalid index.");
	UINT8 ch = PackedSNPs[idx >> 2];
	return (ch >> (2*(idx & 0x03))) & 0x03;
}

void TGenotype::SetSNP(const int idx, int val)
{
	HIBAG_CHECKING((idx<0) || (idx>=HIBAG_MAXNUM_SNP_IN_CLASSIFIER),
		"TGenotype::SetSNP, invalid index.");
	if (val<0 || val>2) val = 3;
	_SetSNP(idx, val);
}

void TGenotype::_SetSNP(const int idx, UINT8 val)
{
	UTYPE &ch = PackedSNPs[idx >> 2];
	UINT8 shift = (idx & 0x03)*2;
	UTYPE mask = ~(0x03 << shift);
	ch = (ch & mask) | (val << shift);
}

string TGenotype::SNPToString(const int Length) const
{
	HIBAG_CHECKING(Length > HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"TGenotype::SNPToString, the length is too large.");
	string rv;
	if (Length > 0)
	{
		rv.resize(Length);
		for (int i=0; i < Length; i++)
		{
			UINT8 ch = GetSNP(i);
			rv[i] = (ch < 3) ? (ch + '0') : '?';
		}
	}
	return rv;
}

void TGenotype::StringToSNP(const string &str)
{
	HIBAG_CHECKING((int)str.size() > HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"TGenotype::StringToSNP, the input string is too long.");
	for (int i=0; i < (int)str.size(); i++)
	{
		char ch = str[i];
		HIBAG_CHECKING(ch!='0' && ch!='1' && ch!='2' && ch!='?',
			"TGenotype::StringToSNP, the input string should be '0', '1', '2' or '?'.");
		SetSNP(i, ch-'0');
	}
}

void TGenotype::SNPToInt(const int Length, int OutArray[]) const
{
	HIBAG_CHECKING((Length<0) || (Length>HIBAG_MAXNUM_SNP_IN_CLASSIFIER),
		"TGenotype::SNPToInt, the length is invalid.");
	for (int i=0; i < Length; i++) OutArray[i] = GetSNP(i);
}

void TGenotype::IntToSNP(int Length, const int InBase[], const int Index[])
{
	HIBAG_CHECKING((Length<0) || (Length > HIBAG_MAXNUM_SNP_IN_CLASSIFIER),
		"TGenotype::IntToSNP, the length is invalid.");
	UTYPE *p = &PackedSNPs[0];
	for (; Length >= 4; Length-=4, Index+=4)
	{
		int g0 = InBase[Index[0]]; if (g0<0 || g0>2) g0 = 3;
		int g1 = InBase[Index[1]]; if (g1<0 || g1>2) g1 = 3;
		int g2 = InBase[Index[2]]; if (g2<0 || g2>2) g2 = 3;
		int g3 = InBase[Index[3]]; if (g3<0 || g3>2) g3 = 3;
		*p++ = g0 | (g1 << 2) | (g2 << 4) | (g3 << 6);
	}
	if (Length > 0)
	{
		*p = 0;
		for (int i=0; i < Length; i++)
		{
			int g = InBase[*Index++]; if (g<0 || g>2) g = 3;
			*p |= (g << (2*i));
		}
	}
}

int TGenotype::HammingDistance(int Length,
	const THaplotype &H1, const THaplotype &H2) const
{
	HIBAG_CHECKING((Length<0) || (Length > HIBAG_MAXNUM_SNP_IN_CLASSIFIER),
		"THaplotype::HammingDistance, the length is too large.");
	return _HamDist(Length, H1, H2);
}

inline int TGenotype::_HamDist(int Length,
	const THaplotype &H1, const THaplotype &H2) const
{
	const UTYPE *p1 = &H1.PackedHaplo[0];
	const UTYPE *p2 = &H2.PackedHaplo[0];
	const UTYPE *s = &PackedSNPs[0];
	int rv = 0;

#ifdef HIBAG_SSE_OPTIMIZE_HAMMING_DISTANCE

	// UTYPE = UINT16
	static UINT16 offset[8] HIBAG_SSE_VAR_ALIGN;

	while (Length > 0)
	{
		// two haplotypes
		__m128i h1 = _mm_loadu_si128((__m128i*)p1); p1 += 8;
		__m128i h2 = _mm_loadu_si128((__m128i*)p2); p2 += 8;
		// h1 = h1 | (h2 << 1)
		h1 = (__m128i)_mm_or_ps((__m128)h1, (__m128)_mm_slli_epi16(h2, 1));

		// SNP genotypes
		__m128i g = _mm_loadu_si128((__m128i*)s); s += 8;
		g = (__m128i)_mm_or_ps((__m128)_mm_slli_si128(g, 1), (__m128)h1);


		#define HAM_BASE	(&_Packed_HammingDistance[0][0])

		// consider Length
		if (Length >= 32)
		{
			// store in offset
			_mm_store_si128((__m128i*)offset, g);

			// sum
			rv += HAM_BASE[offset[0]] + HAM_BASE[offset[1]] +
				HAM_BASE[offset[2]] + HAM_BASE[offset[3]] +
				HAM_BASE[offset[4]] + HAM_BASE[offset[5]] +
				HAM_BASE[offset[6]] + HAM_BASE[offset[7]];

			// decrease Length
			Length -= 32;

		} else {
			// set mask
			g = (__m128i)_mm_and_ps((__m128)g, _Packed_HamLenMask[Length]);
			// stored in offset
			_mm_store_si128((__m128i*)offset, g);

			// sum
			const UINT16 *pl = offset;
			if (Length >= 16)  // 4 bytes
			{
				rv += HAM_BASE[offset[0]] + HAM_BASE[offset[1]] +
					HAM_BASE[offset[2]] + HAM_BASE[offset[3]];
				Length -= 16; pl += 4;
			}
			if (Length >= 8)  // 2 bytes
			{
				rv += HAM_BASE[pl[0]] + HAM_BASE[pl[1]];
				Length -= 8; pl += 2;
			}
			if (Length >= 4)
			{
				rv += HAM_BASE[pl[0]];
				Length -= 4; pl ++;
			}
			if (Length > 0)
			{
				rv += HAM_BASE[pl[0]];
				Length = 0;
			}
		}

		#undef HAM_BASE
	}

#else

	// UTYPE = UINT8
	// potential vectorization for compiler
	for (; Length >= 16; Length -= 16)  // four bytes
	{
		rv += _Packed_HammingDistance[*s++][(*p1++) | (*p2++ << 1)];
		rv += _Packed_HammingDistance[*s++][(*p1++) | (*p2++ << 1)];
		rv += _Packed_HammingDistance[*s++][(*p1++) | (*p2++ << 1)];
		rv += _Packed_HammingDistance[*s++][(*p1++) | (*p2++ << 1)];
	}
	for (; Length >= 4; Length -= 4)  // one byte
	{
		rv += _Packed_HammingDistance[*s++][(*p1++) | (*p2++ << 1)];
	}
	if (Length > 0)
	{
		static const UINT8 MaskAry[4] = { 0x00, 0x03, 0x0F, 0x3F };
		UINT8 mask = MaskAry[Length];
		rv += _Packed_HammingDistance[*s & mask][((*p1) | (*p2 << 1)) & mask];
	}

#endif

	return rv;
}



// -------------------------------------------------------------------------
// The class of haplotype list

CHaplotypeList::CHaplotypeList()
{
	Num_SNP = 0;
}

void CHaplotypeList::DoubleHaplos(CHaplotypeList &OutHaplos) const
{
	HIBAG_CHECKING(Num_SNP >= HIBAG_MAXNUM_SNP_IN_CLASSIFIER,
		"CHaplotypeList::DoubleHaplos, there are too many SNP markers.");

	OutHaplos.Num_SNP = Num_SNP + 1;
	OutHaplos.List.resize(List.size());

	for (int i=0; i < (int)List.size(); i++)
	{
		const vector<THaplotype> &src = List[i];
		vector<THaplotype> &dst = OutHaplos.List[i];
		
		dst.resize(src.size()*2);
		for (int j=0; j < (int)src.size(); j++)
		{
			dst[2*j+0] = src[j];
			dst[2*j+0].SetAllele(Num_SNP, 0);
			dst[2*j+1] = src[j];
			dst[2*j+1].SetAllele(Num_SNP, 1);
		}
	}
}

void CHaplotypeList::DoubleHaplosInitFreq(CHaplotypeList &OutHaplos,
	const TFLOAT AFreq) const
{
	static const char *msg =
		"CHaplotypeList::DoubleHaplosInitFreq, the total number of haplotypes is not correct.";
	HIBAG_CHECKING(List.size() != OutHaplos.List.size(), msg);

	const TFLOAT p0 = 1-AFreq, p1 = AFreq;
	for (int i=0; i < (int)List.size(); i++)
	{
		const vector<THaplotype> &src = List[i];
		vector<THaplotype> &dst = OutHaplos.List[i];
		HIBAG_CHECKING(dst.size() != src.size()*2, msg);

		for (int j=0; j < (int)src.size(); j++)
		{
			dst[2*j+0].Frequency = src[j].Frequency*p0 + EM_INIT_VAL_FRAC;
			dst[2*j+1].Frequency = src[j].Frequency*p1 + EM_INIT_VAL_FRAC;
		}
	}
}

void CHaplotypeList::MergeDoubleHaplos(const TFLOAT RareProb,
	CHaplotypeList &OutHaplos) const
{
	OutHaplos.Num_SNP = Num_SNP;
	OutHaplos.List.resize(List.size());

	for (int i=0; i < (int)List.size(); i++)
	{
		const vector<THaplotype> &src = List[i];
		vector<THaplotype> &dst = OutHaplos.List[i];
		dst.clear();
		dst.reserve(src.size());
		
		for (int j=0; j < (int)src.size(); j += 2)
		{
			const THaplotype &p0 = src[j+0];
			const THaplotype &p1 = src[j+1];

			if ((p0.Frequency < RareProb) || (p1.Frequency < RareProb))
			{
				if (p0.Frequency >= p1.Frequency)
					dst.push_back(p0);
				else
					dst.push_back(p1);
				dst.back().Frequency = p0.Frequency + p1.Frequency;
			} else {
				dst.push_back(p0); dst.push_back(p1);
			}
		}
	}
}

void CHaplotypeList::EraseDoubleHaplos(const TFLOAT RareProb,
	CHaplotypeList &OutHaplos) const
{
	OutHaplos.Num_SNP = Num_SNP;
	OutHaplos.List.resize(List.size());
	TFLOAT sum = 0;

	for (int i=0; i < (int)List.size(); i++)
	{
		const vector<THaplotype> &src = List[i];
		vector<THaplotype> &dst = OutHaplos.List[i];
		dst.clear();
		dst.reserve(src.size());
		
		for (int j=0; j < (int)src.size(); j += 2)
		{
			const THaplotype &p0 = src[j+0];
			const THaplotype &p1 = src[j+1];
			TFLOAT sumfreq = p0.Frequency + p1.Frequency;

			if ((p0.Frequency < RareProb) || (p1.Frequency < RareProb))
			{
				if (sumfreq >= MIN_RARE_FREQ)
				{
					if (p0.Frequency >= p1.Frequency)
						dst.push_back(p0);
					else
						dst.push_back(p1);
					dst.back().Frequency = sumfreq;
					sum += sumfreq;
				}
			} else {
				dst.push_back(p0); dst.push_back(p1);
				sum += sumfreq;
			}
		}
	}

	OutHaplos.ScaleFrequency(1/sum);
}

void CHaplotypeList::SaveClearFrequency()
{
	vector< vector<THaplotype> >::iterator it;
	for (it = List.begin(); it != List.end(); it++)
	{
		vector<THaplotype>::iterator p;
		for (p = it->begin(); p != it->end(); p++)
		{
			p->OldFreq = p->Frequency;
			p->Frequency = 0;
		}
	}
}

void CHaplotypeList::ScaleFrequency(const TFLOAT scale)
{
	vector< vector<THaplotype> >::iterator it;
	for (it = List.begin(); it != List.end(); it++)
	{
		vector<THaplotype>::iterator p;
		for (p = it->begin(); p != it->end(); p++)
		{
			p->Frequency *= scale;
		}
	}
}

int CHaplotypeList::TotalNumOfHaplo() const
{
	vector< vector<THaplotype> >::const_iterator it;
	int cnt = 0;
	for (it = List.begin(); it != List.end(); it++)
		cnt += it->size();
	return cnt;
}

void CHaplotypeList::Print()
{
	for (int i=0; i < (int)List.size(); i++)
	{
		vector<THaplotype> &L = List[i];
		vector<THaplotype>::const_iterator it;
		for (it=L.begin(); it != L.end(); it++)
		{
			Rprintf("%0.7f\tHLA: %4d, %s\n", it->Frequency, i,
				it->HaploToStr(Num_SNP).c_str());
		}
	}
}


#ifdef HIBAG_ALLOW_GPU_SUPPORT

void CHaplotypeList::InitGPUHostData(int *&_HLA_HapIdx, void *&_HapList)
{
	// HLA allele index
	_HLA_HapIdx = new int[nHLA()*2];
	int st = 0;
	for (int i=0; i < (int)List.size(); i++)
	{
		_HLA_HapIdx[2*i + 0] = st;
		_HLA_HapIdx[2*i + 1] = List[i].size();
		st += List[i].size();
	}


	// haplotype list
#if (HIBAG_FLOAT_TYPE_ID == 0)

	if (_gpuRunInDoublePrecision)
	{
		_HapList = new TGPU_Haplotype_F64[TotalNumOfHaplo()];
		TGPU_Haplotype_F64 *pH = (TGPU_Haplotype_F64 *)_HapList;

		vector< vector<THaplotype> >::iterator it;
		vector<THaplotype>::iterator p;
		for (it = List.begin(); it != List.end(); it++)
		{
			for (p = it->begin(); p != it->end(); p++)
			{
			#ifdef HIBAG_SSE_OPTIMIZE_HAMMING_DISTANCE
				// UTYPE = UINT16
				UINT8  *d = &(pH->PackedHaplo[0]);
				UINT16 *s = &(p->PackedHaplo[0]);
				for (size_t n=HIBAG_PACKED_UTYPE_MAXNUM_SNP; n > 0; n--)
					*d ++ = *s ++;
			#else
				memcpy(&(pH->PackedHaplo[0]), &(p->PackedHaplo[0]),
					HIBAG_PACKED_UTYPE_MAXNUM_SNP);
			#endif
				pH->Frequency = p->Frequency;
				pH ++;
			}
		}
	} else {
		_HapList = new TGPU_Haplotype_F32[TotalNumOfHaplo()];
		TGPU_Haplotype_F32 *pH = (TGPU_Haplotype_F32 *)_HapList;

		vector< vector<THaplotype> >::iterator it;
		vector<THaplotype>::iterator p;
		for (it = List.begin(); it != List.end(); it++)
		{
			for (p = it->begin(); p != it->end(); p++)
			{
			#ifdef HIBAG_SSE_OPTIMIZE_HAMMING_DISTANCE
				// UTYPE = UINT16
				UINT8  *d = &(pH->PackedHaplo[0]);
				UINT16 *s = &(p->PackedHaplo[0]);
				for (size_t n=HIBAG_PACKED_UTYPE_MAXNUM_SNP; n > 0; n--)
					*d ++ = *s ++;
			#else
				memcpy(&(pH->PackedHaplo[0]), &(p->PackedHaplo[0]),
					HIBAG_PACKED_UTYPE_MAXNUM_SNP);
			#endif
				pH->Frequency = p->Frequency;
				pH ++;
			}
		}
	}

#elif (HIBAG_FLOAT_TYPE_ID == 1)

	_HapList = new TGPU_Haplotype_F32[TotalNumOfHaplo()];
	TGPU_Haplotype_F32 *pH = (TGPU_Haplotype_F32 *)_HapList;

	vector< vector<THaplotype> >::iterator it;
	vector<THaplotype>::iterator p;
	for (it = List.begin(); it != List.end(); it++)
	{
		for (p = it->begin(); p != it->end(); p++)
		{
		#ifdef HIBAG_SSE_OPTIMIZE_HAMMING_DISTANCE
			// UTYPE = UINT16
			UINT8  *d = &(pH->PackedHaplo[0]);
			UINT16 *s = &(p->PackedHaplo[0]);
			for (size_t n=HIBAG_PACKED_UTYPE_MAXNUM_SNP; n > 0; n--)
				*d ++ = *s ++;
		#else
			memcpy(&(pH->PackedHaplo[0]), &(p->PackedHaplo[0]),
				HIBAG_PACKED_UTYPE_MAXNUM_SNP);
		#endif
			pH->Frequency = p->Frequency;
			pH ++;
		}
	}

#else
#  error "Invalid HIBAG_FLOAT_TYPE_ID"
#endif
}
#endif


#ifdef HIBAG_ALLOW_GPU_SUPPORT
inline void CHaplotypeList::FreqGPUHostData(int *&_HLA_HapIdx, void *&_HapList)
{
	// no worry if _HLA_HapIdx == NULL
	delete[] _HLA_HapIdx;
	_HLA_HapIdx = NULL;

#if (HIBAG_FLOAT_TYPE_ID == 0)

	if (_HapList)
	{
		if (_gpuRunInDoublePrecision)
		{
			TGPU_Haplotype_F64 *tmp = (TGPU_Haplotype_F64*)_HapList;
			delete[] tmp;
		} else {
			TGPU_Haplotype_F32 *tmp = (TGPU_Haplotype_F32*)_HapList;
			delete[] tmp;
		}
		_HapList = NULL;
	}

#elif (HIBAG_FLOAT_TYPE_ID == 1)

	if (_HapList)
	{
		TGPU_Haplotype_F32 *tmp = (TGPU_Haplotype_F32*)_HapList;
		delete[] tmp;
		_HapList = NULL;
	}

#else
#  error "Invalid HIBAG_FLOAT_TYPE_ID"
#endif
}
#endif



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
		if (g<0 || g>2) g = 3;
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

void CGenotypeList::Print()
{
	int idx = 0;
	vector<TGenotype>::const_iterator it;
	for (it = List.begin(); it != List.end(); it++)
	{
		idx ++;
		Rprintf("%d\t[%d]\t%s\n", idx, it->BootstrapCount,
			it->SNPToString(Num_SNP).c_str());
	}
}



// -------------------------------------------------------------------------
// A list of HLA types

CHLATypeList::CHLATypeList() { }

int CHLATypeList::Compare(const THLAType &H1, const THLAType &H2)
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
	for (int i=0; i < m_total; i++) _IdxArray[i] = i;
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
			int I = RandomNum(n_tmp - i);  // int(runif(0, 1)*(n_tmp - i - 1) + 0.5);
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
#if (HIBAG_TIMING == 3)
	_put_timing();
#endif

	HIBAG_CHECKING(GenoList.nSamp() != HLAList.nSamp(),
		"CAlg_EM::PrepareHaplotypes, GenoList and HLAList should have the same number of samples.");

	_SampHaploPair.clear();
	_SampHaploPair.reserve(GenoList.nSamp());
	CurHaplo.DoubleHaplos(NextHaplo);

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

			vector<THaplotype> &pH1 = NextHaplo.List[pHLA.Allele1];
			vector<THaplotype> &pH2 = NextHaplo.List[pHLA.Allele2];
			vector<THaplotype>::iterator p1, p2;
			int MinDiff = GenoList.Num_SNP * 4;

			if (pHLA.Allele1 != pHLA.Allele2)
			{
				for (p1 = pH1.begin(); p1 != pH1.end(); p1++)
				{
					for (p2 = pH2.begin(); p2 != pH2.end(); p2++)
					{
						int d = pG._HamDist(CurHaplo.Num_SNP, *p1, *p2);
						if (d < MinDiff) MinDiff = d;
						if (d == 0)
							HP.PairList.push_back(THaploPair(&(*p1), &(*p2)));
					}
				}

				if (MinDiff > 0)
				{
					for (p1 = pH1.begin(); p1 != pH1.end(); p1++)
					{
						for (p2 = pH2.begin(); p2 != pH2.end(); p2++)
						{
							if (pG._HamDist(CurHaplo.Num_SNP, *p1, *p2) == MinDiff)
								HP.PairList.push_back(THaploPair(&(*p1), &(*p2)));
						}
					}
				}
			} else {
				for (p1 = pH1.begin(); p1 != pH1.end(); p1++)
				{
					for (p2 = p1; p2 != pH1.end(); p2++)
					{
						int d = pG._HamDist(CurHaplo.Num_SNP, *p1, *p2);
						if (d < MinDiff) MinDiff = d;
						if (d == 0)
							HP.PairList.push_back(THaploPair(&(*p1), &(*p2)));
					}
				}

				if (MinDiff > 0)
				{
					for (p1 = pH1.begin(); p1 != pH1.end(); p1++)
					{
						for (p2 = p1; p2 != pH1.end(); p2++)
						{
							if (pG._HamDist(CurHaplo.Num_SNP, *p1, *p2) == MinDiff)
								HP.PairList.push_back(THaploPair(&(*p1), &(*p2)));
						}
					}
				}
			}
		}
	}

#if (HIBAG_TIMING == 3)
	_inc_timing();
#endif
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

	// initialize the haplotype frequencies
	CurHaplo.DoubleHaplosInitFreq(NextHaplo, TFLOAT(allele_cnt)/valid_cnt);

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
#if (HIBAG_TIMING == 2)
	_put_timing();
#endif

	// the converage tolerance
	TFLOAT ConvTol, LogLik = -1e+30;

	// iterate ...
	for (int iter=0; iter <= EM_MaxNum_Iterations; iter++)
	{
		// save old values
		// old log likelihood
		TFLOAT Old_LogLik = LogLik;
		// old haplotype frequencies
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

			TFLOAT psum = 0;
			for (p = s->PairList.begin(); p != s->PairList.end(); p++)
			{
				if (p->Flag)
				{
					p->Freq = (p->H1 != p->H2) ?
						(2 * p->H1->OldFreq * p->H2->OldFreq) : (p->H1->OldFreq * p->H2->OldFreq);
					psum += p->Freq;
				}
			}
			LogLik += s->BootstrapCount * FLOAT_LOG(psum);
			psum = TFLOAT(s->BootstrapCount) / psum;

			// update
			for (p = s->PairList.begin(); p != s->PairList.end(); p++)
			{
				if (p->Flag)
				{
					TFLOAT r = p->Freq * psum;
					p->H1->Frequency += r; p->H2->Frequency += r;
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

#if (HIBAG_TIMING == 2)
	_inc_timing();
#endif
}

void CAlg_EM::THaploPairList::Print(const int Length)
{
	vector<THaploPair>::const_iterator it;
	for (it = PairList.begin(); it != PairList.end(); it++)
	{
		Rprintf("%s, %s\n", it->H1->HaploToStr(Length).c_str(),
			it->H2->HaploToStr(Length).c_str());
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
	memset(&_PostProb[0], 0, _PostProb.size()*sizeof(TFLOAT));
}

void CAlg_Prediction::InitSumPostProbBuffer()
{
	memset(&_SumPostProb[0], 0, _SumPostProb.size()*sizeof(TFLOAT));
	_Sum_Weight = 0;
}

void CAlg_Prediction::AddProbToSum(const TFLOAT weight)
{
	if (weight > 0)
	{
		TFLOAT *p = &_PostProb[0];
		TFLOAT *s = &_SumPostProb[0];
		for (size_t n = _SumPostProb.size(); n > 0; n--, s++, p++)
			*s += (*p) * weight;
		_Sum_Weight += weight;
	}
}

void CAlg_Prediction::NormalizeSumPostProb()
{
	if (_Sum_Weight > 0)
	{
		const TFLOAT scale = 1.0 / _Sum_Weight;
		TFLOAT *s = &_SumPostProb[0];
		for (size_t n = _SumPostProb.size(); n > 0; n--)
			*s++ *= scale;
	}
}

TFLOAT &CAlg_Prediction::IndexPostProb(int H1, int H2)
{
	if (H1 > H2) std::swap(H1, H2);
	return _PostProb[H2 + H1*(2*_nHLA-H1-1)/2];
}

TFLOAT &CAlg_Prediction::IndexSumPostProb(int H1, int H2)
{
	if (H1 > H2) std::swap(H1, H2);
	return _SumPostProb[H2 + H1*(2*_nHLA-H1-1)/2];
}

void CAlg_Prediction::PredictPostProb(const CHaplotypeList &Haplo,
	const TGenotype &Geno)
{
	vector<THaplotype>::const_iterator i1;
	vector<THaplotype>::const_iterator i2;
	TFLOAT *pProb = &_PostProb[0];

	for (int h1=0; h1 < _nHLA; h1++)
	{
		const vector<THaplotype> &L1 = Haplo.List[h1];
		
		// diag value
		*pProb = 0;
		for (i1=L1.begin(); i1 != L1.end(); i1++)
		{
			for (i2=i1; i2 != L1.end(); i2++)
			{
				*pProb += FREQ_MUTANT((i1 != i2) ?
					(2 * i1->Frequency * i2->Frequency) : (i1->Frequency * i2->Frequency),
					Geno._HamDist(Haplo.Num_SNP, *i1, *i2));
			}
		}
		pProb ++;

		// off-diag value
		for (int h2=h1+1; h2 < _nHLA; h2++)
		{
			const vector<THaplotype> &L2 = Haplo.List[h2];
			*pProb = 0;
			for (i1=L1.begin(); i1 != L1.end(); i1++)
			{
				for (i2=L2.begin(); i2 != L2.end(); i2++)
				{
					*pProb += FREQ_MUTANT(2 * i1->Frequency * i2->Frequency,
						Geno._HamDist(Haplo.Num_SNP, *i1, *i2));
				}
			}
			pProb ++;
		}
	}

	// normalize
	TFLOAT sum = 0;
	TFLOAT *p = &_PostProb[0];
	for (size_t n = _PostProb.size(); n > 0; n--) sum += *p++;
	sum = 1.0 / sum;
	p = &_PostProb[0];
	for (size_t n = _PostProb.size(); n > 0; n--) *p++ *= sum;
}

THLAType CAlg_Prediction::_PredBestGuess(const CHaplotypeList &Haplo,
	const TGenotype &Geno)
{
	THLAType rv;
	rv.Allele1 = rv.Allele2 = NA_INTEGER;
	TFLOAT max=0, prob;

	vector<THaplotype>::const_iterator i1;
	vector<THaplotype>::const_iterator i2;

	for (int h1=0; h1 < _nHLA; h1++)
	{
		const vector<THaplotype> &L1 = Haplo.List[h1];

		// diag value
		prob = 0;
		for (i1=L1.begin(); i1 != L1.end(); i1++)
		{
			for (i2=i1; i2 != L1.end(); i2++)
			{
				prob += FREQ_MUTANT((i1 != i2) ?
					(2 * i1->Frequency * i2->Frequency) : (i1->Frequency * i2->Frequency),
					Geno._HamDist(Haplo.Num_SNP, *i1, *i2));
			}
		}
		if (max < prob)
		{
			max = prob;
			rv.Allele1 = rv.Allele2 = h1;
		}

		// off-diag value
		for (int h2=h1+1; h2 < _nHLA; h2++)
		{
			const vector<THaplotype> &L2 = Haplo.List[h2];
			prob = 0;
			for (i1=L1.begin(); i1 != L1.end(); i1++)
			{
				for (i2=L2.begin(); i2 != L2.end(); i2++)
				{
					prob += FREQ_MUTANT(2 * i1->Frequency * i2->Frequency,
						Geno._HamDist(Haplo.Num_SNP, *i1, *i2));
				}
			}
			if (max < prob)
			{
				max = prob;
				rv.Allele1 = h1; rv.Allele2 = h2;
			}
		}
	}

	return rv;
}

TFLOAT CAlg_Prediction::_PredPostProb(const CHaplotypeList &Haplo,
	const TGenotype &Geno, const THLAType &HLA)
{
	int H1=HLA.Allele1, H2=HLA.Allele2;
	if (H1 > H2) std::swap(H1, H2);
	int IxHLA = H2 + H1*(2*_nHLA-H1-1)/2;
	int idx = 0;

	TFLOAT sum=0, hlaProb=0, prob;
	vector<THaplotype>::const_iterator i1;
	vector<THaplotype>::const_iterator i2;

	for (int h1=0; h1 < _nHLA; h1++)
	{
		const vector<THaplotype> &L1 = Haplo.List[h1];

		// diag value
		prob = 0;
		for (i1=L1.begin(); i1 != L1.end(); i1++)
		{
			for (i2=i1; i2 != L1.end(); i2++)
			{
				prob += FREQ_MUTANT((i1 != i2) ?
					(2 * i1->Frequency * i2->Frequency) : (i1->Frequency * i2->Frequency),
					Geno._HamDist(Haplo.Num_SNP, *i1, *i2));
			}
		}
		if (IxHLA == idx) hlaProb = prob;
		idx ++; sum += prob;

		// off-diag value
		for (int h2=h1+1; h2 < _nHLA; h2++)
		{
			const vector<THaplotype> &L2 = Haplo.List[h2];
			prob = 0;
			for (i1=L1.begin(); i1 != L1.end(); i1++)
			{
				for (i2=L2.begin(); i2 != L2.end(); i2++)
				{
					prob += FREQ_MUTANT(2 * i1->Frequency * i2->Frequency,
						Geno._HamDist(Haplo.Num_SNP, *i1, *i2));
				}
			}
			if (IxHLA == idx) hlaProb = prob;
			idx ++; sum += prob;
		}
	}

	return hlaProb / sum;
}

THLAType CAlg_Prediction::BestGuess()
{
	THLAType rv;
	rv.Allele1 = rv.Allele2 = NA_INTEGER;

	TFLOAT *p = &_PostProb[0];
	TFLOAT max = 0;
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

	TFLOAT *p = &_SumPostProb[0];
	TFLOAT max = 0;
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
		_GenoList.List[i].BootstrapCount = _BootstrapCnt[i];
	_GenoList.Num_SNP = 0;

	_Predict.InitPrediction(nHLA());
}

void CVariableSelection::_InitHaplotype(CHaplotypeList &Haplo)
{
	vector<int> tmp(_HLAList->Num_HLA_Allele(), 0);
	int SumCnt = 0;
	for (int i=0; i < nSamp(); i++)
	{
		int cnt = _GenoList.List[i].BootstrapCount;
		tmp[_HLAList->List[i].Allele1] += cnt;
		tmp[_HLAList->List[i].Allele2] += cnt;
		SumCnt += cnt;
	}

	const TFLOAT scale = 0.5 / SumCnt;
	Haplo.Num_SNP = 0;
	Haplo.List.clear();
	Haplo.List.resize(_HLAList->Num_HLA_Allele());
	for (int i=0; i < (int)tmp.size(); i++)
	{
		if (tmp[i] > 0)
			Haplo.List[i].push_back(THaplotype(tmp[i] * scale));
	}
}

TFLOAT CVariableSelection::_OutOfBagAccuracy(CHaplotypeList &Haplo)
{
#if (HIBAG_TIMING == 1)
	_put_timing();
#endif

	HIBAG_CHECKING(Haplo.Num_SNP != _GenoList.Num_SNP,
		"CVariableSelection::_OutOfBagAccuracy, Haplo and GenoList should have the same number of SNP markers.");

	int TotalCnt, CorrectCnt;

#ifdef HIBAG_ALLOW_GPU_SUPPORT
	if (_gpuOutOfBagAcc_F64 || _gpuOutOfBagAcc_F32)
	{
		int *_HLA_HapIdx = NULL;
		void *_HapList = NULL;
		TGPU_Genotype *_GList = NULL;

		Haplo.InitGPUHostData(_HLA_HapIdx, _HapList);
		TotalCnt = InitGPUHostData_OutOfBag(_GList);

		if (_gpuRunInDoublePrecision)
		{
			CorrectCnt = (*_gpuOutOfBagAcc_F64)(nHLA(), _HLA_HapIdx,
				Haplo.TotalNumOfHaplo(), (TGPU_Haplotype_F64*)_HapList,
				TotalCnt, _GList, _GenoList.Num_SNP);
		} else {
			CorrectCnt = (*_gpuOutOfBagAcc_F32)(nHLA(), _HLA_HapIdx,
				Haplo.TotalNumOfHaplo(), (TGPU_Haplotype_F32*)_HapList,
				TotalCnt, _GList, _GenoList.Num_SNP);
		}

		Haplo.FreqGPUHostData(_HLA_HapIdx, _HapList);
		FreqGPUHostData(_GList);

		TotalCnt *= 2;
	} else {
#endif

		TotalCnt = CorrectCnt = 0;
		vector<TGenotype>::const_iterator it  = _GenoList.List.begin();
		vector<THLAType>::const_iterator pHLA = _HLAList->List.begin();

		for (; it != _GenoList.List.end(); it++, pHLA++)
		{
			if (it->BootstrapCount <= 0)
			{
				CorrectCnt += CHLATypeList::Compare(
					_Predict._PredBestGuess(Haplo, *it), *pHLA);
				TotalCnt += 2;
			}
		}

#ifdef HIBAG_ALLOW_GPU_SUPPORT
	}
#endif

#if (HIBAG_TIMING == 1)
	_inc_timing();
#endif

	return (TotalCnt>0) ? TFLOAT(CorrectCnt)/TotalCnt : 1;
}

TFLOAT CVariableSelection::_InBagLogLik(CHaplotypeList &Haplo)
{
#if (HIBAG_TIMING == 1)
	_put_timing();
#endif

	HIBAG_CHECKING(Haplo.Num_SNP != _GenoList.Num_SNP,
		"CVariableSelection::_InBagLogLik, Haplo and GenoList should have the same number of SNP markers.");

#ifdef HIBAG_ALLOW_GPU_SUPPORT

#endif

	vector<TGenotype>::const_iterator it   = _GenoList.List.begin();
	vector<THLAType>::const_iterator  pHLA = _HLAList->List.begin();
	TFLOAT LogLik = 0;

	for (; it != _GenoList.List.end(); it++, pHLA++)
	{
		if (it->BootstrapCount > 0)
		{
			LogLik += it->BootstrapCount *
				FLOAT_LOG(_Predict._PredPostProb(Haplo, *it, *pHLA));
		}
	}

#if (HIBAG_TIMING == 1)
	_inc_timing();
#endif
	return -2 * LogLik;
}

void CVariableSelection::Search(CBaseSampling &VarSampling,
	CHaplotypeList &OutHaplo, vector<int> &OutSNPIndex,
	TFLOAT &Out_Global_Max_OutOfBagAcc, int mtry, bool prune,
	bool verbose, bool verbose_detail)
{
	// rare probability
	const TFLOAT RARE_PROB = std::max(FRACTION_HAPLO/(2*nSamp()), MIN_RARE_FREQ);

	// initialize output
	_InitHaplotype(OutHaplo);
	OutSNPIndex.clear();

	// initialize internal variables
	TFLOAT Global_Max_OutOfBagAcc = 0;
	TFLOAT Global_Min_Loss = 1e+30;

	CHaplotypeList NextHaplo, NextReducedHaplo, MinHaplo;

	while ((VarSampling.TotalNum()>0) &&
		((int)OutSNPIndex.size() < HIBAG_MAXNUM_SNP_IN_CLASSIFIER))
	{
		// prepare for growing the individual classifier
		_EM.PrepareHaplotypes(OutHaplo, _GenoList, *_HLAList, NextHaplo);

		TFLOAT max_OutOfBagAcc = Global_Max_OutOfBagAcc;
		TFLOAT min_loss = Global_Min_Loss;
		int min_i = -1;

		// sample mtry from all candidate SNP markers
		VarSampling.RandomSelect(mtry);
		for (int i=0; i < VarSampling.NumOfSelection(); i++)
		{
			if (_EM.PrepareNewSNP(VarSampling[i], OutHaplo, *_SNPMat, _GenoList, NextHaplo))
			{
				// run EM algorithm
				_EM.ExpectationMaximization(NextHaplo);
				NextHaplo.EraseDoubleHaplos(RARE_PROB, NextReducedHaplo);

				// evaluate losses
				_GenoList.AddSNP(VarSampling[i], *_SNPMat);
				TFLOAT loss = 0;
				TFLOAT acc = _OutOfBagAccuracy(NextReducedHaplo);
				if (acc >= max_OutOfBagAcc)
					loss = _InBagLogLik(NextReducedHaplo);
				_GenoList.ReduceSNP();

				// compare
				if (acc > max_OutOfBagAcc)
				{
					min_i = i;
					min_loss = loss; max_OutOfBagAcc = acc;
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
				Rprintf("\t%-3d, added snp: %d, loss: %g, out-of-bag acc: %0.2f%%, # of haplo: %d\n",
					OutSNPIndex.size(), OutSNPIndex.back()+1,
					Global_Min_Loss, Global_Max_OutOfBagAcc*100, OutHaplo.TotalNumOfHaplo());
			}
		} else {
			// only keep "n_tmp - m" predictors
			VarSampling.RemoveSelection();
		}
	}
	
	Out_Global_Max_OutOfBagAcc = Global_Max_OutOfBagAcc;
}

#ifdef HIBAG_ALLOW_GPU_SUPPORT

int CVariableSelection::InitGPUHostData_OutOfBag(TGPU_Genotype *&_GList)
{
	int Cnt = 0;
	vector<TGenotype>::iterator it;
	for (it=_GenoList.List.begin(); it != _GenoList.List.end(); it++)
	{
		if (it->BootstrapCount <= 0)
			Cnt ++;
	}

	if (Cnt > 0)
	{
		_GList = new TGPU_Genotype[Cnt];
		TGPU_Genotype *p = _GList;
		vector<THLAType>::const_iterator pHLA = _HLAList->List.begin();

		for (it=_GenoList.List.begin(); it != _GenoList.List.end(); it++, pHLA++)
		{
			if (it->BootstrapCount <= 0)
			{
			#ifdef HIBAG_SSE_OPTIMIZE_HAMMING_DISTANCE
				UTYPE *s = it->PackedSNPs;
				UINT8 *d = p->PackedSNPs;
				for (size_t n=HIBAG_PACKED_UTYPE_MAXNUM_SNP; n > 0; n--)
					*d ++ = *s ++;
			#else
				memcpy(p->PackedSNPs, it->PackedSNPs, HIBAG_PACKED_UTYPE_MAXNUM_SNP);
			#endif
				p->BootstrapCount = 0;
				p->HLA = *pHLA;
				p ++;
			}
		}
	} else
		_GList = NULL;

	return Cnt;
}

int CVariableSelection::InitGPUHostData_InBag(TGPU_Genotype *&_GList)
{
	int Cnt = 0;
	vector<TGenotype>::iterator it;
	for (it=_GenoList.List.begin(); it != _GenoList.List.end(); it++)
	{
		if (it->BootstrapCount > 0)
			Cnt ++;
	}

	if (Cnt > 0)
	{
		_GList = new TGPU_Genotype[Cnt];
		TGPU_Genotype *p = _GList;
		vector<THLAType>::const_iterator pHLA = _HLAList->List.begin();

		for (it=_GenoList.List.begin(); it != _GenoList.List.end(); it++, pHLA++)
		{
			if (it->BootstrapCount > 0)
			{
			#ifdef HIBAG_SSE_OPTIMIZE_HAMMING_DISTANCE
				UTYPE *s = it->PackedSNPs;
				UINT8 *d = p->PackedSNPs;
				for (size_t n=HIBAG_PACKED_UTYPE_MAXNUM_SNP; n > 0; n--)
					*d ++ = *s ++;
			#else
				memcpy(p->PackedSNPs, it->PackedSNPs, HIBAG_PACKED_UTYPE_MAXNUM_SNP);
			#endif
				p->BootstrapCount = it->BootstrapCount;
				p->HLA = *pHLA;
				p ++;
			}
		}
	} else
		_GList = NULL;

	return Cnt;
}

inline void CVariableSelection::FreqGPUHostData(TGPU_Genotype *&_GList)
{
	// no worry if _GList == NULL
	delete []_GList;
	_GList = NULL;
}

#endif



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
	_Haplo.List.clear();
	_SNPIndex.clear();
	_OutOfBag_Accuracy = 0;
}

void CAttrBag_Classifier::Assign(int n_snp, const int snpidx[],
	const int samp_num[], int n_haplo, const TFLOAT *freq, const int *hla,
	char *const haplo[], TFLOAT *_acc)
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
	_Haplo.List.clear();
	_Haplo.List.resize(_Owner->nHLA());
	_Haplo.Num_SNP = n_snp;
	for (int i=0; i < n_haplo; i++)
	{
		_Haplo.List[hla[i]].push_back(THaplotype(haplo[i], freq[i]));
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
#if (HIBAG_TIMING > 0)
	_timing_ = 0;
	clock_t _start_time = clock();
#endif

	CSamplingWithoutReplace VarSampling;

	for (int k=0; k < nclassifier; k++)
	{
		VarSampling.Init(nSNP());

		CAttrBag_Classifier *I = NewClassifierBootstrap();
		I->Grow(VarSampling, mtry, prune, verbose, verbose_detail);
		if (verbose)
		{
			time_t tm; time(&tm);
			string s(ctime(&tm));
			s.erase(s.size()-1, 1);
			Rprintf(
				"%s, %3d individual classifier, out-of-bag acc: %0.2f%%, # of SNPs: %d, # of haplo: %d\n",
				s.c_str(), k+1, I->OutOfBag_Accuracy()*100, I->nSNP(), I->nHaplo());
		}
	}

#if (HIBAG_TIMING > 0)
	Rprintf("It took %0.2f seconds, in %0.2f%%.\n",
		((TFLOAT)_timing_)/CLOCKS_PER_SEC,
		((TFLOAT)_timing_) / (clock() - _start_time) * 100.0);
#endif
}

void CAttrBag_Model::PredictHLA(const int *genomat, int n_samp, int vote_method,
	int OutH1[], int OutH2[], TFLOAT OutMaxProb[],
	TFLOAT OutProbArray[], bool ShowInfo)
{
	if ((vote_method < 1) || (vote_method > 2))
		throw ErrHLA("Invalid 'vote_method'.");

	const int nPairHLA = nHLA()*(nHLA()+1)/2;

	_Predict.InitPrediction(nHLA());
	Progress.Info = "Predicting:";
	Progress.Init(n_samp, ShowInfo);

	vector<int> Weight(nSNP());
	_GetSNPWeights(&Weight[0]);

	for (int i=0; i < n_samp; i++, genomat+=nSNP())
	{
		_PredictHLA(genomat, &Weight[0], vote_method);

		THLAType HLA = _Predict.BestGuessEnsemble();
		OutH1[i] = HLA.Allele1; OutH2[i] = HLA.Allele2;

		if ((HLA.Allele1 != NA_INTEGER) && (HLA.Allele2 != NA_INTEGER))
			OutMaxProb[i] = _Predict.IndexSumPostProb(HLA.Allele1, HLA.Allele2);
		else
			OutMaxProb[i] = 0;

		if (OutProbArray)
		{
			for (int j=0; j < nPairHLA; j++)
				*OutProbArray++ = _Predict.SumPostProb()[j];
		}

		Progress.Forward(1, ShowInfo);
	}
}

void CAttrBag_Model::PredictHLA_Prob(const int *genomat, int n_samp,
	int vote_method, TFLOAT OutProb[], bool ShowInfo)
{
	if ((vote_method < 1) || (vote_method > 2))
		throw ErrHLA("Invalid 'vote_method'.");

	const int n = nHLA()*(nHLA()+1)/2;
	_Predict.InitPrediction(nHLA());
	Progress.Info = "Predicting:";
	Progress.Init(n_samp, ShowInfo);

	vector<int> Weight(nSNP());
	_GetSNPWeights(&Weight[0]);

	for (int i=0; i < n_samp; i++, genomat+=nSNP())
	{
		_PredictHLA(genomat, &Weight[0], vote_method);
		for (int j=0; j < n; j++)
			*OutProb++ = _Predict.SumPostProb()[j];
		Progress.Forward(1, ShowInfo);
	}
}

void CAttrBag_Model::_PredictHLA(const int *geno, const int weights[],
	int vote_method)
{
	TGenotype Geno;
	_Predict.InitSumPostProbBuffer();

	// missing proportion
	vector<CAttrBag_Classifier>::const_iterator it;
	for (it = _ClassifierList.begin(); it != _ClassifierList.end(); it++)
	{
		const int n = it->nSNP();
		int nWeight=0, SumWeight=0;
		for (int i=0; i < n; i++)
		{
			int k = it->_SNPIndex[i];
			SumWeight += weights[k];
			if ((0 <= geno[k]) && (geno[k] <= 2))
				nWeight += weights[k];
		}

		/// set weight with respect to missing SNPs
		if (nWeight > 0)
		{
			Geno.IntToSNP(n, geno, &(it->_SNPIndex[0]));
			_Predict.PredictPostProb(it->_Haplo, Geno);

			if (vote_method == 1)
			{
				// predicting based on the averaged posterior probabilities
				_Predict.AddProbToSum(TFLOAT(nWeight) / SumWeight);
			} else if (vote_method == 2)
			{
				// predicting by class majority voting
				THLAType pd = _Predict.BestGuess();
				if ((pd.Allele1 != NA_INTEGER) && (pd.Allele2 != NA_INTEGER))
				{
					_Predict.InitPostProbBuffer();  // fill by ZERO
					_Predict.IndexPostProb(pd.Allele1, pd.Allele2) = 1.0;

					// _Predict.AddProbToSum(TFLOAT(nWeight) / SumWeight);
					_Predict.AddProbToSum(1.0);
				}
			}
		}
	}

	_Predict.NormalizeSumPostProb();
}

void CAttrBag_Model::_GetSNPWeights(int OutWeight[])
{
	// ZERO
	memset(OutWeight, 0, sizeof(int)*nSNP());
	// for each classifier
	vector<CAttrBag_Classifier>::const_iterator it;
	for (it = _ClassifierList.begin(); it != _ClassifierList.end(); it++)
	{
		const int n = it->nSNP();
		for (int i=0; i < n; i++)
			OutWeight[ it->_SNPIndex[i] ] ++;
	}
}



// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

template <typename TO, typename FROM> TO nasty_cast(FROM f)
{
	union {
		FROM f; TO t;
	} u;
	u.f = f;
	return u.t;
}

#if defined(HIBAG_SYS_UNIX)

	static void* GDS_Handle = NULL;
	#define LOAD(var, type, name)	\
		var = nasty_cast<type>(dlsym(GDS_Handle, name)); \
		if ((err = dlerror()) != NULL) \
			{ Done_GPU_Support(); throw ErrHLA(err); }

#elif defined(HIBAG_SYS_WIN)

	static HMODULE GDS_Handle = NULL;
	#define LOAD(var, type, name)	\
		var = nasty_cast<type>(GetProcAddress(GDS_Handle, name)); \
		if (var == NULL) \
			{ Done_GPU_Support(); throw ErrHLA("No %s function.", name); }

#else
#  error "not supported"
#endif


//  throw errors if it fails
void HLA_LIB::Init_GPU_Support(const char *lib_fn)
{
	if (!GDS_Handle)
	{
	#if defined(HIBAG_SYS_UNIX)
		// open dll
		GDS_Handle = dlopen(lib_fn, RTLD_LAZY);
		if (!GDS_Handle) throw ErrHLA(dlerror());
		dlerror();  // Clear any existing error
		const char *err;
	#elif defined(HIBAG_SYS_WIN)
		GDS_Handle = LoadLibrary(lib_fn);
		if (GDS_Handle == NULL)
		{
			char buf[1024];
			memset((void*)buf, 0, sizeof(buf));
			FormatMessage(FORMAT_MESSAGE_FROM_SYSTEM |
				FORMAT_MESSAGE_IGNORE_INSERTS | FORMAT_MESSAGE_ARGUMENT_ARRAY,
				NULL, GetLastError(), 0, buf, sizeof(buf), NULL);
			throw ErrHLA(buf);
		}
	#else
	#  error "not supported"
	#endif

		// load functions
		LOAD(_gpuInitialize, TGPUInitFunc, "hlaGPU_Initialize");
		LOAD(_gpuFinalize, TGPUFunc, "hlaGPU_Finalize");
		LOAD(_gpuDeviceQuery, TGPUFunc, "hlaGPU_DeviceQuery");
		LOAD(_gpuOutOfBagAcc_F32, TGPU_OutOfBagAcc_F32, "hlaGPU_OutOfBagAcc_F32");
		LOAD(_gpuOutOfBagAcc_F64, TGPU_OutOfBagAcc_F64, "hlaGPU_OutOfBagAcc_F64");

		// initialize GPU
		(*_gpuInitialize)(&_Packed_HammingDistance[0][0]);
		(*_gpuDeviceQuery)();
	}
}


// Finalize the GPU computing library
void HLA_LIB::Done_GPU_Support()
{
	if (_gpuFinalize) (*_gpuFinalize)();

	_gpuInitialize = NULL;
	_gpuFinalize = NULL;
	_gpuDeviceQuery = NULL;

#if defined(HIBAG_SYS_UNIX)

	if (GDS_Handle)
	{
		dlclose(GDS_Handle);
		GDS_Handle = NULL;
	}

#elif defined(HIBAG_SYS_WIN)

	if (GDS_Handle != NULL)
	{
		FreeLibrary(GDS_Handle);
		GDS_Handle = NULL;
	}

#else
#  error "not supported"
#endif
}

