// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
// Name        : LibHLA
// Author      : Xiuwen Zheng
// Version     : 1.0.0.0
// Copyright   : Xiuwen Zheng (GPL v3.0)
// Created     : 06/21/2012
// Description : An Attribute Bagging Method for Imputing HLA Types with SNP Data
// ===========================================================

#ifndef _LibHLA_H_
#define _LibHLA_H_

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <vector>
#include <list>
#include <string>
#include <algorithm>


#define HLA_CHECK_CODE


namespace HLA_LIB
{
	using namespace std;


	// ******************************************************************************* //
	// ********                          description                          ********
	// SNP allele: 0, 1
	// SNP genotype: 0, 1, 2, 3 (missing, or other values)
	// HLA allele: start from 0
	//
	// Packed SNP alleles: s8 s7 s6 s5 s4 s3 s2 s1
	//     the 1st allele: (s1), the 2nd allele: (s3)
	//     the 3rd allele: (s5), the 4th allele: (s7)
	// Packed SNP genotype: s8 s7 s6 s5 s4 s3 s2 s1
	//     the 1st genotype: (s2 s1), the 2nd genotype: (s4 s3)
	//     the 3rd genotype: (s6 s5), the 4th genotype: (s8 s7)
	// ********                                                               ********
	// ******************************************************************************* //



	// ******************************************************************************* //
	// ******************************************************************************* //

	/// macro for checking error
	#define CHECKING(x, msg)	{ if (x) throw ErrHLA(msg); }



	// ******************************************************************************* //
	// ********                           container                           ********

	/// Define an unsigned integer of 8 bits
	typedef unsigned char	UINT8;

	/// The number of bytes for packed SNP genotypes
	#define PACKED_NUM_IN_BYTE(x)	( ((x)/4) + (((x) % 4)!=0 ? 1 : 0) )

	/// The max number of SNP markers in an individual classifier
	static const int MAXNUM_SNP_IN_CLASSIFIER = 128;
	/// The max number of bytes for packed SNP genotypes
	static const int PACKED_BYTE_MAXNUM_SNP = PACKED_NUM_IN_BYTE(MAXNUM_SNP_IN_CLASSIFIER);


	/// Packed SNP haplotype structure: 4 SNPs in a byte
	struct THaplotype
	{
	public:
		UINT8 PackedSNPs[PACKED_BYTE_MAXNUM_SNP];  //< packed SNP genotypes
		double Frequency;  //< haplotype frequency
		double OldFreq;  //< old haplotype frequency

		THaplotype() {}
		THaplotype(const double _freq);
		THaplotype(const char *str, const double _freq);

		UINT8 GetSNP(int idx) const;
		void SetSNP(int idx, UINT8 val);
		string SNPToString(int Length) const;
		void StringToSNP(const string &str);
	};

	/// A list of haplotypes
	class CHaplotypeList
	{
	public:	
		CHaplotypeList();

		void DoubleHaplos(CHaplotypeList &OutHaplos) const;
		void DoubleHaplosInitFreq(CHaplotypeList &OutHaplos, const double AFreq) const;
		void MergeDoubleHaplos(const double RareProb, CHaplotypeList &OutHaplos) const;
		void EraseDoubleHaplos(const double RareProb, CHaplotypeList &OutHaplos) const;
		void SaveClearFrequency();
		void ScaleFrequency(const double scale);
		int TotalNumOfHaplo() const;
		void Print();

		vector< vector<THaplotype> > List; //< haplotype list
		int Num_SNP; //< the number of SNP markers
	};


	/// Packed SNP genotype structure: 4 SNPs in a byte
	class TGenotype
	{
	public:
		friend class CGenotypeList;
		friend class CAlg_EM;
		friend class CAlg_Prediction;
	
		UINT8 PackedSNPs[PACKED_BYTE_MAXNUM_SNP];  //< packed SNP genotypes
		int BootstrapCount;  //< the count in the bootstrapped data

		UINT8 GetSNP(const int idx) const;
		void SetSNP(const int idx, int val);
		string SNPToString(const int Length) const;
		void StringToSNP(const string &str);
		void SNPToInt(const int Length, int OutArray[]) const;
		void IntToSNP(int Length, const int InArray[]);
		void IntToSNP(int Length, const int InBase[], const int Index[]);

		int Diff(int Length, const THaplotype &H1, const THaplotype &H2) const;

	protected:
		void _SetSNP(const int idx, UINT8 val);
		int _Diff(int Length, const THaplotype &H1, const THaplotype &H2) const;
	};

	/// SNP genotype container
	class CSNPGenoMatrix
	{
	public:
		CSNPGenoMatrix();

		const int Get(const int IdxSamp, const int IdxSNP) const;
		int *Get(const int IdxSamp);

		int Num_Total_SNP, Num_Total_Samp;
		int *pGeno;
	};

	/// A list of genotypes
	class CGenotypeList
	{
	public:
		CGenotypeList();

		void AddSNP(int IdxSNP, const CSNPGenoMatrix &SNPMat);
		void ReduceSNP();
		void Print();
		/// return the total number of samples
		inline const int nSamp() const { return List.size(); }

		vector<TGenotype> List;  //< genotype list
		int Num_SNP;  //< the number of SNP markers
	protected:
	};


	/// A pair of HLA alleles
	struct THLAType
	{
	public:
		int Allele1, Allele2;
	};

	/// A list of HLA types
	class CHLATypeList
	{
	public:
		CHLATypeList();

		/// return the total number of samples
		inline const int nSamp() const { return List.size(); }
		/// return the total number of HLA alleles
		inline const int Num_HLA_Allele() const { return Str_HLA_Allele.size(); }
		
		/// return how many alleles are the same
		static int Compare(const THLAType &H1, const THLAType &H2);

		vector<THLAType> List;  //< a list of HLA types
		vector<string> Str_HLA_Allele;  //< HLA alleles
	protected:
	};



	// ******************************************************************************* //
	// ********                           algorithm                           ********

	// the parameter of EM algorithm for estimating haplotype frequencies
	/// The max number of iterations
	extern int EM_MaxNum_Iterations;  // 500
	/// The reltol convergence tolerance, sqrt(machine.epsilon) by default, used in EM algorithm
	extern double EM_FuncRelTol;  // sqrt(DBL_EPSILON)


	/// variable sampling
	class CBaseSampling
	{
	public:
		virtual int TotalNum() const = 0;
		virtual void RandomSelect(int m_try) = 0;
		virtual int NumOfSelection() const = 0;
		virtual void Remove(int idx) = 0;
		virtual void RemoveSelection() = 0;
		virtual void RemoveFlag() = 0;
		virtual int &operator[] (int idx) = 0;
	};

	/// variable sampling with flat prior probability
	class CSamplingWithoutReplace: public CBaseSampling
	{
	public:
		CSamplingWithoutReplace();
		CBaseSampling *Init(int m_total);

		virtual int TotalNum() const;
		virtual void RandomSelect(int m_try);
		virtual int NumOfSelection() const;
		virtual void Remove(int idx);
		virtual void RemoveSelection();
		virtual void RemoveFlag();
		virtual int &operator[] (int idx);
	protected:
		vector<int> _IdxArray;
		int _m_try;
	};

	/// variable sampling with prior probability
	class CSamplingWithoutReplaceWithProb: public CBaseSampling
	{
	public:
		CSamplingWithoutReplaceWithProb();
		CBaseSampling *Init(int m_total, double *var_prob);

		virtual int TotalNum() const;
		virtual void RandomSelect(int m_try);
		virtual int NumOfSelection() const;
		virtual void Remove(int idx);
		virtual void RemoveSelection();
		virtual void RemoveFlag();
		virtual int &operator[] (int idx);
	protected:
		struct type { int index; double prob; };
		vector<type> _IdxArray;
		int _m_try;
		double _prob_sum;
	};


	/// Expectation Maximization algorithm for estimating haplotype frequencies
	class CAlg_EM
	{
	public:
		CAlg_EM();

		// call PrepareHaplotypes first, and then call PrepareNewSNP

		void PrepareHaplotypes(const CHaplotypeList &CurHaplo,
			const CGenotypeList &GenoList, const CHLATypeList &HLAList,
			CHaplotypeList &NextHaplo);
		/// , return true if the new SNP is not monomorphic
		bool PrepareNewSNP(const int NewSNP, const CHaplotypeList &CurHaplo,
			const CSNPGenoMatrix &SNPMat, CGenotypeList &GenoList, CHaplotypeList &NextHaplo);

		void ExpectationMaximization(CHaplotypeList &NextHaplo);

	protected:
		/// A pair of haplotypes
		struct THaploPair
		{
			bool Flag;  //< if true, the haplotype pair exists in the sample
			THaplotype *H1;  //< the first haplotype
			THaplotype *H2;  //< the second haplotype
			double Freq;  //< genotype frequency

			THaploPair() { Flag = true; H1 = H2 = NULL; }
			THaploPair(THaplotype *i1, THaplotype *i2) { Flag = true; H1 = i1; H2 = i2; }
		};
		/// A list of haplotype pairs
		struct THaploPairList
		{
			int BootstrapCount;  //< the count in the bootstrapped data
			int SampIndex;  //< the sample index in the source data
			vector<THaploPair> PairList;  //< a list of haplotype pairs
			
			void Print(const int Length);
		};

		vector<THaploPairList> _SampHaploPair;
	};

	/// the algorithm prediction
	class CAlg_Prediction
	{
	public:
		CAlg_Prediction();

		void InitPrediction(int n_hla);

		void InitSumPostProb();
		void AddProbToSum(const double weight);
		void NormalizeSumPostProb();
		
		double &IndexPostProb(int H1, int H2);
		double &IndexSumPostProb(int H1, int H2);

		void Predict(const CHaplotypeList &Haplo, const TGenotype &Geno);
		THLAType MaxProb();
		THLAType MaxSumProb();

		inline const int nHLA() const { return _nHLA; }
		inline const vector<double> &PostProb() const { return _PostProb; }
		inline const vector<double> &SumPostProb() const { return _SumPostProb; }

	protected:
		int _nHLA;  //< the number of different HLA alleles
		double _Sum_Weight;  //< +1 after calling AddProbToSum()
		vector<double> _PostProb;
		vector<double> _SumPostProb;
	};


	/// the algorithm prediction
	class CVariableSelection
	{
	public:
		CVariableSelection();

		void InitSelection(CSNPGenoMatrix &snpMat, CHLATypeList &hlaList, const int _BootstrapCnt[]);
		void Search(CBaseSampling &VarSampling, CHaplotypeList &OutHaplo, vector<int> &OutSNPIndex,
			double &Out_Global_Max_OutOfBagAcc,
			int mtry, bool prune, bool verbose, bool verbose_detail, bool debug);

		inline const int nSamp() const { return _SNPMat->Num_Total_Samp; }
		inline const int nSNP() const { return _SNPMat->Num_Total_SNP; }
		inline const int nHLA() const { return _HLAList->Num_HLA_Allele(); }

	protected:
		CSNPGenoMatrix *_SNPMat;
		CHLATypeList *_HLAList;
		
		CGenotypeList _GenoList;
		CAlg_EM _EM;
		CAlg_Prediction _Predict;
	
		void _InitHaplotype(CHaplotypeList &Haplo);
		double _OutOfBagAccuracy(CHaplotypeList &Haplo);
		double _InBagLogLik(CHaplotypeList &Haplo);
	};




	// ******************************************************************************* //
	// ********                             model                             ********

	class CAttrBag_Model;

	/// the individual classifier
	class CAttrBag_Classifier
	{
	public:
		friend class CAttrBag_Model;

		CAttrBag_Classifier(CAttrBag_Model &_owner);

		void InitBootstrapCount(int SampCnt[]);
		void Assign(int n_snp, const int snpidx[], const int samp_num[], int n_haplo,
			const double *freq, const int *hla, char *const haplo[], double *_acc=NULL);
		void Grow(CBaseSampling &VarSampling, int mtry, bool prune,
			bool verbose, bool verbose_detail, bool debug);

		inline CAttrBag_Model &Owner() { return *_Owner; }
		inline const int nSNP() const { return _SNPIndex.size(); }
		inline const int nHaplo() const { return _Haplo.TotalNumOfHaplo(); }
		inline const double OutOfBag_Accuracy() const { return _OutOfBag_Accuracy; }
		inline const vector<int> &SNPIndex() const { return _SNPIndex; }
		inline const vector<int> &BootstrapCount() const { return _BootstrapCount; }
		inline const CHaplotypeList &Haplotype() const { return _Haplo; }

	protected:
		CAttrBag_Model *_Owner;
		CHaplotypeList _Haplo;
		vector<int> _BootstrapCount;
		vector<int> _SNPIndex;
		double _OutOfBag_Accuracy;
	};

	/// the attribute bagging model
	class CAttrBag_Model
	{
	public:
		friend class CAttrBag_Classifier;

		CAttrBag_Model();

		void InitTraining(int n_snp, int n_samp, int n_hla);
		void InitTraining(int n_snp, int n_samp, int *snp_geno, int n_hla, int *H1, int *H2);

		CAttrBag_Classifier *NewClassifier();
		CAttrBag_Classifier *NewClassifierAllSamp();

		void RunNewClassifiers(int nclassifier, int mtry, double *var_prob,
			bool prune, bool verbose, bool verbose_detail=false, bool debug=false);

		void PredictHLA(const int *genomat, int n_samp,
			int OutH1[], int OutH2[], double OutProb[], bool ShowInfo);
		void PredictHLA_Prob(const int *genomat, int n_samp,
			double OutProb[], bool ShowInfo);

		inline const int nSamp() const { return _SNPMat.Num_Total_Samp; }
		inline const int nSNP() const { return _SNPMat.Num_Total_SNP; }
		inline const int nHLA() const { return _HLAList.Num_HLA_Allele(); }

		inline const CHLATypeList &HLAList() const { return _HLAList; }
		inline const vector<CAttrBag_Classifier> &ClassifierList() const { return _ClassifierList; }

	protected:
		CSNPGenoMatrix _SNPMat;
		CHLATypeList _HLAList;
		vector<CAttrBag_Classifier> _ClassifierList;
		CVariableSelection _VarSelect;
		CAlg_Prediction _Predict;

		void _PredictHLA(const int *geno, const int weights[]);
		void _GetSNPWeights(int OutWeight[]);
	};



	// ******************************************************************************* //
	// ******************************************************************************* //

	/// The basic class for progress object
	class CdProgression
	{
	public:
		static const int TotalPercent = 10;
		static const int StepPercent = 10;
		/// The associated information
		std::string Info;

		/// Constructor
		CdProgression();

		void Init(long TotalCnt, bool ShowInit);
		bool Forward(long step, bool Show);
		virtual void ShowProgress();

        /// Return the current percentile
		inline int Percent() const { return fPercent; }
		/// Return the total number
		inline long Total() const { return fTotal; }
		/// Return the current position
		inline long Current() const { return fCurrent; }
	protected:
		long fTotal, fCurrent;
		int fPercent;
		clock_t OldTime;
	};

	/// The progression information
	extern CdProgression Progress;


	/// Exceptions for HLA imputation
	class ErrHLA: public std::exception
	{
	public:
		ErrHLA() {};
		ErrHLA(const std::string &msg) { fMessage = msg; }
		ErrHLA(const char *fmt, ...)
			{
				va_list args;
				va_start(args, fmt);
				char buf[1024];
				vsnprintf(buf, sizeof(buf), fmt, args);
				fMessage = buf;
				va_end(args);
			}
		virtual const char *what() const throw() { return fMessage.c_str(); }
		virtual ~ErrHLA() throw() {};
	protected:
		std::string fMessage;
	};
}


#ifdef  __cplusplus
extern "C" {
#endif
	void Rprintf(const char *, ...);
#ifdef  __cplusplus
}
#endif


#endif /* _LibHLA_H_ */
