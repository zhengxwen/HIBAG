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
// This file is part of HIBAG package.
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
// Name        : LibHLA
// Author      : Xiuwen Zheng
// Version     : 1.1.0
// Copyright   : Xiuwen Zheng (GPL v3.0)
// Created     : 11/14/2011
// Description : HLA Genotype Imputation with Attribute Bagging
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
	// Packed SNP storage strategy is used for faster matching
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
	static const int MAXNUM_SNP_IN_CLASSIFIER = 256;
	/// The max number of bytes for packed SNP genotypes
	static const int PACKED_BYTE_MAXNUM_SNP = PACKED_NUM_IN_BYTE(MAXNUM_SNP_IN_CLASSIFIER);


	/// Packed SNP haplotype structure: 4 SNPs in a byte
	struct THaplotype
	{
	public:
		/// packed SNP alleles
		UINT8 PackedSNPs[PACKED_BYTE_MAXNUM_SNP];
		/// haplotype frequency
		double Frequency;
		/// old haplotype frequency
		double OldFreq;

		THaplotype() {}
		THaplotype(const double _freq);
		THaplotype(const char *str, const double _freq);

		/// get SNP allele, idx starts from ZERO
		UINT8 GetSNP(int idx) const;
		/// set SNP allele, idx starts from ZERO
		void SetSNP(int idx, UINT8 val);
		/// get a string of "0" and "1" from packed SNP alleles
		string SNPToString(const int Length) const;
		/// set packed SNP alleles from a string of "0" and "1"
		void StringToSNP(const string &str);
	};


	/// A list of haplotypes
	class CHaplotypeList
	{
	public:	
		CHaplotypeList();

		/// initialize haplotypes for EM algorithm
		void DoubleHaplos(CHaplotypeList &OutHaplos) const;
		/// initialize haplotypes for EM algorithm with the allele frequency of the new SNP
		void DoubleHaplosInitFreq(CHaplotypeList &OutHaplos, const double AFreq) const;
		/// merge rare haplotypes
		void MergeDoubleHaplos(const double RareProb, CHaplotypeList &OutHaplos) const;
		/// remove rare haplotypes
		void EraseDoubleHaplos(const double RareProb, CHaplotypeList &OutHaplos) const;
		/// save frequency to old.frequency, and set current freq. to be 0
		void SaveClearFrequency();
		/// scale the haplotype frequencies by a factor
		void ScaleFrequency(const double scale);
		/// the total number of haplotypes
		int TotalNumOfHaplo() const;
		/// print all haplotypes
		void Print();

		/// haplotype list
		vector< vector<THaplotype> > List;
		/// the number of SNP markers
		int Num_SNP;
	};


	/// Packed SNP genotype structure: 4 SNPs in a byte
	class TGenotype
	{
	public:
		friend class CGenotypeList;
		friend class CAlg_EM;
		friend class CAlg_Prediction;

		/// packed SNP genotypes
		UINT8 PackedSNPs[PACKED_BYTE_MAXNUM_SNP];
		/// the count in the bootstrapped data
		int BootstrapCount;

		/// get SNP genotype (0, 1, 2) at the specified locus, idx starts from ZERO
		UINT8 GetSNP(const int idx) const;
		/// set SNP genotype (0, 1, 2) at the specified locus, idx starts from ZERO
		void SetSNP(const int idx, int val);
		/// get a string of SNP genotypes, consisting of "0" or "1"
		string SNPToString(const int Length) const;
		/// set SNP genotypes by a string of "0" and "1"
		void StringToSNP(const string &str);
		/// export SNPs to a vector of integers
		void SNPToInt(const int Length, int OutArray[]) const;
		/// import SNPs from a vector of integers
		void IntToSNP(int Length, const int InArray[]);
		/// import SNPs from an integer vector 'InBase' with 'Index'
		void IntToSNP(int Length, const int InBase[], const int Index[]);
		/// compute the Hamming distance between SNPs and H1+H2
		int Diff(int Length, const THaplotype &H1, const THaplotype &H2) const;

	protected:
		/// set SNP genotype (0, 1, 2) without checking
		void _SetSNP(const int idx, UINT8 val);
		/// compute the Hamming distance between SNPs and H1+H2 without checking
		int _Diff(int Length, const THaplotype &H1, const THaplotype &H2) const;
	};


	/// SNP genotype container
	class CSNPGenoMatrix
	{
	public:
		CSNPGenoMatrix();

		/// get SNP genotype of 'IdxSamp' individual and 'IdxSNP' SNP
		const int Get(const int IdxSamp, const int IdxSNP) const;
		/// get the pointer to SNP genotypes of 'IdxSamp' individual
		int *Get(const int IdxSamp);

		/// the total number of SNPs
		int Num_Total_SNP;
		/// the total number of samples
		int Num_Total_Samp;
		/// the pointer to SNP genotypes
		int *pGeno;
	};


	/// A list of genotypes
	class CGenotypeList
	{
	public:
		CGenotypeList();

		/// add a new SNP
		void AddSNP(int IdxSNP, const CSNPGenoMatrix &SNPMat);
		/// remove the last SNP
		void ReduceSNP();
		/// print all SNP gentoypes
		void Print();

		/// return the total number of samples
		inline const int nSamp() const { return List.size(); }

		/// genotype list
		vector<TGenotype> List;
		/// the number of SNPs
		int Num_SNP;
	protected:
	};


	/// A pair of HLA alleles
	struct THLAType
	{
	public:
		int Allele1;  //< the first HLA allele
		int Allele2;  //< the second HLA allele
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

		/// a list of HLA types
		vector<THLAType> List;
		/// HLA alleles
		vector<string> Str_HLA_Allele;
	protected:
	};




	// ******************************************************************************* //
	// ********                           algorithm                           ********

	// the parameter of EM algorithm for estimating haplotype frequencies

	/// The max number of iterations
	extern int EM_MaxNum_Iterations;  // = 500

	/// The reltol convergence tolerance, sqrt(machine.epsilon) by default, used in EM algorithm
	extern double EM_FuncRelTol;  // = sqrt(DBL_EPSILON)


	/// variable sampling
	class CBaseSampling
	{
	public:
		/// the total number of candidate SNPs
		virtual int TotalNum() const = 0;
		/// randomly select 'm_try' SNPs for further searching
		virtual void RandomSelect(int m_try) = 0;
		/// the number of selected SNPs
		virtual int NumOfSelection() const = 0;
		/// remove the 'idx' SNP
		virtual void Remove(int idx) = 0;
		/// remove the selected SNPs
		virtual void RemoveSelection() = 0;
		/// remove the SNPs with flag
		virtual void RemoveFlag() = 0;
		/// get SNP index
		virtual int &operator[] (int idx) = 0;
	};


	/// variable sampling without replacement
	class CSamplingWithoutReplace: public CBaseSampling
	{
	public:
		CSamplingWithoutReplace();
		CBaseSampling *Init(int m_total);

		/// the total number of candidate SNPs
		virtual int TotalNum() const;
		/// randomly select 'm_try' SNPs for further searching
		virtual void RandomSelect(int m_try);
		/// the number of selected SNPs
		virtual int NumOfSelection() const;
		/// remove the 'idx' SNP
		virtual void Remove(int idx);
		/// remove the selected SNPs
		virtual void RemoveSelection();
		/// remove the SNPs with flag
		virtual void RemoveFlag();
		/// get SNP index
		virtual int &operator[] (int idx);

	protected:
		/// saving the SNP indices
		vector<int> _IdxArray;
		/// the number of selected SNPs
		int _m_try;
	};


	/// Expectation Maximization algorithm for estimating haplotype frequencies
	class CAlg_EM
	{
	public:
		CAlg_EM();

		// call PrepareHaplotypes first, and then call PrepareNewSNP

		/// 
		void PrepareHaplotypes(const CHaplotypeList &CurHaplo,
			const CGenotypeList &GenoList, const CHLATypeList &HLAList,
			CHaplotypeList &NextHaplo);

		/// , return true if the new SNP is not monomorphic
		bool PrepareNewSNP(const int NewSNP, const CHaplotypeList &CurHaplo,
			const CSNPGenoMatrix &SNPMat, CGenotypeList &GenoList, CHaplotypeList &NextHaplo);

		/// call EM algorithm to estimate haplotype frequencies
		void ExpectationMaximization(CHaplotypeList &NextHaplo);

	protected:
		/// A pair of haplotypes
		struct THaploPair
		{
			bool Flag;       //< if true, the haplotype pair exists in the sample
			THaplotype *H1;  //< the first haplotype
			THaplotype *H2;  //< the second haplotype
			double Freq;     //< genotype frequency

			THaploPair() { Flag = true; H1 = H2 = NULL; }
			THaploPair(THaplotype *i1, THaplotype *i2) { Flag = true; H1 = i1; H2 = i2; }
		};

		/// A list of haplotype pairs
		struct THaploPairList
		{
			int BootstrapCount;           //< the count in the bootstrapped data
			int SampIndex;                //< the sample index in the source data
			vector<THaploPair> PairList;  //< a list of haplotype pairs
			
			/// print information
			void Print(const int Length);
		};

		/// pairs of haplotypes for individuals
		vector<THaploPairList> _SampHaploPair;
	};


	/// the algorithm prediction
	class CAlg_Prediction
	{
	public:
		CAlg_Prediction();

		/// initialize
		/** \param n_hla    the number of unique HLA alleles **/
		void InitPrediction(int n_hla);

		/// initialize the posterior probabilities by setting ZERO
		void InitPostProbBuffer();
		/// initialize the sums of posterior probabilities by setting ZERO
		void InitSumPostProbBuffer();
		/// add the posterior probabilities of a classifier with a weight to the variable for summing up
		void AddProbToSum(const double weight);
		/// average over all classifiers
		void NormalizeSumPostProb();

		/// 
		double &IndexPostProb(int H1, int H2);
		/// 
		double &IndexSumPostProb(int H1, int H2);

		/// predict based on SNP profiles and haplotype list
		void Predict(const CHaplotypeList &Haplo, const TGenotype &Geno);
		/// the best-guess HLA type for the current classifier
		THLAType MaxProb();
		/// the best-guess HLA type from all classifiers
		THLAType MaxSumProb();

		/// get the number of unique HLA alleles
		inline const int nHLA() const
			{ return _nHLA; }
		/// the posterior probabilities for the current classifier
		inline const vector<double> &PostProb() const
			{ return _PostProb; }
		/// the average posterior probabilities for all classifiers
		inline const vector<double> &SumPostProb() const
			{ return _SumPostProb; }

	protected:
		/// the number of different HLA alleles
		int _nHLA;
		/// plus weight after calling AddProbToSum()
		double _Sum_Weight;
		/// a vector of posterior probabilities
		vector<double> _PostProb;
		/// a vector of posterior probabilities for summing up
		vector<double> _SumPostProb;
	};


	/// variable selection algorithm
	class CVariableSelection
	{
	public:
		CVariableSelection();

		/// initialize
		void InitSelection(CSNPGenoMatrix &snpMat, CHLATypeList &hlaList,
			const int _BootstrapCnt[]);
		/// searching algorithm
		void Search(CBaseSampling &VarSampling, CHaplotypeList &OutHaplo,
			vector<int> &OutSNPIndex, double &Out_Global_Max_OutOfBagAcc,
			int mtry, bool prune, bool verbose, bool verbose_detail, bool debug);

		/// the number of samples
		inline const int nSamp() const { return _SNPMat->Num_Total_Samp; }
		/// the number of SNPs
		inline const int nSNP() const { return _SNPMat->Num_Total_SNP; }
		/// the number of unique HLA alleles
		inline const int nHLA() const { return _HLAList->Num_HLA_Allele(); }

	protected:
		/// store the genotype matrix
		CSNPGenoMatrix *_SNPMat;
		/// a list of HLA types
		CHLATypeList *_HLAList;
		
		/// a list of genotypes
		CGenotypeList _GenoList;
		/// EM algorithm
		CAlg_EM _EM;
		/// the prediction algorithm
		CAlg_Prediction _Predict;

		/// initialize the haplotype list
		void _InitHaplotype(CHaplotypeList &Haplo);
		/// compute the out-of-bag accuracy using the haplotypes 'Haplo'
		double _OutOfBagAccuracy(CHaplotypeList &Haplo);
		/// compute the in-bag log likelihood using the haplotypes 'Haplo'
		double _InBagLogLik(CHaplotypeList &Haplo);
	};




	// ******************************************************************************* //
	// ********                        HIBAG -- model                         ********

	class CAttrBag_Model;

	/// the individual classifier of HIBAG
	class CAttrBag_Classifier
	{
	public:
		friend class CAttrBag_Model;

		CAttrBag_Classifier(CAttrBag_Model &_owner);

		/// initialize the bootstrap sample
		void InitBootstrapCount(int SampCnt[]);
		/// assign the haplotype frequencies
		void Assign(int n_snp, const int snpidx[], const int samp_num[], int n_haplo,
			const double *freq, const int *hla, char *const haplo[], double *_acc=NULL);
		/// grow this classifier by adding SNPs
		void Grow(CBaseSampling &VarSampling, int mtry, bool prune,
			bool verbose, bool verbose_detail, bool debug);

		/// the owner
		inline CAttrBag_Model &Owner() { return *_Owner; }
		/// the number of SNPs
		inline const int nSNP() const { return _SNPIndex.size(); }
		/// the number of haplotyeps
		inline const int nHaplo() const { return _Haplo.TotalNumOfHaplo(); }
		/// the out-of-bag accuracy
		inline const double OutOfBag_Accuracy() const { return _OutOfBag_Accuracy; }
		/// the SNP selection
		inline const vector<int> &SNPIndex() const { return _SNPIndex; }
		/// the bootstrapped individuals
		inline const vector<int> &BootstrapCount() const { return _BootstrapCount; }
		/// the haplotype list
		inline const CHaplotypeList &Haplotype() const { return _Haplo; }

	protected:
		/// the owner
		CAttrBag_Model *_Owner;
		/// the haplotype list
		CHaplotypeList _Haplo;
		/// the bootstrapped individuals
		vector<int> _BootstrapCount;
		/// the SNP selection
		vector<int> _SNPIndex;
		/// the out-of-bag accuracy
		double _OutOfBag_Accuracy;
	};


	/// HIBAG -- the attribute bagging model
	class CAttrBag_Model
	{
	public:
		friend class CAttrBag_Classifier;

		CAttrBag_Model();

		/// initialize the training model
		void InitTraining(int n_snp, int n_samp, int n_hla);
		/// initialize the training model
		void InitTraining(int n_snp, int n_samp, int *snp_geno, int n_hla, int *H1, int *H2);

		/// get a new individual classifier
		CAttrBag_Classifier *NewClassifierBootstrap();
		/// get a new individual classifier
		CAttrBag_Classifier *NewClassifierAllSamp();

		/// build n individual classifiers with the specified parameters
		void BuildClassifiers(int nclassifier, int mtry, bool prune,
			bool verbose, bool verbose_detail=false, bool debug=false);

		/** get the best-guess HLA types
		 *  \param genomat
		 *  \param n_samp
		 *  \param vote_method  1: average posterior prob, 2: majority voting
		 *  \param OutH1
		 *  \param OutH2
		 *  \param OutProb
		 *  \param ShowInfo
		**/
		void PredictHLA(const int *genomat, int n_samp, int vote_method,
			int OutH1[], int OutH2[], double OutProb[], bool ShowInfo);

		/** get the posterior probabilities of HLA type
		 *  \param genomat
		 *  \param n_samp
		 *  \param vote_method  1: average posterior prob, 2: majority voting
		 *  \param OutProb
		 *  \param ShowInfo
		**/
		void PredictHLA_Prob(const int *genomat, int n_samp, int vote_method,
			double OutProb[], bool ShowInfo);

		/// the number of samples
		inline const int nSamp() const { return _SNPMat.Num_Total_Samp; }
		/// the number of SNPs
		inline const int nSNP() const { return _SNPMat.Num_Total_SNP; }
		/// the number of unique HLA alleles
		inline const int nHLA() const { return _HLAList.Num_HLA_Allele(); }

		/// a list of HLA types
		inline const CHLATypeList &HLAList() const
			{ return _HLAList; }
		/// a list of individual classifiers
		inline const vector<CAttrBag_Classifier> &ClassifierList() const
			{ return _ClassifierList; }

	protected:
		/// the SNP genotype matrix
		CSNPGenoMatrix _SNPMat;
		/// a list of HLA alleles
		CHLATypeList _HLAList;
		/// a list of individual classifiers
		vector<CAttrBag_Classifier> _ClassifierList;
		/// variable selection algorithm
		CVariableSelection _VarSelect;
		/// prediction algorithm
		CAlg_Prediction _Predict;

		/// prediction HLA types internally
		void _PredictHLA(const int *geno, const int weights[], int vote_method);
		/// get weight with respect to missing SNPs
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

		/// initialize
		void Init(long TotalCnt, bool ShowInit);
		/// move forward
		bool Forward(long step, bool Show);
		/// show progress information
		virtual void ShowProgress();

        /// Return the current percentile
		inline int Percent() const { return fPercent; }
		/// Return the total number
		inline long Total() const { return fTotal; }
		/// Return the current position
		inline long Current() const { return fCurrent; }

	protected:
		long fTotal;      //< the total number
		long fCurrent;    //< the current number
		int fPercent;     //< the corresponding percent
		clock_t OldTime;  //< the old time point
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

#endif /* _LibHLA_H_ */
