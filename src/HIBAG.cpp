// ===============================================================
//
// HIBAG R package, HLA Genotype Imputation with Attribute Bagging
// Copyright (C) 2011 - 2015   Xiuwen Zheng (zhengx@u.washington.edu)
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

#include <string>
#include <memory>
#include <limits>
#include <map>
#include <algorithm>
#include <fstream>
#include <vector>

#include "LibHLA.h"
#include <R.h>
#include <Rdefines.h>


namespace HLA_LIB
{
	// the last error information
	std::string _LastError;
}


using namespace std;
using namespace HLA_LIB;


#define LongBool int
#define DLLEXPORT


extern "C"
{
	#define CORETRY			try {
	#define CORECATCH(cmd)	} \
		catch (exception &E) { \
			_LastError = E.what(); \
			cmd; \
		} \
		catch (const char *E) { \
			_LastError = E; \
			cmd; \
		}


// ===========================================================
// the private functions
// ===========================================================

/// assign a string
inline static void RStrAgn(const char *Text, char **rstr)
{
	*rstr = R_alloc(strlen(Text)+1, 1);
	if (*rstr == NULL)
		throw "R_alloc return NULL!";
	strcpy(*rstr, Text);
}



// ===========================================================
// the public functions
// ===========================================================

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
//
// The functions for HLA alleles
//

struct TAlleleItem
{
	vector<int> Index;
	vector<string> Idx_Suffix;
	int list_index;
	
	TAlleleItem(const char *str, int _idx)
	{
		string num, suffix;
		bool prefix_num = true;
		while (true)
		{
			char ch = *str++;
			if (prefix_num)
			{
				if (('0' <= ch) && (ch <= '9'))
					num.push_back(ch);
				else
					prefix_num = false;
			}
			if (!prefix_num)
			{
				if ((ch == ':') || (ch == 0))
				{
					int m = numeric_limits<int>::max();
					if (!num.empty()) m = atoi(num.c_str());
					Index.push_back(m);
					Idx_Suffix.push_back(suffix);
					num.clear(); suffix.clear();
					prefix_num = true;
					if (ch == 0) break;
				} else {
					suffix.push_back(ch);
				}
			}
		}

		list_index = _idx;
	}
};

static bool sortfn(const TAlleleItem *I1, const TAlleleItem *I2)
{
	int smin = min((int)I1->Index.size(), (int)I2->Index.size());
	for (int i=0; i < smin; i++)
	{
		if (I1->Index[i] < I2->Index[i])
		{
			return true;
		} else if (I1->Index[i] > I2->Index[i])
		{
			return false;
		} else {
			if (I1->Idx_Suffix[i] < I2->Idx_Suffix[i])
			{
				return true;
			} else if (I1->Idx_Suffix[i] > I2->Idx_Suffix[i])
			{
				return false;
			}		
		}
	}

	return (I1->Index.size() <= I2->Index.size());
}

/**
 *  to sort the HLA alleles
 *
 *  \param n_hla      the number of HLA alleles
 *  \param hlastr     the pointer to HLA allele strings
 *  \param outstr     the pointer to output allele strings
 *  \param out_err    output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void HIBAG_SortAlleleStr(int *n_hla, char *const hlastr[],
	char *outstr[], LongBool *out_err)
{
	CORETRY
		// HLA alleles
		vector<TAlleleItem> HLA;
		for (int i=0; i < *n_hla; i++)
			HLA.push_back(TAlleleItem(hlastr[i], i));

		vector<TAlleleItem*> list;
		for (int i=0; i < *n_hla; i++)
			list.push_back(&HLA[i]);

		// sort
		sort(list.begin(), list.end(), sortfn);

		// output
		for (int i=0; i < *n_hla; i++)
			RStrAgn(hlastr[list[i]->list_index], &outstr[i]);

		*out_err = 0;
	CORECATCH(*out_err = 1)
}



// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
//
//  SNP functions
//

/// to detect and correct strand problem

static inline bool ATGC(const string &s)
{
	return (s=="A") || (s=="T") || (s=="G") || (s=="C");
}
static inline int ALLELE_MINOR(double freq)
{
	return (freq <= 0.5) ? 0 : 1;
}
static inline void split_allele(const char *txt, string &allele1, string &allele2)
{
	const char *p = strchr(txt, '/');
	if (p != NULL)
	{
		// the first allele
		allele1.assign(txt, p);
		for (unsigned int i=0; i < allele1.size(); i++)
			allele1[i] = toupper(allele1[i]);

		// the second allele	
		allele2 = p + 1;
		for (unsigned int i=0; i < allele2.size(); i++)
			allele2[i] = toupper(allele2[i]);
	} else {
		// no a second allele
		allele1 = txt;
		for (unsigned int i=0; i < allele1.size(); i++)
			allele1[i] = toupper(allele1[i]);
		allele2.clear();
	}
}

DLLEXPORT void HIBAG_AlleleStrand(char *allele1[], double afreq1[], int I1[],
	char *allele2[], double afreq2[], int I2[],
	int *if_same_strand, int *n,
	LongBool out_flag[], int *out_n_stand_amb, int *out_n_mismatch,
	LongBool *out_err)
{
	CORETRY
		// initialize: A-T pair, C-G pair
		map<string, string> MAP;
		MAP["A"] = "T"; MAP["C"] = "G"; MAP["G"] = "C"; MAP["T"] = "A";

		*out_n_stand_amb = 0;
		*out_n_mismatch = 0;

		const bool check_strand = (if_same_strand[0] == 0);

		// loop for each SNP
		for (int i=0; i < *n; i++)
		{
			// if true, need switch strand
			bool switch_flag = false;

			// if true, need to compare the allele frequencies
			//   0 -- no switch detect
			//   1 -- detect whether switch or not for stand ambiguity
			//   2 -- detect whether switch or not for mismatching alleles
			int switch_freq_detect = 0;

			// ``ref / nonref alleles''
			string s1, s2;
			string p1, p2;
			split_allele(allele1[I1[i]-1], s1, s2);
			split_allele(allele2[I2[i]-1], p1, p2);

			// allele frequency
			const double F1 = afreq1[I1[i]-1];
			const double F2 = afreq2[I2[i]-1];

			if (ATGC(s1) && ATGC(s2) && ATGC(p1) && ATGC(p2))
			{
				// check
				if ( (s1 == p1) && (s2 == p2) )
				{
					if (check_strand)
					{
						// for example, + C/G <---> - C/G, strand ambi
						if (s1 == MAP[p2])
							switch_freq_detect = 1;
					}
				} else if ( (s1 == p2) && (s2 == p1) )
				{
					if (check_strand)
					{
						// for example, + C/G <---> - G/C, strand ambi
						if (s1 == MAP[p1])
							switch_freq_detect = 1;
						else
							switch_flag = true;
					} else
						switch_flag = true;
				} else {
					if (check_strand)
					{
						if ( (s1 == MAP[p1]) && (s2 == MAP[p2]) )
						{
							// for example, + C/G <---> - G/C, strand ambi
							if (s1 == p2)
								switch_freq_detect = 1;
						} else if ( (s1 == MAP[p2]) && (s2 == MAP[p1]) )
							switch_flag = true;
						else
							switch_freq_detect = 2;
					} else
						switch_freq_detect = 2;
				}
			} else {
				if ((s1 == p1) && (s2 == p2))
				{
					if (s1 == s2)
						switch_freq_detect = 1;  // ambiguous
				} else if ((s1 == p2) && (s2 == p1))
				{
					if (s1 == s2)
						switch_freq_detect = 1;  // ambiguous
					else
						switch_flag = true;
				} else
					switch_freq_detect = 2;
			}

			if (switch_freq_detect != 0)
			{
				switch_flag = (ALLELE_MINOR(F1) != ALLELE_MINOR(F2));
				if (switch_freq_detect == 1)
					(*out_n_stand_amb) ++;
				else
					(*out_n_mismatch) ++;
			}

			out_flag[i] = switch_flag;
		}

		*out_err = 0;

	CORECATCH(*out_err = 1)
}



// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
//
// HIBAG: Attribute Bagging method
//

/// the number of attribute bagging models allowed
#define MODEL_NUM_LIMIT 256
/// the model container
static CAttrBag_Model* _HIBAG_MODELS_[MODEL_NUM_LIMIT];

/// need a new model
static int _Need_New_HIBAG_Model()
{
	for (int i=0; i < MODEL_NUM_LIMIT; i++)
		if (_HIBAG_MODELS_[i] == NULL) return i;
	throw ErrHLA("No memory space to store a new HIBAG model, "
		"please call \"hlaClose\" to release unused HIBAG models.");
}

/// check the model
static void _Check_HIBAG_Model(int model)
{
	if ((0 <= model) && (model < MODEL_NUM_LIMIT))
	{
		if (_HIBAG_MODELS_[model] != NULL)
			return;
	}
	throw ErrHLA("The handle of HIBAG model has been closed.");
}


/**
 *  to build a HIBAG model
 *
 *  \param nSNP       the number of SNPs
 *  \param nSamp      the number of samples
 *  \param nHLA       the number of different HLA alleles
 *  \param out_Model  output the model index
 *  \param out_err    output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void HIBAG_New(int *nSamp, int *nSNP, int *nHLA, int *out_Model,
	LongBool *out_err)
{
	CORETRY
		int model = _Need_New_HIBAG_Model();
		_HIBAG_MODELS_[model] = new CAttrBag_Model;
		_HIBAG_MODELS_[model]->InitTraining(*nSNP, *nSamp, *nHLA);
		*out_Model = model;
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/**
 *  to build a HIBAG model
 *  \param nSNP       the number of SNPs
 *  \param nSamp      the number of samples
 *  \param snp_geno   the SNP genotypes, (0 -- BB, 1 -- AB, 2 -- AA, other -- missing value)
 *  \param nHLA       the number of different HLA alleles
 *  \param H1         the first HLA allele of a HLA type
 *  \param H2         the second HLA allele of a HLA type
 *  \param out_Model  output the model index
 *  \param out_err    output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void HIBAG_Training(int *nSNP, int *nSamp, int *snp_geno, int *nHLA,
	int *H1, int *H2, int *out_Model, LongBool *out_err)
{
	CORETRY
		int model = _Need_New_HIBAG_Model();
		_HIBAG_MODELS_[model] = new CAttrBag_Model;
		_HIBAG_MODELS_[model]->InitTraining(*nSNP, *nSamp, snp_geno, *nHLA, H1, H2);
		*out_Model = model;
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/**
 *  to close an existing HIBAG model 
 *
 *  \param model      the model index
 *  \param out_err    output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void HIBAG_Close(int *model, LongBool *out_err)
{
	CORETRY
		_Check_HIBAG_Model(*model);
		CAttrBag_Model* m = _HIBAG_MODELS_[*model];
		_HIBAG_MODELS_[*model] = NULL;
		delete m;
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/**
 *  to add individual classifiers
 *
 *  \param model           the model index
 *  \param nclassifier     the total number of individual classifiers to be created
 *  \param mtry            the number of variables randomly sampled as candidates for selection
 *  \param prune           if TRUE, perform a parsimonious forward variable selection
 *  \param verbose         show information if TRUE
 *  \param verbose_detail  show more information if TRUE
 *  \param out_err         output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void HIBAG_NewClassifiers(int *model, int *nclassifier, int *mtry,
	LongBool *prune, LongBool *verbose, LongBool *verbose_detail,
	LongBool *out_err)
{
	GetRNGstate();
	CORETRY
		_Check_HIBAG_Model(*model);
		_HIBAG_MODELS_[*model]->BuildClassifiers(*nclassifier, *mtry,
			*prune, *verbose, *verbose_detail);
		*out_err = 0;
	CORECATCH(*out_err = 1)
	PutRNGstate();
}


/**
 *  to predict HLA types
 *
 *  \param model        the model index
 *  \param GenoMat      the pointer to the SNP genotypes
 *  \param nSamp        the number of samples in GenoMat
 *  \param out_H1       the first HLA alleles
 *  \param out_H2       the second HLA alleles
 *  \param out_Prob     the posterior probabilities
 *  \param out_err      output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void HIBAG_Predict_Resp(int *model, int *GenoMat, int *nSamp,
	int *vote_method, LongBool *ShowInfo,
	int out_H1[], int out_H2[], double out_Prob[], LongBool *out_err)
{
	CORETRY
		_Check_HIBAG_Model(*model);
		CAttrBag_Model &M = *_HIBAG_MODELS_[*model];

	#if (HIBAG_FLOAT_TYPE_ID == 0)

		M.PredictHLA(GenoMat, *nSamp, *vote_method, out_H1, out_H2, out_Prob,
			NULL, *ShowInfo);

	#elif (HIBAG_FLOAT_TYPE_ID == 1)

		vector<float> tmp(*nSamp);
		M.PredictHLA(GenoMat, *nSamp, *vote_method, out_H1, out_H2, &tmp[0],
			NULL, *ShowInfo);
		for (int i=0; i < *nSamp; i++) out_Prob[i] = tmp[i];

	#else
	#  error "Invalid HIBAG_FLOAT_TYPE_ID"
	#endif

		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/**
 *  to predict HLA types, output posterior probabilities
 *
 *  \param model        the model index
 *  \param GenoMat      the pointer to the SNP genotypes
 *  \param nSamp        the number of samples in GenoMat
 *  \param out_H1       the first HLA alleles
 *  \param out_H2       the second HLA alleles
 *  \param out_Prob     the posterior probabilities
 *  \param out_err      output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void HIBAG_Predict_Prob(int *model, int *GenoMat, int *nSamp,
	int *vote_method, LongBool *ShowInfo, double out_Prob[], LongBool *out_err)
{
	CORETRY
		_Check_HIBAG_Model(*model);
		CAttrBag_Model &M = *_HIBAG_MODELS_[*model];

	#if (HIBAG_FLOAT_TYPE_ID == 0)

		M.PredictHLA_Prob(GenoMat, *nSamp, *vote_method, out_Prob, *ShowInfo);

	#elif (HIBAG_FLOAT_TYPE_ID == 1)

		const int n = M.nHLA()*(M.nHLA()+1) / 2;
		vector<float> tmp(n);
		M.PredictHLA_Prob(GenoMat, *nSamp, *vote_method, &tmp[0], *ShowInfo);
		for (int i=0; i < n; i++) out_Prob[i] = tmp[i];

	#else
	#  error "Invalid HIBAG_FLOAT_TYPE_ID"
	#endif

		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/**
 *  to predict HLA types
 *
 *  \param model        the model index
 *  \param GenoMat      the pointer to the SNP genotypes
 *  \param nSamp        the number of samples in GenoMat
 *  \param out_H1       the first HLA alleles
 *  \param out_H2       the second HLA alleles
 *  \param out_Prob     the posterior probabilities
 *  \param out_err      output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void HIBAG_Predict_Resp_Prob(int *model, int *GenoMat, int *nSamp,
	int *vote_method, LongBool *ShowInfo,
	int out_H1[], int out_H2[], double out_MaxProb[],
	double out_Prob[], LongBool *out_err)
{
	CORETRY
		_Check_HIBAG_Model(*model);
		CAttrBag_Model &M = *_HIBAG_MODELS_[*model];

	#if (HIBAG_FLOAT_TYPE_ID == 0)

		M.PredictHLA(GenoMat, *nSamp, *vote_method, out_H1, out_H2,
			out_MaxProb, out_Prob, *ShowInfo);

	#elif (HIBAG_FLOAT_TYPE_ID == 1)

		vector<float> tmp(*nSamp);
		const int n = M.nHLA()*(M.nHLA()+1) / 2;
		vector<float> tmp_d(n);
		M.PredictHLA(GenoMat, *nSamp, *vote_method, out_H1, out_H2, &tmp[0],
			&tmp_d[0], *ShowInfo);
		for (int i=0; i < *nSamp; i++) out_Prob[i] = tmp[i];
		for (int i=0; i < n; i++) out_Prob[i] = tmp_d[i];

	#else
	#  error "Invalid HIBAG_FLOAT_TYPE_ID"
	#endif

		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/**
 *  to create a new individual classifier with specified parameters
 *
 *  \param model        the model index
 *  \param n_snp        the number of selected SNP markers
 *  \param snpidx       the indices of SNP markers
 *  \param n_haplo      the number of haplotypes
 *  \param freq         the haplotype frequencies
 *  \param hla          the HLA alleles corresponding to the hapltype list
 *  \param haplo        the vector of characters specifying the SNP haplotype list
 *  \param acc          the out-of-bag accuracy
 *  \param out_err      output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void HIBAG_NewClassifierHaplo(int *model, int *n_snp, int snpidx[],
	int samp_num[], int *n_haplo, double *freq, int *hla, char *haplo[],
	double *acc, LongBool *out_err)
{
	CORETRY
		_Check_HIBAG_Model(*model);
		CAttrBag_Classifier *I = _HIBAG_MODELS_[*model]->NewClassifierAllSamp();

	#if (HIBAG_FLOAT_TYPE_ID == 0)

		I->Assign(*n_snp, snpidx, samp_num, *n_haplo, freq, hla, haplo, acc);

	#elif (HIBAG_FLOAT_TYPE_ID == 1)

		float acc_float = *acc;
		vector<float> tmp(*n_haplo);
		for (int i=0; i < *n_haplo; i++) tmp[i] = freq[i];
		I->Assign(*n_snp, snpidx, samp_num, *n_haplo, &tmp[0], hla, haplo,
			&acc_float);

	#else
	#  error "Invalid HIBAG_FLOAT_TYPE_ID"
	#endif

		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/**
 *  to get the number of individual component classifiers
 *
 *  \param model        the model index
 *  \param out_Num      output the number of individual classifiers
 *  \param out_err      output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void HIBAG_GetNumClassifiers(int *model, int *out_Num, LongBool *out_err)
{
	CORETRY
		_Check_HIBAG_Model(*model);
		*out_Num = _HIBAG_MODELS_[*model]->ClassifierList().size();
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/**
 *  to get the number of haplotypes and selected SNP markers in a specified individual classifier
 *
 *  \param model         the model index
 *  \param idx           the index of individual classifier
 *  \param out_NumHaplo  output the number of haplotypes
 *  \param out_NumSNP    output the number of selected SNP markers
 *  \param out_err       output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void HIBAG_Idv_GetNumHaplo(int *model, int *idx,
	int *out_NumHaplo, int *out_NumSNP, LongBool *out_err)
{
	CORETRY
		_Check_HIBAG_Model(*model);
		CAttrBag_Model *AB = _HIBAG_MODELS_[*model];
		*out_NumHaplo = AB->ClassifierList()[*idx - 1].nHaplo();
		*out_NumSNP   = AB->ClassifierList()[*idx - 1].nSNP();
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/**
 *  to get the details of a specified individual classifier
 *
 *  \param model         the model index
 *  \param idx           the index of individual classifier
 *  \param out_freq      output the haplotype frequencies
 *  \param out_hla       output the HLA alleles
 *  \param out_haplo     output the haplotype list
 *  \param out_snpidx    output the indices of SNP markers
 *  \param out_acc       output the out-of-bag accuracy
 *  \param out_err       output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void HIBAG_Classifier_GetHaplos(int *model, int *idx,
	double out_freq[], int out_hla[], char *out_haplo[], int out_snpidx[], int out_samp_num[],
	double *out_acc, LongBool *out_err)
{
	CORETRY
		_Check_HIBAG_Model(*model);
		CAttrBag_Model *AB = _HIBAG_MODELS_[*model];

		const CAttrBag_Classifier &Voter = AB->ClassifierList()[*idx - 1];
		const vector< vector<THaplotype> > &List = Voter.Haplotype().List;

		int idx = 0;
		for (int i=0; i < (int)List.size(); i++)
		{
			vector<THaplotype>::const_iterator it;
			for (it=List[i].begin(); it != List[i].end(); it++)
			{
				out_freq[idx] = it->Frequency;
				out_hla[idx] = i + 1;
				RStrAgn(it->HaploToStr(Voter.nSNP()).c_str(), &out_haplo[idx]);
				idx ++;
			}
		}
		for (int i=0; i < (int)Voter.SNPIndex().size(); i++)
			out_snpidx[i] = Voter.SNPIndex()[i] + 1;
		const vector<int> &N = Voter.BootstrapCount();
		for (int i=0; i < (int)N.size(); i++)
			out_samp_num[i] = N[i];

		*out_acc = Voter.OutOfBag_Accuracy();
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/**
 *  to estimate the confusion matrix
 *
 *  \param n_hla         the number of different HLA alleles
 *  \param init_mat      the initial confusion matrix without any ambiguous state
 *  \param n_DConfusion  the number of double confusions
 *  \param D_mat
 *  \param out_mat       the output confusion matrix
 *  \param out_err       output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void HIBAG_Confusion(int *n_hla, double *init_mat,
	int *n_DConfusion, int *D_mat, double *out_mat, double *tmp_mat,
	LongBool *out_err)
{
	// the max number of iterations
	const int N_MAX_ITERATION = 100;

	CORETRY
		const int nHLA = *n_hla;
		const int SIZE_MAT = sizeof(double)*nHLA*(nHLA+1);
		#define INDEX(T, P, var) var[(nHLA+1)*T + P]

		// initial values
		memcpy(out_mat, init_mat, SIZE_MAT);
		for (int i=0; i < *n_DConfusion; i++)
		{
			int *T = &D_mat[i*4], *P = &D_mat[i*4+2];
			INDEX(T[0], P[0], out_mat) += 0.5;
			INDEX(T[0], P[1], out_mat) += 0.5;
			INDEX(T[1], P[0], out_mat) += 0.5;
			INDEX(T[1], P[1], out_mat) += 0.5;
		}

		// EM update ...
		for (int iter=0; iter < N_MAX_ITERATION; iter++)
		{
			double f1, f2, s;
			// copy the current probabilities to the old ones
			memcpy(tmp_mat, out_mat, SIZE_MAT);
			// update ...
			memcpy(out_mat, init_mat, SIZE_MAT);
			for (int i=0; i < *n_DConfusion; i++)
			{
				int *T = &D_mat[i*4], *P = &D_mat[i*4+2];

				f1 = INDEX(T[0], P[0], tmp_mat);
				f2 = INDEX(T[0], P[1], tmp_mat);
				s = 1.0 / (f1 + f2);
				INDEX(T[0], P[0], out_mat) += f1 * s;
				INDEX(T[0], P[1], out_mat) += f2 * s;

				f1 = INDEX(T[1], P[0], tmp_mat);
				f2 = INDEX(T[1], P[1], tmp_mat);
				s = 1.0 / (f1 + f2);
				INDEX(T[1], P[0], out_mat) += f1 * s;
				INDEX(T[1], P[1], out_mat) += f2 * s;
			}
		}

		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/**
 *  to get the algorithm parameters
 *
 *  \param EM_MaxNum        The max number of iterations for EM algorithm
 *  \param EM_RelTol        The reltol convergence tolerance
 *  \param Search_MaxNum    The max number of SNP markers in an individual classifier
**/
DLLEXPORT void HIBAG_GetParam(int *EM_MaxNum, double *EM_RelTol, int *Search_MaxNum)
{
	*EM_MaxNum = EM_MaxNum_Iterations;
	*EM_RelTol = EM_FuncRelTol;
	*Search_MaxNum = HIBAG_MAXNUM_SNP_IN_CLASSIFIER;
}


/**
 *  to set the algorithm parameters
 *
 *  \param EM_MaxNum        The max number of iterations for EM algorithm
 *  \param EM_RelTol        The reltol convergence tolerance
 *  \param Search_MaxNum    The max number of SNP markers in an individual classifier
**/
DLLEXPORT void HIBAG_SetParam(int *EM_MaxNum, LongBool *If_EM_MaxNum,
	double *EM_RelTol, LongBool *If_EM_RelTol,
	int *Search_MaxNum, LongBool *If_Search_MaxNum)
{
	if (*If_EM_MaxNum)
		EM_MaxNum_Iterations = *EM_MaxNum;
	if (*If_EM_RelTol)
		EM_FuncRelTol = *EM_RelTol;
//	if (*If_Search_MaxNum)
//		SEARCH_SNP_MAXNUM = *Search_MaxNum;
}


/**
 *  to detect the storage mode of a PLINK BED file
 *
 *  \param bedfn         the file name of PLINK BED file
 *  \param out_SNPOrder  output the storage mode
 *  \param out_err       output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void HIBAG_BEDFlag(char **bedfn, int *out_SNPOrder, LongBool *out_err)
{
	CORETRY
		ifstream file(*bedfn, ios::binary);
		if (!file.good())
			throw ErrHLA("Cannot open the file %s.", *bedfn);
		char prefix[3];
		file.read(prefix, 3);
		if ((prefix[0] != 0x6C) || (prefix[1] != 0x1B))
			throw ErrHLA("Invalid prefix in the bed file.");
		*out_SNPOrder = (unsigned char)(prefix[2]);
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/**
 *  to convert from a PLINK BED file
 *
 *  \param bedfn         the file name of PLINK BED file
 *  \param n_samp        the number of samples
 *  \param n_snp         the number of SNPs
 *  \param mode          the storage mode
 *  \param verbose       if TRUE, show information
 *  \param out_geno      output genotypes
 *  \param out_err       output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void HIBAG_ConvBED(char **bedfn, int *n_samp, int *n_snp, int *n_save_snp,
	LongBool *mode, LongBool *snp_flag, LongBool *verbose, int *out_geno,
	LongBool *out_err)
{
	CORETRY
		// open file
		ifstream file(*bedfn, ios::binary);
		if (!file.good())
			throw ErrHLA("Fail to open the file \"%s\".", *bedfn);

		// read prefix
		char prefix[3];
		file.read(prefix, 3);

		// determine the values of packed genotypes
		int nRe, nPack, nNumPack, nNum;
		if (*mode)
		{
			// the individual-major mode
			nRe = (*n_snp) % 4; nNumPack = (*n_snp)/4;
			nPack = (nRe > 0) ? (nNumPack + 1) : nNumPack;
			nNum = (*n_samp);
		} else {
			// the SNP-major mode
			nRe = (*n_samp) % 4; nNumPack = (*n_samp)/4;
			nPack = (nRe > 0) ? (nNumPack + 1) : nNumPack;
			nNum = (*n_snp);
		}

		vector<char> srcgeno(nPack);
		vector<int> dstgeno((nNumPack+1) * 4);
		static const int cvt[4] = { 2, INT_MIN, 1, 0 };
		int I_SNP = 0;

		// for - loop
		for (int i=0; i < nNum; i++)
		{
			// read genotypes
			file.read(&srcgeno[0], nPack);
			// unpacked
			int *p = &dstgeno[0];
			for (int k=0; k < nNumPack; k++)
			{
				unsigned char g = srcgeno[k];
				*p++ = cvt[g & 0x03]; g >>= 2;
				*p++ = cvt[g & 0x03]; g >>= 2;
				*p++ = cvt[g & 0x03]; g >>= 2;
				*p++ = cvt[g & 0x03];
			}
			if (nRe > 0)
			{
				unsigned char g = srcgeno[nNumPack];
				for (long k=0; k < nRe; k++)
				{
					*p++ = cvt[g & 0x03]; g >>= 2;
				}
			}

			// write
			if (*mode)
			{
				// the individual-major mode
				int *pI = (out_geno + i * (*n_save_snp));
				for (int j=0; j < *n_snp; j++)
				{
					if (snp_flag[j])
						*pI++ = dstgeno[j];
				}
			} else {
				// the SNP-major mode
				if (snp_flag[i])
				{
					int *pI = (out_geno + I_SNP);
					I_SNP ++;
					for (int j=0; j < *n_samp; j++)
					{
						*pI = dstgeno[j];
						pI += *n_save_snp;
					}
				}
			}
		}
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
/**
 *  to get an error message
 *  \param Msg           the last error information
**/
DLLEXPORT SEXP HIBAG_ErrMsg()
{
	return mkString(_LastError.c_str());
}



// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
//
// HIBAG: an attribute bagging method
//

/// initialize the package
DLLEXPORT SEXP HIBAG_Init()
{
	memset((void*)_HIBAG_MODELS_, 0, sizeof(_HIBAG_MODELS_));

	SEXP ans;
	#ifdef HIBAG_SSE_OPTIMIZE_HAMMING_DISTANCE
	#   ifdef __SSE4_2__
		ans = ScalarInteger(2);
	#   else
		ans = ScalarInteger(1);
	#   endif
	#else
		ans = ScalarInteger(0);
	#endif

	return ans;
}

/// finalize the package
DLLEXPORT SEXP HIBAG_Done()
{
	try {
		for (int i=0; i < MODEL_NUM_LIMIT; i++)
		{
			CAttrBag_Model* m = _HIBAG_MODELS_[i];
			_HIBAG_MODELS_[i] = NULL;
			try {
				delete m;
			} catch(...) {}
		}
	} catch(...) {}

	return R_NilValue;
}

} // extern "C"
