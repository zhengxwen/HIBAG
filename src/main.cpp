// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// main.cpp: HLA Genotype Imputation with Attribute Bagging
//
// Copyright (C) 2012	Xiuwen Zheng
//
// This file is part of HIBAG package.
//
// HIBAG is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License Version 3 as
// published by the Free Software Foundation.
//
// HIBAG is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with CoreArray.
// If not, see <http://www.gnu.org/licenses/>.

#include <R.h>
#include <string>
#include <memory>
#include <limits>
#include <map>
#include <algorithm>
#include <fstream>
#include <LibHLA.h>

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
	
	TAlleleItem(const char *str)
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
	}
};

static bool sortfn(const pair<TAlleleItem*, int> &I1, const pair<TAlleleItem*, int> &I2)
{
	const TAlleleItem &p1 = *I1.first;
	const TAlleleItem &p2 = *I2.first;

	int smin = min((int)p1.Index.size(), (int)p2.Index.size());
	for (int i=0; i < smin; i++)
	{
		if (p1.Index[i] < p2.Index[i])
		{
			return true;
		} else if (p1.Index[i] > p2.Index[i])
		{
			return false;
		} else {
			if (p1.Idx_Suffix[i] < p2.Idx_Suffix[i])
			{
				return true;
			} else if (p1.Idx_Suffix[i] > p2.Idx_Suffix[i])
			{
				return false;
			}		
		}
	}

	return (p1.Index.size() <= p2.Index.size());
}

/**
 *  to sort the HLA alleles
 *
 *  \param n_hla      the number of HLA alleles
 *  \param hlastr     the pointer to HLA allele strings
 *  \param outstr     the pointer to output allele strings
 *  \param out_err    output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void hlaSortAlleleStr(int *n_hla, char *const hlastr[], char *outstr[],
	LongBool *out_err)
{
	CORETRY
		// HLA alleles
		vector<TAlleleItem> HLA;
		for (int i=0; i < *n_hla; i++)
			HLA.push_back(TAlleleItem(hlastr[i]));

		vector< pair<TAlleleItem*, int> > list;
		for (int i=0; i < *n_hla; i++)
			list.push_back(pair<TAlleleItem*, int>(&HLA[i], i));

		// sort
		sort(list.begin(), list.end(), sortfn);

		// output
		for (int i=0; i < *n_hla; i++)
			RStrAgn(hlastr[list[i].second], &outstr[i]);

		*out_err = 0;
	CORECATCH(*out_err = 1)
}



// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
//
//  SNP functions
//

/// to detect strand problem

static inline bool ATGC(char ch)
{
	return (ch=='A') || (ch=='T') || (ch=='G') || (ch=='C');
}
static inline bool ATGC_atgc(char ch)
{
	return (ch=='A') || (ch=='T') || (ch=='G') || (ch=='C') ||
		(ch=='a') || (ch=='t') || (ch=='g') || (ch=='c');
}
static inline int ALLELE_MINOR(double freq)
{
	return (freq <= 0.5) ? 0 : 1;
}
static inline void split_allele(const char *txt, char &allele1, char &allele2)
{
	const char *p = strchr(txt, '/');
	if (p != NULL)
	{
		// the first allele
		if ((p == (txt+1)) && ATGC_atgc(*txt))
			allele1 = toupper(*txt);
		else
			allele1 = '.';

		// the second allele	
		txt = p + 1;
		if ((strlen(txt) == 1) && ATGC_atgc(*txt))
			allele2 = toupper(*txt);
		else
			allele2 = '.';
	} else {
		// no a second allele
		if ((strlen(txt) == 1) && ATGC_atgc(*txt))
			allele1 = toupper(*txt);
		else
			allele1 = '.';
		allele2 = '.';
	}
}

DLLEXPORT void hlaAlleleStrand(char *allele1[], double afreq1[], int I1[],
	char *allele2[], double afreq2[], int I2[], int *n,
	LongBool out_flag[], int *out_n_stand_amb, int *out_n_mismatch,
	LongBool *out_err)
{
	CORETRY
		// initialize: A-T pair, C-G pair
		map<char, char> MAP;
		MAP['A'] = 'T'; MAP['C'] = 'G'; MAP['G'] = 'C'; MAP['T'] = 'A';

		*out_n_stand_amb = 0;
		*out_n_mismatch = 0;

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
			char s1, s2;
			char p1, p2;
			split_allele(allele1[I1[i]-1], s1, s2);
			split_allele(allele2[I2[i]-1], p1, p2);

			// allele frequency
			double F1=afreq1[I1[i]-1], F2=afreq2[I2[i]-1];

			if (ATGC(s1) && ATGC(s2) && ATGC(p1) && ATGC(p2))
			{
				// check
				if ( (s1 == p1) && (s2 == p2) )
				{
					// for example, + C/G <---> - C/G, strand ambi
					if (s1 == MAP[p2])
						switch_freq_detect = 1;
				} else if ( (s1 == p2) && (s2 == p1) )
				{
					// for example, + C/G <---> - G/C, strand ambi
					if (s1 == MAP[p1])
						switch_freq_detect = 1;
					else
						switch_flag = true;
				} else if ( (s1 == MAP[p1]) && (s2 == MAP[p2]) )
				{
					// for example, + C/G <---> - G/C, strand ambi
					if (s1 == p2)
						switch_freq_detect = 1;
				} else if ( (s1 == MAP[p2]) && (s2 == MAP[p1]) )
				{
					switch_flag = true;
				} else {
					switch_freq_detect = 2;
					// throw ErrHLA("Invalid strand in SNP %d: %c/%c <--> %c/%c",
					//	i+1, s1, s2, p1, p2);
				}
			} else {
				switch_freq_detect = 2;
				// throw ErrHLA("Invalid alleles in the sample: %d", i+1);
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
static CAttrBag_Model* _AB_Model_[MODEL_NUM_LIMIT];

/// need a new model
static int _Need_New_AB_Model()
{
	for (int i=0; i < MODEL_NUM_LIMIT; i++)
		if (_AB_Model_[i] == NULL) return i;
	throw ErrHLA("No memory space to store a new attribute bagging model, "
		"please call \"hlaClose\" to release unused models.");
}

/// check the model
static void _Check_Model(int model)
{
	if ((0 <= model) && (model < MODEL_NUM_LIMIT))
	{
		if (_AB_Model_[model] != NULL)
			return;
	}
	throw ErrHLA("The attribute bagging model has been closed.");
}


/**
 *  to build an attribute bagging (AB) model
 *
 *  \param nSNP       the number of SNPs
 *  \param nSamp      the number of samples
 *  \param nHLA       the number of different HLA alleles
 *  \param out_Model  output the model index
 *  \param out_err    output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void hlaAB_Model_New(int *nSamp, int *nSNP, int *nHLA, int *out_Model, LongBool *out_err)
{
	CORETRY
		int model = _Need_New_AB_Model();
		_AB_Model_[model] = new CAttrBag_Model;
		_AB_Model_[model]->InitTraining(*nSNP, *nSamp, *nHLA);
		*out_Model = model;
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/**
 *  to build an attribute bagging (AB) model
 *  \param nSNP       the number of SNPs
 *  \param nSamp      the number of samples
 *  \param snp_geno   the SNP genotypes, (0 -- BB, 1 -- AB, 2 -- AA, other -- missing value)
 *  \param nHLA       the number of different HLA alleles
 *  \param H1         the first HLA allele of a HLA type
 *  \param H2         the second HLA allele of a HLA type
 *  \param out_Model  output the model index
 *  \param out_err    output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void hlaAB_Model_Training(int *nSNP, int *nSamp, int *snp_geno, int *nHLA,
	int *H1, int *H2, int *out_Model, LongBool *out_err)
{
	CORETRY
		int model = _Need_New_AB_Model();
		_AB_Model_[model] = new CAttrBag_Model;
		_AB_Model_[model]->InitTraining(*nSNP, *nSamp, snp_geno, *nHLA, H1, H2);
		*out_Model = model;
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/**
 *  to close the specified attribute bagging model
 *
 *  \param model      the model index
 *  \param out_err    output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void hlaAB_Close(int *model, LongBool *out_err)
{
	CORETRY
		_Check_Model(*model);
		CAttrBag_Model* m = _AB_Model_[*model];
		_AB_Model_[*model] = NULL;
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
 *  \param verbose         show information if TRUE
 *  \param verbose_detail  show more information if TRUE
 *  \param out_err         output the error information, 0 -- no error, 1 -- an error exists
**/
DLLEXPORT void hlaAB_NewClassifiers(int *model, int *nclassifier, int *mtry,
	double *var_prob, LongBool *var_prob_flag,
	LongBool *prune, LongBool *verbose, LongBool *verbose_detail, LongBool *Debug,
	LongBool *out_err)
{
	CORETRY
		_Check_Model(*model);
		_AB_Model_[*model]->RunNewClassifiers(*nclassifier, *mtry,
			(*var_prob_flag) ? var_prob : NULL,
			*prune, *verbose, *verbose_detail, *Debug);
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
DLLEXPORT void hlaAB_Predict(int *model, int *GenoMat, int *nSamp,
	LongBool *ShowInfo, int out_H1[], int out_H2[], double out_Prob[], LongBool *out_err)
{
	CORETRY
		_Check_Model(*model);
		_AB_Model_[*model]->PredictHLA(GenoMat, *nSamp, out_H1, out_H2, out_Prob, *ShowInfo);
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
DLLEXPORT void hlaAB_Predict_Prob(int *model, int *GenoMat, int *nSamp,
	LongBool *ShowInfo, double out_Prob[], LongBool *out_err)
{
	CORETRY
		_Check_Model(*model);
		_AB_Model_[*model]->PredictHLA_Prob(GenoMat, *nSamp, out_Prob, *ShowInfo);
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
DLLEXPORT void hlaAB_NewClassifierHaplo(int *model, int *n_snp, int snpidx[], int samp_num[],
	int *n_haplo, double *freq, int *hla, char *haplo[], double *acc, LongBool *out_err)
{
	CORETRY
		_Check_Model(*model);
		CAttrBag_Classifier *I = _AB_Model_[*model]->NewClassifierAllSamp();
		I->Assign(*n_snp, snpidx, samp_num, *n_haplo, freq, hla, haplo, acc);
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
DLLEXPORT void hlaAB_GetNumClassifiers(int *model, int *out_Num, LongBool *out_err)
{
	CORETRY
		_Check_Model(*model);
		*out_Num = _AB_Model_[*model]->ClassifierList().size();
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
DLLEXPORT void hlaAB_Idv_GetNumHaplo(int *model, int *idx,
	int *out_NumHaplo, int *out_NumSNP, LongBool *out_err)
{
	CORETRY
		_Check_Model(*model);
		CAttrBag_Model *AB = _AB_Model_[*model];
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
DLLEXPORT void hlaAB_Tree_GetHaplos(int *model, int *idx,
	double out_freq[], int out_hla[], char *out_haplo[], int out_snpidx[], int out_samp_num[],
	double *out_acc, LongBool *out_err)
{
	CORETRY
		_Check_Model(*model);
		CAttrBag_Model *AB = _AB_Model_[*model];

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
				RStrAgn(it->SNPToString(Voter.nSNP()).c_str(), &out_haplo[idx]);
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
DLLEXPORT void hlaConfusion(int *n_hla, double *init_mat,
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
DLLEXPORT void hlaAB_GetParam(int *EM_MaxNum, double *EM_RelTol, int *Search_MaxNum)
{
	*EM_MaxNum = EM_MaxNum_Iterations;
	*EM_RelTol = EM_FuncRelTol;
	*Search_MaxNum = MAXNUM_SNP_IN_CLASSIFIER;
}


/**
 *  to set the algorithm parameters
 *
 *  \param EM_MaxNum        The max number of iterations for EM algorithm
 *  \param EM_RelTol        The reltol convergence tolerance
 *  \param Search_MaxNum    The max number of SNP markers in an individual classifier
**/
DLLEXPORT void hlaAB_SetParam(int *EM_MaxNum, LongBool *If_EM_MaxNum,
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
DLLEXPORT void hlaBEDFlag(char **bedfn, int *out_SNPOrder, LongBool *out_err)
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
DLLEXPORT void hlaConvBED(char **bedfn, int *n_samp, int *n_snp, int *n_save_snp,
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

		auto_ptr<char> srcgeno(new char[nPack]);
		auto_ptr<int> dstgeno(new int[(nNumPack+1)*4]);
		static const int cvt[4] = { 2, INT_MIN, 1, 0 };
		int I_SNP = 0;

		// for - loop
		for (int i=0; i < nNum; i++)
		{
			// read genotypes
			file.read(srcgeno.get(), nPack);
			// unpacked
			int *p = dstgeno.get();
			for (int k=0; k < nNumPack; k++)
			{
				unsigned char g = srcgeno.get()[k];
				*p++ = cvt[g & 0x03]; g >>= 2;
				*p++ = cvt[g & 0x03]; g >>= 2;
				*p++ = cvt[g & 0x03]; g >>= 2;
				*p++ = cvt[g & 0x03];
			}
			if (nRe > 0)
			{
				unsigned char g = srcgeno.get()[nNumPack];
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
						{ *pI++ = dstgeno.get()[j]; }
				}
			} else {
				// the SNP-major mode
				if (snp_flag[i])
				{
					int *pI = (out_geno + I_SNP);
					I_SNP ++;
					for (int j=0; j < *n_samp; j++)
					{
						*pI = dstgeno.get()[j];
						pI += *n_save_snp;
					}
				}
			}
		}
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/**
 *  to get an error message
 *
 *  \param Msg           the last error information
**/
DLLEXPORT void hlaErrMsg(char **Msg)
{
	RStrAgn(_LastError.c_str(), Msg);
}



// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
//
// Attribute Bagging (AB) method
//

/// initialize the package
DLLEXPORT void hlaInit()
{
	memset((void*)_AB_Model_, 0, sizeof(_AB_Model_));
}

/// finalize the package
DLLEXPORT void hlaDone()
{
	try {
		for (int i=0; i < MODEL_NUM_LIMIT; i++)
		{
			CAttrBag_Model* m = _AB_Model_[i];
			_AB_Model_[i] = NULL;
			try {
				delete m;
			} catch(...) {}
		}
	} catch(...) {}
}

} // extern "C"
