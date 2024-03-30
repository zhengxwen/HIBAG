// ===============================================================
//
// HIBAG R package (HLA Genotype Imputation with Attribute Bagging)
// Copyright (C) 2011-2024   Xiuwen Zheng (zhengx@u.washington.edu)
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

#include "LibHLA.h"
#include <string>
#include <memory>
#include <limits>
#include <map>
#include <algorithm>
#include <fstream>
#include <vector>

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <RcppParallel.h>
#include <tbb/parallel_for.h>

using namespace std;
using namespace HLA_LIB;


extern "C"
{
/// try block
#define CORE_TRY    \
	bool has_error = false; \
	SEXP rv_ans = R_NilValue; \
	try {

/// catch block
#define CORE_CATCH    \
	} \
	catch (exception &E) { \
		_LastError = E.what(); has_error = true; \
	} \
	catch (const char *E) { \
		_LastError = E; has_error = true; \
	} \
	catch (...) { \
		_LastError = "unknown error!"; has_error = true; \
	} \
	if (has_error) Rf_error("%s", _LastError.c_str()); \
	return rv_ans;


/// the last error information
static std::string _LastError;

static SEXP hibag_data_frame = R_NilValue;
static SEXP hibag_clr_nm = R_NilValue;
static SEXP hibag_clr_haplo_nm = R_NilValue;


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
	/// 1-field --> 2-digit, 2-field --> 4-digit
	vector<int> Field;
	/// the suffix for each field
	vector<string> Field_Suffix;
	/// the order in the container (starting from 0)
	int Index;
	
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
					Field.push_back(m);
					Field_Suffix.push_back(suffix);
					num.clear(); suffix.clear();
					prefix_num = true;
					if (ch == 0) break;
				} else {
					suffix.push_back(ch);
				}
			}
		}

		Index = _idx;
	}
};

static bool sortfn(const TAlleleItem *I1, const TAlleleItem *I2)
{
	int smin = min((int)I1->Field.size(), (int)I2->Field.size());
	for (int i=0; i < smin; i++)
	{
		if (I1->Field[i] < I2->Field[i])
		{
			return true;
		} else if (I1->Field[i] > I2->Field[i])
		{
			return false;
		} else {
			if (I1->Field_Suffix[i] < I2->Field_Suffix[i])
			{
				return true;
			} else if (I1->Field_Suffix[i] > I2->Field_Suffix[i])
			{
				return false;
			}		
		}
	}

	return (I1->Field.size() <= I2->Field.size());
}


/**
 *  Sort the HLA alleles
 *
 *  \param hlastr     HLA allele strings
**/
SEXP HIBAG_SortAlleleStr(SEXP hlastr)
{
	CORE_TRY
		const int n = Rf_length(hlastr);

		// HLA alleles
		vector<TAlleleItem> HLA;
		for (int i=0; i < n; i++)
			HLA.push_back(TAlleleItem(CHAR(STRING_ELT(hlastr, i)), i));

		vector<TAlleleItem*> lst;
		for (int i=0; i < n; i++)
			lst.push_back(&HLA[i]);

		// sort
		sort(lst.begin(), lst.end(), sortfn);

		// output
		rv_ans = PROTECT(NEW_CHARACTER(n));
		for (int i=0; i < n; i++)
			SET_STRING_ELT(rv_ans, i, STRING_ELT(hlastr, lst[i]->Index));
		UNPROTECT(1);
	CORE_CATCH
}



// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
//
//  SNP functions
//

/// Detect and correct strand problem

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

SEXP HIBAG_AlleleStrand(SEXP allele1, SEXP afreq1, SEXP I1,
	SEXP allele2, SEXP afreq2, SEXP I2, SEXP if_same_strand, SEXP num)
{
	double *pAF1 = REAL(afreq1);
	double *pAF2 = REAL(afreq2);
	int *pI1 = INTEGER(I1);
	int *pI2 = INTEGER(I2);
	const bool check_strand = (Rf_asLogical(if_same_strand) != TRUE);
	const int n = Rf_asInteger(num);

	CORE_TRY
		// initialize: A-T pair, C-G pair
		map<string, string> MAP;
		MAP["A"] = "T"; MAP["C"] = "G"; MAP["G"] = "C"; MAP["T"] = "A";

		rv_ans = PROTECT(NEW_LIST(4));

		SEXP Flag = PROTECT(NEW_LOGICAL(n));
		SET_ELEMENT(rv_ans, 0, Flag);
		int *out_flag = LOGICAL(Flag);

		int out_n_stand_amb = 0;
		int out_n_mismatch = 0;
		int out_n_swap_stand = 0;

		// loop for each SNP
		for (int i=0; i < n; i++)
		{
			// if true, need flipping A&B alleles
			bool switch_flag = false;

			// if true, need to compare the allele frequencies
			//   0 -- no switch detect
			//   1 -- detect whether switch or not for stand ambiguity
			//   2 -- detect whether switch or not for mismatching alleles
			int switch_freq_detect = 0;

			// ``ref / nonref alleles''
			string s1, s2;
			string p1, p2;
			split_allele(CHAR(STRING_ELT(allele1, pI1[i]-1)), s1, s2);
			split_allele(CHAR(STRING_ELT(allele2, pI2[i]-1)), p1, p2);

			// allele frequency
			const double F1 = pAF1[pI1[i]-1];
			const double F2 = pAF2[pI2[i]-1];

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
							else
								out_n_swap_stand++;
						} else if ( (s1 == MAP[p2]) && (s2 == MAP[p1]) )
						{
							switch_flag = true;
							out_n_swap_stand++;
						} else
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
					out_n_stand_amb ++;
				else
					out_n_mismatch ++;
			}

			out_flag[i] = switch_flag;
		}

		SET_ELEMENT(rv_ans, 1, ScalarInteger(out_n_stand_amb));
		SET_ELEMENT(rv_ans, 2, ScalarInteger(out_n_mismatch));
		SET_ELEMENT(rv_ans, 3, ScalarInteger(out_n_swap_stand));
		UNPROTECT(2);

	CORE_CATCH
}


SEXP HIBAG_AlleleStrand2(SEXP allele1, SEXP allele2)
{
	if (XLENGTH(allele1) != XLENGTH(allele2))
		error("'allele1' and 'allele2' should have the same length.");

	CORE_TRY
		// initialize: A-T pair, C-G pair
		map<string, string> MAP;
		MAP["A"] = "T"; MAP["C"] = "G"; MAP["G"] = "C"; MAP["T"] = "A";

		const int n = XLENGTH(allele1);
		rv_ans = PROTECT(NEW_LOGICAL(n));
		int *pValid = LOGICAL(rv_ans);

		// loop for each SNP
		for (int i=0; i < n; i++)
		{
			// if true, need switch strand
			bool valid = false;

			// ``ref / nonref alleles''
			string s1, s2;
			string p1, p2;
			split_allele(CHAR(STRING_ELT(allele1, i)), s1, s2);
			split_allele(CHAR(STRING_ELT(allele2, i)), p1, p2);

			if (ATGC(s1) && ATGC(s2) && ATGC(p1) && ATGC(p2))
			{
				// check
				if ( (s1 == p1) && (s2 == p2) )
				{
					valid = true;
				} else if ( (s1 == p2) && (s2 == p1) )
				{
					valid = true;
				} else {
					if ( (s1 == MAP[p1]) && (s2 == MAP[p2]) )
					{
						// for example, + C/G <---> - G/C, strand ambi
						valid = true;
					} else if ( (s1 == MAP[p2]) && (s2 == MAP[p1]) )
					{
						valid = true;
					}
				}
			}

			pValid[i] = valid;
		}

		UNPROTECT(1);

	CORE_CATCH
}



// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
//
// HIBAG: Attribute Bagging method
//

/// the number of attribute bagging models allowed
#define MODEL_NUM_LIMIT    256
/// the model container
static CAttrBag_Model* _HIBAG_MODELS_[MODEL_NUM_LIMIT];

/// need a new model
static int _Need_New_HIBAG_Model()
{
	for (int i=0; i < MODEL_NUM_LIMIT; i++)
		if (_HIBAG_MODELS_[i] == NULL) return i;
	throw ErrHLA("No memory space to store a new HIBAG model, "
		"please call \"hlaClose()\" to release unused HIBAG models.");
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


static void model_free(SEXP ptr_obj)
{
	// pointer
	void *ptr = R_ExternalPtrAddr(ptr_obj);
	if (!ptr) return;
	R_ClearExternalPtr(ptr_obj);
	// file ID
	SEXP ID = R_ExternalPtrProtected(ptr_obj);
	int i = Rf_asInteger(ID);
	if (TYPEOF(ID)==INTSXP && Rf_length(ID)>=1)
		INTEGER(ID)[0] = -1;
	// free model
	if ((0 <= i) && (i < MODEL_NUM_LIMIT))
	{
		CAttrBag_Model *m = _HIBAG_MODELS_[i];
		if (m != NULL)
		{
			bool has_error = false;
			try {
				_HIBAG_MODELS_[i] = NULL;
				delete m;
			} catch (exception &E) {
				_LastError = E.what(); has_error = true;
			} catch (const char *E) {
				_LastError = E; has_error = true;
			} catch (...) {
				_LastError = "unknown error!"; has_error = true;
			}
			if (has_error) Rf_error("%s", _LastError.c_str());
		}
	}
}

static SEXP new_model_id(int id, CAttrBag_Model *mod_ptr)
{
	SEXP ans = PROTECT(ScalarInteger(id));
	SEXP ptr = PROTECT(R_MakeExternalPtr(mod_ptr, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(ptr, model_free, (Rboolean)TRUE);
	Rf_setAttrib(ans, install("handle_ptr"), ptr);
	UNPROTECT(2);
	return ans;
}


/**
 *  Build a HIBAG model
 *
 *  \param nSNP       the number of SNPs
 *  \param nSamp      the number of samples
 *  \param nHLA       the number of different HLA alleles
 *  \return a model index
**/
SEXP HIBAG_New(SEXP nSamp, SEXP nSNP, SEXP nHLA)
{
	int n_samp = Rf_asInteger(nSamp);
	if (n_samp <= 0)
		error("Invalid number of samples: %d.", n_samp);
	int n_snp  = Rf_asInteger(nSNP);
	if (n_snp <= 0)
		error("Invalid number of SNPs: %d.", n_snp);
	int n_hla  = Rf_asInteger(nHLA);
	if (n_hla <= 0)
		error("Invalid number of unique HLA alleles: %d.", n_hla);
	CORE_TRY
		int id = _Need_New_HIBAG_Model();
		CAttrBag_Model *obj = _HIBAG_MODELS_[id] = new CAttrBag_Model;
		obj->InitTraining(n_snp, n_samp, n_hla);
		rv_ans = new_model_id(id, obj);
	CORE_CATCH
}


/**
 *  Build a HIBAG model
 *
 *  \param nSNP       the number of SNPs
 *  \param nSamp      the number of samples
 *  \param snp_geno   the SNP genotypes, (0 -- BB, 1 -- AB, 2 -- AA, other -- missing value)
 *  \param nHLA       the number of different HLA alleles
 *  \param H1         the first HLA allele of a HLA type
 *  \param H2         the second HLA allele of a HLA type
 *  \return the model index
**/
SEXP HIBAG_Training(SEXP nSNP, SEXP nSamp, SEXP snp_geno, SEXP nHLA,
	SEXP H1, SEXP H2)
{
	int n_samp = Rf_asInteger(nSamp);
	if (n_samp <= 0)
		error("Invalid number of samples: %d.", n_samp);
	int n_snp  = Rf_asInteger(nSNP);
	if (n_snp <= 0)
		error("Invalid number of SNPs: %d.", n_snp);
	int n_hla  = Rf_asInteger(nHLA);
	if (n_hla <= 0)
		error("Invalid number of unique HLA alleles: %d.", n_hla);
	CORE_TRY
		int id = _Need_New_HIBAG_Model();
		CAttrBag_Model *obj = _HIBAG_MODELS_[id] = new CAttrBag_Model;
		obj->InitTraining(n_snp, n_samp, INTEGER(snp_geno),
			n_hla, INTEGER(H1), INTEGER(H2));
		rv_ans = new_model_id(id, obj);
	CORE_CATCH
}


/**
 *  Close an existing HIBAG model 
 *
 *  \param model      the index in the model list
**/
SEXP HIBAG_Close(SEXP model)
{
	int midx = Rf_asInteger(model);
	if (midx < 0) return R_NilValue;
	CORE_TRY
		_Check_HIBAG_Model(midx);
		INTEGER(model)[0] = -1;
		CAttrBag_Model *m = _HIBAG_MODELS_[midx];
		_HIBAG_MODELS_[midx] = NULL;
		delete m;
	CORE_CATCH
}


/// set and clear GPUExtProcPtr for a try-final block
struct set_gpu_ptr
{
	/// constructor
	set_gpu_ptr(SEXP proc_ptr)
	{
		GPUExtProcPtr = NULL;
		if (!Rf_isNull(proc_ptr))
			GPUExtProcPtr = (TypeGPUExtProc *)R_ExternalPtrAddr(proc_ptr);
	}
	/// destructor
	~set_gpu_ptr()
	{
		GPUExtProcPtr = NULL;
	}
};

/// show the number of threads
inline static void verbose_num_thread(bool verbose)
{
	if (verbose)
	{
	#if RCPP_PARALLEL_USE_TBB
		const int n = tbb::this_task_arena::max_concurrency();
	#else
		const int n = 1;
	#endif
		if (!GPUExtProcPtr)
			Rprintf("# of threads: %d\n", n);
	}
}

/**
 *  Add individual classifiers
 *
 *  \param model           the model index
 *  \param nclassifier     the total number of individual classifiers to be created
 *  \param mtry            the number of variables randomly sampled as candidates for selection
 *  \param prune           if TRUE, perform a parsimonious forward variable selection
 *  \param verbose         show information if TRUE
 *  \param verbose_detail  show more information if TRUE
 *  \param proc_ptr        pointer to functions for an extensible component
**/
SEXP HIBAG_NewClassifiers(SEXP model, SEXP NClassifier, SEXP MTry,
	SEXP Prune, SEXP NThread, SEXP Verbose, SEXP VerboseDetail, SEXP proc_ptr)
{
	const int midx = Rf_asInteger(model);
	const int nclassifier = Rf_asInteger(NClassifier);
	const int mtry = Rf_asInteger(MTry);
	const bool prune = Rf_asLogical(Prune) == TRUE;
	const bool verbose = Rf_asLogical(Verbose) == TRUE;
	const bool verbose_detail = Rf_asLogical(VerboseDetail) == TRUE;

	CORE_TRY
		_Check_HIBAG_Model(midx);
		GetRNGstate();
		set_gpu_ptr set(proc_ptr);

	#if RCPP_PARALLEL_USE_TBB
		const int nthread = Rf_asInteger(NThread);
		tbb::task_arena arena(abs(nthread));
		// tbb::task_scheduler_init init(abs(nthread)); deprecated!
		arena.execute([&]{
	#endif
		verbose_num_thread(verbose && nthread>0);
		if (verbose && nthread>0)
			Rprintf("[-] %s\n", date_text());

		_HIBAG_MODELS_[midx]->BuildClassifiers(nclassifier, mtry,
			prune, verbose, verbose_detail);
	#if RCPP_PARALLEL_USE_TBB
		});
	#endif

		PutRNGstate();
	CORE_CATCH
}


/**
 *  Predict HLA types, output the best-guess and their prob.
 *
 *  \param Model        the model index
 *  \param GenoMat      the pointer to the SNP genotypes
 *  \param NumSamp      the number of samples in GenoMat
 *  \param VoteMethod   the voting method
 *  \param NThread      the number of threads
 *  \param Verbose      whether showing information
 *  \param proc_ptr     pointer to functions for an extensible component
 *  \return H1, H2 and posterior prob.
**/
SEXP HIBAG_Predict_Resp(SEXP Model, SEXP GenoMat, SEXP NumSamp,
	SEXP VoteMethod, SEXP NThread, SEXP Verbose, SEXP proc_ptr)
{
	const int midx = Rf_asInteger(Model);
	const int nSamp = Rf_asInteger(NumSamp);
	const int vote_method = Rf_asInteger(VoteMethod);
	const bool verbose = Rf_asLogical(Verbose)==TRUE;

	CORE_TRY
		_Check_HIBAG_Model(midx);
		CAttrBag_Model &M = *_HIBAG_MODELS_[midx];
		set_gpu_ptr set(proc_ptr);

		rv_ans = PROTECT(NEW_LIST(4));
		SEXP out_H1 = NEW_INTEGER(nSamp);
		SET_ELEMENT(rv_ans, 0, out_H1);
		SEXP out_H2 = NEW_INTEGER(nSamp);
		SET_ELEMENT(rv_ans, 1, out_H2);
		SEXP out_Prob = NEW_NUMERIC(nSamp);
		SET_ELEMENT(rv_ans, 2, out_Prob);
		SEXP out_Matching = NEW_NUMERIC(nSamp);
		SET_ELEMENT(rv_ans, 3, out_Matching);

	#if RCPP_PARALLEL_USE_TBB
		const int nthread = Rf_asInteger(NThread);
		tbb::task_arena arena(nthread);
		// tbb::task_scheduler_init init(nthread); deprecated!
		arena.execute([&]{
	#endif
		verbose_num_thread(verbose);
		M.PredictHLA(INTEGER(GenoMat), nSamp, vote_method,
			INTEGER(out_H1), INTEGER(out_H2), REAL(out_Prob),
			REAL(out_Matching), NULL, NULL, verbose);
	#if RCPP_PARALLEL_USE_TBB
		});
	#endif

		UNPROTECT(1);
	CORE_CATCH
}


/**
 *  Predict HLA types, output the best-guess and the dosages.
 *
 *  \param Model        the model index
 *  \param GenoMat      the pointer to the SNP genotypes
 *  \param NumSamp      the number of samples in GenoMat
 *  \param VoteMethod   the voting method
 *  \param NThread      the number of threads
 *  \param Verbose      whether showing information
 *  \param proc_ptr     pointer to functions for an extensible component
 *  \return H1, H2 and dosages.
**/
SEXP HIBAG_Predict_Dosage(SEXP Model, SEXP GenoMat, SEXP NumSamp,
	SEXP VoteMethod, SEXP NThread, SEXP Verbose, SEXP proc_ptr)
{
	const int midx = Rf_asInteger(Model);
	const int nSamp = Rf_asInteger(NumSamp);
	const int vote_method = Rf_asInteger(VoteMethod);
	const bool verbose = Rf_asLogical(Verbose)==TRUE;

	CORE_TRY
		_Check_HIBAG_Model(midx);
		CAttrBag_Model &M = *_HIBAG_MODELS_[midx];
		set_gpu_ptr set(proc_ptr);

		rv_ans = PROTECT(NEW_LIST(5));
		SEXP out_H1 = NEW_INTEGER(nSamp);
		SET_ELEMENT(rv_ans, 0, out_H1);
		SEXP out_H2 = NEW_INTEGER(nSamp);
		SET_ELEMENT(rv_ans, 1, out_H2);
		SEXP out_Prob = NEW_NUMERIC(nSamp);
		SET_ELEMENT(rv_ans, 2, out_Prob);
		SEXP out_Matching = NEW_NUMERIC(nSamp);
		SET_ELEMENT(rv_ans, 3, out_Matching);
		SEXP out_MatDosage = Rf_allocMatrix(REALSXP, M.nHLA(), nSamp);
		SET_ELEMENT(rv_ans, 4, out_MatDosage);

	#if RCPP_PARALLEL_USE_TBB
		const int nthread = Rf_asInteger(NThread);
		tbb::task_arena arena(nthread);
		// tbb::task_scheduler_init init(nthread); deprecated!
		arena.execute([&]{
	#endif
		verbose_num_thread(verbose);
		M.PredictHLA(INTEGER(GenoMat), nSamp, vote_method,
			INTEGER(out_H1), INTEGER(out_H2), REAL(out_Prob),
			REAL(out_Matching), REAL(out_MatDosage), NULL, verbose);
	#if RCPP_PARALLEL_USE_TBB
		});
	#endif

		UNPROTECT(1);
	CORE_CATCH
}


/**
 *  Predict HLA types, output the best-guess, their prob. and a matrix
 *      of all posterior probabilities
 *
 *  \param model        the model index
 *  \param GenoMat      the pointer to the SNP genotypes
 *  \param nSamp        the number of samples in GenoMat
 *  \param vote_method  the voting method
 *  \param ShowInfo     whether showing information
 *  \param proc_ptr     pointer to functions for an extensible component
 *  \return H1, H2, prob. and a matrix of all probabilities
**/
SEXP HIBAG_Predict_Resp_Prob(SEXP Model, SEXP GenoMat, SEXP NumSamp,
	SEXP VoteMethod, SEXP NThread, SEXP Verbose, SEXP proc_ptr)
{
	const int midx = Rf_asInteger(Model);
	const int nSamp = Rf_asInteger(NumSamp);
	const int vote_method = Rf_asInteger(VoteMethod);
	const bool verbose = Rf_asLogical(Verbose)==TRUE;

	CORE_TRY
		_Check_HIBAG_Model(midx);
		CAttrBag_Model &M = *_HIBAG_MODELS_[midx];
		set_gpu_ptr set(proc_ptr);

		rv_ans = PROTECT(NEW_LIST(6));
		SEXP out_H1 = NEW_INTEGER(nSamp);
		SET_ELEMENT(rv_ans, 0, out_H1);
		SEXP out_H2 = NEW_INTEGER(nSamp);
		SET_ELEMENT(rv_ans, 1, out_H2);
		SEXP out_Prob = NEW_NUMERIC(nSamp);
		SET_ELEMENT(rv_ans, 2, out_Prob);
		SEXP out_Matching = NEW_NUMERIC(nSamp);
		SET_ELEMENT(rv_ans, 3, out_Matching);
		SEXP out_MatDosage = Rf_allocMatrix(REALSXP, M.nHLA(), nSamp);
		SET_ELEMENT(rv_ans, 4, out_MatDosage);
		SEXP out_MatProb = Rf_allocMatrix(REALSXP, M.nHLA()*(M.nHLA()+1)/2, nSamp);
		SET_ELEMENT(rv_ans, 5, out_MatProb);

	#if RCPP_PARALLEL_USE_TBB
		const int nthread = Rf_asInteger(NThread);
		tbb::task_arena arena(nthread);
		// tbb::task_scheduler_init init(nthread); deprecated!
		arena.execute([&]{
	#endif
		verbose_num_thread(verbose);
		M.PredictHLA(INTEGER(GenoMat), nSamp, vote_method,
			INTEGER(out_H1), INTEGER(out_H2), REAL(out_Prob),
			REAL(out_Matching), REAL(out_MatDosage), REAL(out_MatProb),
			verbose);
	#if RCPP_PARALLEL_USE_TBB
		});
	#endif

		UNPROTECT(1);
	CORE_CATCH
}


/**
 *  Create a new individual classifier with specified parameters
 *
 *  \param model        the model index
 *  \param snpidx       the indices of SNP markers
 *  \param samp_num     the bootstrap count
 *  \param freq         the haplotype frequencies
 *  \param hla          the HLA alleles corresponding to the hapltype list
 *  \param haplo        the vector of characters specifying the SNP haplotype list
 *  \param acc          the out-of-bag accuracy
**/
SEXP HIBAG_NewClassifierHaplo(SEXP model, SEXP snpidx,
	SEXP samp_num, SEXP freq, SEXP hla, SEXP haplo, SEXP acc)
{
	int midx = Rf_asInteger(model);
	int nHaplo = Rf_length(freq);
	if (nHaplo != Rf_length(hla))
		error("Invalid length of 'hla'.");
	if (nHaplo != Rf_length(haplo))
		error("Invalid length of 'haplo'.");
	double Acc = Rf_isNull(acc) ? 0.0 : Rf_asReal(acc);

	CORE_TRY
		_Check_HIBAG_Model(midx);

		vector<const char*> HapList(nHaplo);
		for (int i=0; i < nHaplo; i++)
			HapList[i] = CHAR(STRING_ELT(haplo, i));

		CAttrBag_Classifier *I =
			_HIBAG_MODELS_[midx]->NewClassifierAllSamp();
		I->Assign(Rf_length(snpidx), INTEGER(snpidx),
			INTEGER(samp_num), nHaplo, REAL(freq), INTEGER(hla),
			&HapList[0], &Acc);
	CORE_CATCH
}


/**
 *  Clear the individual classifier(s) in the model
 *
 *  \param model        the model index
**/
SEXP HIBAG_ClearClassifier(SEXP model)
{
	int midx = Rf_asInteger(model);
	CORE_TRY
		_Check_HIBAG_Model(midx);
		_HIBAG_MODELS_[midx]->ClearClassifierList();
	CORE_CATCH
}


/**
 *  Get the number of individual component classifiers
 *
 *  \param model        the model index
 *  \return the number of individual classifiers
**/
SEXP HIBAG_GetNumClassifiers(SEXP model)
{
	int midx = Rf_asInteger(model);
	CORE_TRY
		_Check_HIBAG_Model(midx);
		rv_ans = ScalarInteger(_HIBAG_MODELS_[midx]->ClassifierList().size());
	CORE_CATCH
}


/**
 *  Get the details of a specified individual classifier
 *
 *  \param model         the model index
 *  \return a list of individual classifiers
**/
SEXP HIBAG_GetClassifierList(SEXP model, SEXP hla_str)
{
	int midx = Rf_asInteger(model);
	CORE_TRY
		_Check_HIBAG_Model(midx);
		CAttrBag_Model *AB = _HIBAG_MODELS_[midx];

		const size_t nc = AB->ClassifierList().size();
		rv_ans = PROTECT(NEW_LIST(nc));
		for (size_t k=0; k < nc; k++)
		{
			const CAttrBag_Classifier &M = AB->ClassifierList()[k];
			SEXP lst = NEW_LIST(4);
			SET_ELEMENT(rv_ans, k, lst);
			SET_NAMES(lst, hibag_clr_nm);
			{ // samp.num
				const vector<int> &ns = M.BootstrapCount();
				SEXP nl = NEW_INTEGER(ns.size());
				SET_ELEMENT(lst, 0, nl);
				memcpy(INTEGER(nl), &ns[0], sizeof(int)*ns.size());
			}
			{ // haplos
				SEXP dt = NEW_LIST(3);
				SET_ELEMENT(lst, 1, dt);
				const size_t nHaplo = M.nHaplo();
				const CHaplotypeList &Haplo = M.Haplotype();

				// haplos$freq
				SEXP freq = NEW_NUMERIC(nHaplo);
				SET_ELEMENT(dt, 0, freq);
				double *pF = REAL(freq);
				for (size_t i=0; i < nHaplo; i++)
					pF[i] = Haplo.List[i].Freq;

				// haplos$hla
				SEXP hla = NEW_CHARACTER(nHaplo);
				SET_ELEMENT(dt, 1, hla);
				size_t hla_i = 0;
				for (size_t i=0; i < Haplo.LenPerHLA.size(); i++)
				{
					SEXP s = STRING_ELT(hla_str, i);
					for (size_t j=Haplo.LenPerHLA[i]; j > 0; j--)
						SET_STRING_ELT(hla, hla_i++, s);
				}

				// haplos$haplo
				SEXP haplo = NEW_CHARACTER(nHaplo);
				SET_ELEMENT(dt, 2, haplo);
				const size_t nsnp = M.nSNP();
				for (size_t i=0; i < nHaplo; i++)
				{
					SET_STRING_ELT(haplo, i,
						mkChar(Haplo.List[i].HaploToStr(nsnp).c_str()));
				}

				// convert to data.frame
				SET_NAMES(dt, hibag_clr_haplo_nm);
				SET_CLASS(dt, hibag_data_frame);
				SEXP rs = NEW_INTEGER(nHaplo);
				Rf_setAttrib(dt, R_RowNamesSymbol, rs);
				int *pR = INTEGER(rs);
				for (size_t i=0; i < nHaplo; i++) pR[i] = i+1;
			}
			{ // snpidx
				const vector<int> &I = M.SNPIndex();
				SEXP idx = NEW_INTEGER(I.size());
				SET_ELEMENT(lst, 2, idx);
				int *p = INTEGER(idx);
				for (size_t i=0; i < I.size(); i++) p[i] = I[i] + 1;
			}
			{ // outofbag.acc
				SET_ELEMENT(lst, 3, ScalarReal(M.OutOfBag_Accuracy()));
			}
		}
		UNPROTECT(1);

	CORE_CATCH
}


/**
 *  Get the number of SNPs, the number of haplotypes and out-of-bag accuracy
 *
 *  \param model         the model index
 *  \return a real vector
**/
SEXP HIBAG_GetLastClassifierInfo(SEXP model)
{
	int midx = Rf_asInteger(model);
	CORE_TRY
		_Check_HIBAG_Model(midx);
		CAttrBag_Model *AB = _HIBAG_MODELS_[midx];
		rv_ans = PROTECT(NEW_NUMERIC(3));
		double *p = REAL(rv_ans);
		const size_t nc = AB->ClassifierList().size();
		if (nc > 0)
		{
			const CAttrBag_Classifier &M = AB->ClassifierList()[nc-1];
			p[0] = M.nSNP();
			p[1] = M.nHaplo();
			p[2] = M.OutOfBag_Accuracy();
		} else {
			p[0] = p[1] = p[2] = R_NaN;
		}
		UNPROTECT(1);
	CORE_CATCH
}


/**
 *  Estimate the confusion matrix
 *
 *  \param n_hla         the number of different HLA alleles
 *  \param init_mat      the initial confusion matrix without any ambiguous state
 *  \param n_DConfusion  the number of double confusions
 *  \param D_mat
 *  \return a confusion matrix
**/
SEXP HIBAG_Confusion(SEXP n_hla, SEXP init_mat, SEXP n_DConfusion,
	SEXP D_mat)
{
	// the max number of iterations
	const int N_MAX_ITERATION = 100;
	// the number of unique HLA alleles
	const int nHLA = Rf_asInteger(n_hla);
	// the number of double confusions
	const int nDConf = Rf_asInteger(n_DConfusion);

	CORE_TRY
		const size_t SIZE_MAT = sizeof(double)*nHLA*(nHLA+1);

		#define INDEX(T, P, var) var[(nHLA+1)*T + P]

		rv_ans = Rf_allocMatrix(REALSXP, nHLA+1, nHLA);
		double *out_mat = REAL(rv_ans);

		vector<double> TmpMat(nHLA*(nHLA+1));
		double *tmp_mat = &TmpMat[0];

		// initial values
		memcpy(out_mat, REAL(init_mat), SIZE_MAT);
		for (int i=0; i < nDConf; i++)
		{
			int *T = INTEGER(D_mat) + i*4;
			int *P = INTEGER(D_mat) + i*4 + 2;
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
			memcpy(out_mat, REAL(init_mat), SIZE_MAT);
			for (int i=0; i < nDConf; i++)
			{
				int *T = INTEGER(D_mat) + i*4;
				int *P = INTEGER(D_mat) + i*4 + 2;

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

	CORE_CATCH
}


/**
 *  Detect the storage mode of a PLINK BED file
 *
 *  \param bedfn         the file name of PLINK BED file
**/
SEXP HIBAG_BEDFlag(SEXP bedfn)
{
	const char *fn = CHAR(STRING_ELT(bedfn, 0));
	CORE_TRY
		ifstream file(fn, ios::binary);
		if (!file.good())
			throw ErrHLA("Cannot open the file %s.", fn);
		char prefix[3];
		file.read(prefix, 3);
		if ((prefix[0] != 0x6C) || (prefix[1] != 0x1B))
			throw ErrHLA("Invalid prefix in the PLINK BED file.");
		rv_ans = ScalarInteger((unsigned char)prefix[2]);
	CORE_CATCH
}


/**
 *  Convert from a PLINK BED file
 *
 *  \param bedfn         the file name of PLINK BED file
 *  \param n_samp        the number of samples
 *  \param n_snp         the number of SNPs
 *  \param n_save_snp    the number of imported SNPs
 *  \param snp_flag      the flag of whether it is selected
 *  \return a matrix of genotypes
**/
SEXP HIBAG_ConvBED(SEXP bedfn, SEXP n_samp, SEXP n_snp, SEXP n_save_snp,
	SEXP snp_flag)
{
	const char *fn     = CHAR(STRING_ELT(bedfn, 0));
	const int NumSamp  = Rf_asInteger(n_samp);
	const int NumSNP   = Rf_asInteger(n_snp);
	const int NumSvSNP = Rf_asInteger(n_save_snp);
	const int *pflag   = LOGICAL(snp_flag);

	CORE_TRY
		// open file
		ifstream file(fn, ios::binary);
		if (!file.good())
			throw ErrHLA("Fail to open the file \"%s\".", fn);

		// read prefix
		char prefix[3];
		file.read(prefix, 3);
		if ((prefix[0] != 0x6C) || (prefix[1] != 0x1B))
			throw ErrHLA("Invalid prefix in the PLINK BED file.");
		int mode = prefix[2];

		// determine the values of packed genotypes
		int nRe, nPack, nNumPack, nNum;
		if (mode == 0)
		{
			// the individual-major mode
			nRe = NumSNP % 4;
			nNumPack = NumSNP / 4;
			nPack = (nRe > 0) ? (nNumPack + 1) : nNumPack;
			nNum = NumSamp;
		} else {
			// the SNP-major mode
			nRe = NumSamp % 4;
			nNumPack = NumSamp / 4;
			nPack = (nRe > 0) ? (nNumPack + 1) : nNumPack;
			nNum = NumSNP;
		}

		vector<char> srcgeno(nPack);
		vector<int> dstgeno((nNumPack+1) * 4);
		static const int cvt[4] = { 2, NA_INTEGER, 1, 0 };
		int I_SNP = 0;

		// output
		rv_ans = Rf_allocMatrix(INTSXP, NumSvSNP, NumSamp);

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
			if (mode == 0)
			{
				// the individual-major mode
				int *pI = INTEGER(rv_ans) + i * NumSvSNP;
				for (int j=0; j < NumSNP; j++)
				{
					if (pflag[j])
						*pI++ = dstgeno[j];
				}
			} else {
				// the SNP-major mode
				if (pflag[i])
				{
					int *pI = INTEGER(rv_ans) + I_SNP;
					I_SNP ++;
					for (int j=0; j < NumSamp; j++)
					{
						*pI = dstgeno[j];
						pI += NumSvSNP;
					}
				}
			}
		}

	CORE_CATCH
}


/**
 *  Merge multiple sequences with asterisk
**/
SEXP HIBAG_SeqMerge(SEXP seq)
{
	if (!Rf_isNull(seq))
	{
		const int len = XLENGTH(seq);
		if (len <= 0)
			error("Internal error in 'HIBAG_SeqMerge()'.");
		// get the maximum length
		int nmax = -1;
		for (int i=0; i < len; i++)
		{
			int m = strlen(CHAR(STRING_ELT(seq, i)));
			if (m > nmax) nmax = m;
		}
		// the first sequence
		string ss(nmax, '-');
		const char *p = CHAR(STRING_ELT(seq, 0));
		int j = 0;
		for (; (j < nmax) && (*p); j++) ss[j] = *p++;
		for (; j < nmax; j++) ss[j] = '*';
		// other sequences
		for (int i=1; i < len; i++)
		{
			p = CHAR(STRING_ELT(seq, i));
			j = 0;
			for (; (j < nmax) && (*p); j++, p++)
			{
				if (*p != ss[j]) ss[j] = '*';
			}
			for (; j < nmax; j++) ss[j] = '*';
		}
		// output
		return mkString(ss.c_str());
	} else
		return Rf_ScalarString(NA_STRING);
}


/**
 *  Remove dots in the sequences
**/
SEXP HIBAG_SeqRmDot(SEXP ref, SEXP seq)
{
	const char *s = CHAR(STRING_ELT(ref, 0));
	bool has_dot = false;
	for (; *s; s++)
	{
		if (*s == '.')
			{ has_dot = true; break; }
	}

	if (has_dot)
	{
		PROTECT(ref); PROTECT(seq);
		// remove dots of reference
		string ref_str;
		for (s = CHAR(STRING_ELT(ref, 0)); *s; s++)
		{
			if (*s != '.')
				ref_str.push_back(*s);
		}
		// sequences
		int n = XLENGTH(seq) / 2;
		for (int i=n; i < 2*n; i++)
		{
			string ss;
			s = CHAR(STRING_ELT(ref, 0));
			const char *p = CHAR(STRING_ELT(seq, i));
			for (; *p && *s; p++, s++)
			{
				if (*s != '.')
					ss.push_back(*p);
			}
			SET_STRING_ELT(seq, i, mkChar(ss.c_str()));
		}
		// set
		SET_STRING_ELT(ref, 0, mkChar(ref_str.c_str()));
		UNPROTECT(2);
	}

	return R_NilValue;
}


/**
 *  Calculate the distances among different HLA alleles
**/
SEXP HIBAG_Distance(SEXP NumHLA, SEXP Idx, SEXP Freq, SEXP HaploStr)
{
	int num_hla = Rf_asInteger(NumHLA);
	int n = LENGTH(Idx);
	int *I = INTEGER(Idx);
	double *freq = REAL(Freq);

	SEXP freq_sum = Rf_allocMatrix(REALSXP, num_hla, num_hla);
	PROTECT(freq_sum);
	double *p_freq_sum = REAL(freq_sum);
	memset(p_freq_sum, 0, sizeof(double)*num_hla*num_hla);

	SEXP dist_sum = Rf_allocMatrix(REALSXP, num_hla, num_hla);
	PROTECT(dist_sum);
	double *p_dist_sum = REAL(dist_sum);
	memset(p_dist_sum, 0, sizeof(double)*num_hla*num_hla);

	for (int i=0; i < n; i++)
	{
		for (int j=i; j < n; j++)
		{
			if (I[i]!=NA_INTEGER && I[j]!=NA_INTEGER)
			{
				const char *s1 = CHAR(STRING_ELT(HaploStr, i));
				const char *s2 = CHAR(STRING_ELT(HaploStr, j));
				int d = 0;
				for (; *s1 && *s2; s1++, s2++)
					if (*s1 != *s2) d++;
				double f = freq[i] * freq[j];
				int ii = (I[i]-1) * num_hla + (I[j]-1);
				p_freq_sum[ii] += f;
				p_dist_sum[ii] += f*d;
			}
		}
	}

	for (int i=0; i < num_hla; i++)
	{
		for (int j=i; j < num_hla; j++)
		{
			int ii = i * num_hla + j;
			int jj = j * num_hla + i;
			p_dist_sum[jj] = p_dist_sum[ii] = p_dist_sum[ii] / p_freq_sum[ii];
		}
	}

	UNPROTECT(2);
	return dist_sum;
}


// whether compile the algorithm with specified targets or not
extern const bool HIBAG_ALGORITHM_SSE2;
extern const bool HIBAG_ALGORITHM_SSE2_POPCNT;
extern const bool HIBAG_ALGORITHM_SSE4_2;
extern const bool HIBAG_ALGORITHM_AVX;
extern const bool HIBAG_ALGORITHM_AVX2;
extern const bool HIBAG_ALGORITHM_AVX512F;
extern const bool HIBAG_ALGORITHM_AVX512BW;
extern const bool HIBAG_ALGORITHM_AVX512VPOPCNTDQ;

static void target_add(string &s, const char *code, bool work)
{
	if (work)
		s.append(" ").append(code);
}


/**
 *  Get the version and SSE information
**/
SEXP HIBAG_Kernel_Version()
{
	SEXP ans = PROTECT(NEW_LIST(3));
	// version
	SEXP I = NEW_INTEGER(2);
	SET_ELEMENT(ans, 0, I);
	INTEGER(I)[0] = HIBAG_KERNEL_VERSION >> 8;
	INTEGER(I)[1] = HIBAG_KERNEL_VERSION & 0xFF;

	// CPU information
	SEXP info = NEW_CHARACTER(3);
	SET_ELEMENT(ans, 1, info);
	SET_STRING_ELT(info, 0, mkChar(CPU_Info()));
	// compiler information
#ifdef __VERSION__
	string version = __VERSION__;
#else
	string version;
#endif
#ifdef __GNUC__
	char buf_compiler[128] = { 0 };
	#ifndef __GNUC_PATCHLEVEL__
	#   define __GNUC_PATCHLEVEL__    0
	#endif
	snprintf(buf_compiler, sizeof(buf_compiler), "GNUG_v%d.%d.%d",
		__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
	string gnug = buf_compiler;
#else
	string gnug;
#endif
	string s;
	if (!version.empty())
	{
		if (!gnug.empty())
			s = version + ", " + gnug;
		else
			s = version;
	} else {
		if (!gnug.empty())
			s = gnug;
		else
			s = "Unknown compiler";
	}
	SET_STRING_ELT(info, 1, mkChar(s.c_str()));
	// supported implementation
	s = "Algorithm SIMD:";
	target_add(s, "SSE2", HIBAG_ALGORITHM_SSE2);
	target_add(s, "SSE4.2", HIBAG_ALGORITHM_SSE4_2);
	target_add(s, "AVX", HIBAG_ALGORITHM_AVX);
	target_add(s, "AVX2", HIBAG_ALGORITHM_AVX2);
	target_add(s, "AVX512F", HIBAG_ALGORITHM_AVX512F);
	target_add(s, "AVX512BW", HIBAG_ALGORITHM_AVX512BW);
	target_add(s, "AVX512VPOPCNTDQ", HIBAG_ALGORITHM_AVX512VPOPCNTDQ);
	SET_STRING_ELT(info, 2, mkChar(s.c_str()));

	// using Intel TBB or not
#if RCPP_PARALLEL_USE_TBB
#if defined(__APPLE__) && (defined(__arm__) || defined(__arm64__))
	int ntot = 1;
#else
	int ntot = tbb::this_task_arena::max_concurrency();
#endif
	SET_ELEMENT(ans, 2, ScalarInteger(ntot));
#else
	SET_ELEMENT(ans, 2, ScalarInteger(NA_INTEGER));
#endif
	// output
	UNPROTECT(1);
	return ans;
}


/**
 *  Get the version and SSE information
**/
SEXP HIBAG_Kernel_SetTarget(SEXP CPU)
{
	const char *cpu = CHAR(STRING_ELT(CPU, 0));
	CAlg_Prediction::Init_Target_IFunc(cpu);
	return R_NilValue;
}


/**
 *  Get the version and SSE information
**/
SEXP HIBAG_Clear_GPU()
{
	GPUExtProcPtr = NULL;
	return R_NilValue;
}


/**
 *  Called by .onLoad()
**/
SEXP HIBAG_Init(SEXP data_lst)
{
	hibag_data_frame = VECTOR_ELT(data_lst, 0);
	hibag_clr_nm = VECTOR_ELT(data_lst, 1);
	hibag_clr_haplo_nm = VECTOR_ELT(data_lst, 2);
	return R_NilValue;
}


/**
 *  Sum the vectors in the list with weights
**/
SEXP HIBAG_SumList(SEXP weight, SEXP list)
{
	const int n  = Rf_length(weight);
	const int nl = Rf_length(list);
	if (n != nl)
		Rf_error("HIBAG_SumList error.");
	if (n > 0)
	{
		const int nx = Rf_length(VECTOR_ELT(list, 0));
		SEXP rv_ans = PROTECT(NEW_NUMERIC(nx));
		double *pv = REAL(rv_ans);
		memset(pv, 0, sizeof(double)*nx);  // zero fill
		const double *pw = REAL(weight);
		for (int i=0; i < n; i++)
		{
			SEXP y = VECTOR_ELT(list, i);
			if (nx != Rf_length(y))
				Rf_error("HIBAG_SumList, list length error.");
			const double *py = REAL(y);
			const double w = pw[i];
			for (int j=0; j < nx; j++) pv[j] += w * py[j];
		}
		UNPROTECT(1);
		return rv_ans;
	} else {
		return Rf_ScalarLogical(FALSE);
	}
}


/**
 *  Update out_prob by adding weighted in_p
**/
SEXP HIBAG_UpdateAddProbW(SEXP out_prob, SEXP ii, SEXP in_p, SEXP weight,
	SEXP matching)
{
	// check
	if (!Rf_isMatrix(out_prob))
		Rf_error("HIBAG_UpdateAddProbW out_prob error.");
	if (!Rf_isMatrix(in_p))
		Rf_error("HIBAG_UpdateAddProbW in_p error.");
	const int *n1 = INTEGER(GET_DIM(out_prob));
	const int *n2 = INTEGER(GET_DIM(in_p));
	if (n1[1] != n2[1])
		Rf_error("HIBAG_UpdateAddProbW dim(prob) error.");
	if (Rf_length(ii) != n2[0])
		Rf_error("HIBAG_UpdateAddProbW ii error.");
	if (Rf_length(matching) != n2[1])
		Rf_error("HIBAG_UpdateAddProbW matching error.");
	// update
	const double w = Rf_asReal(weight);
	const double *pm = !Rf_isNull(matching) ? REAL(matching) : NULL;
	double *po = REAL(out_prob);
	double *pi = REAL(in_p);
	const int *I = INTEGER(ii);
	for (int i=0; i < n2[1]; i++)
	{
		const double w2 = pm ? w * pm[i] : w;
		for (int j=0; j < n2[0]; j++)
			po[I[j] - 1] += pi[j] * w2;
		po += n1[0];
		pi += n2[0];
	}
	return out_prob;
}


/**
 *  Normalize the prob matrix for each column
**/
SEXP HIBAG_NormalizeProb(SEXP prob)
{
	if (!Rf_isMatrix(prob))
		Rf_error("HIBAG_NormalizeProb prob error.");
	const int *dm = INTEGER(GET_DIM(prob));
	double *p = REAL(prob);
	for (int i=0; i < dm[1]; i++)
	{
		double sum = 0;
		for (int j=0; j < dm[0]; j++) sum += p[j];
		for (int j=0; j < dm[0]; j++) p[j] /= sum;
		p += dm[0];
	}
	return prob;
}


// -----------------------------------------------------------------------
// -----------------------------------------------------------------------

/// Initialize the package
void R_init_HIBAG(DllInfo *info)
{
	#define CALL(name, num)    { #name, (DL_FUNC)&name, num }

	static R_CallMethodDef callMethods[] =
	{
		CALL(HIBAG_Init, 1),
		CALL(HIBAG_AlleleStrand, 8),
		CALL(HIBAG_AlleleStrand2, 2),
		CALL(HIBAG_BEDFlag, 1),
		CALL(HIBAG_GetNumClassifiers, 1),
		CALL(HIBAG_GetClassifierList, 2),
		CALL(HIBAG_Close, 1),
		CALL(HIBAG_Confusion, 4),
		CALL(HIBAG_ConvBED, 5),
		CALL(HIBAG_Distance, 4),
		CALL(HIBAG_Kernel_Version, 0),
		CALL(HIBAG_New, 3),
		CALL(HIBAG_NewClassifierHaplo, 7),
		CALL(HIBAG_NewClassifiers, 8),
		CALL(HIBAG_Predict_Resp, 7),
		CALL(HIBAG_Predict_Dosage, 7),
		CALL(HIBAG_Predict_Resp_Prob, 7),
		CALL(HIBAG_Training, 6),
		CALL(HIBAG_SortAlleleStr, 1),
		CALL(HIBAG_SeqMerge, 1),
		CALL(HIBAG_SeqRmDot, 2),
		CALL(HIBAG_Clear_GPU, 0),
		{ NULL, NULL, 0 }
	};

	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	memset((void*)_HIBAG_MODELS_, 0, sizeof(_HIBAG_MODELS_));
}

/// Finalize the package
void R_unload_HIBAG(DllInfo *info)
{
	try
	{
		for (int i=0; i < MODEL_NUM_LIMIT; i++)
		{
			CAttrBag_Model *m = _HIBAG_MODELS_[i];
			try
			{
				if (m != NULL)
				{
					_HIBAG_MODELS_[i] = NULL;
					delete m;
				}
			} catch(...) {}
		}
	} catch(...) {}
}

} // extern "C"
