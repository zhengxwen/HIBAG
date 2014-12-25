\name{HIBAG-package}
\alias{HIBAG-package}
\alias{HIBAG}
\docType{package}
\title{HLA Genotype Imputation with Attribute Bagging}
\description{
	To impute HLA types from unphased SNP data using an attribute bagging
method.
}
\details{
	\tabular{ll}{
		Package: \tab HIBAG\cr
		Type:    \tab Package\cr
		Version: \tab 1.2.5\cr
		License: \tab GPL version 3\cr
	}

	HIBAG is a state of the art software package for imputing HLA types using
SNP data, and it uses the R statistical programming language. HIBAG is highly
accurate, computationally tractable, and can be used by researchers with
published parameter estimates instead of requiring access to large training
sample datasets. It combines the concepts of attribute bagging, an ensemble
classifier method, with haplotype inference for SNPs and HLA types. Attribute
bagging is a technique which improves the accuracy and stability of classifier
ensembles using bootstrap aggregating and random variable selection.

\bold{Features:} \cr
1) HIBAG can be used by researchers with published parameter estimates
	(\url{http://www.biostat.washington.edu/~bsweir/HIBAG/}) instead of
	requiring access to large training sample datasets. \cr
2) A typical HIBAG parameter file contains only haplotype frequencies at
	different SNP subsets rather than individual training genotypes. \cr
3) SNPs within the xMHC region (chromosome 6) are used for imputation. \cr
4) HIBAG employs unphased genotypes of unrelated individuals as a training
	set. \cr
5) HIBAG supports parallel computing with R. \cr
}
\references{
	Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS;
	HIBAG -- HLA Genotype Imputation with Attribute Bagging;
	(Abstract 294, Platform/Oral Talk);
	Present at the 62nd Annual Meeting of the American Society of Human Genetics,
	November 9, 2012 in San Francisco, California.

	Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS;
	HIBAG -- HLA Genotype Imputation with Attribute Bagging.
	The Pharmacogenomics Journal. doi: 10.1038/tpj.2013.18.
	\url{http://dx.doi.org/10.1038/tpj.2013.18}
}
\author{
	Xiuwen Zheng [aut, cre] \email{zhengx@u.washington.edu},
	Bruce S. Weir [ctb, ths] \email{bsweir@u.washington.edu}
}

\examples{
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

head(HLA_Type_Table)
dim(HLA_Type_Table)  # 60 13

summary(HapMap_CEU_Geno)


######################################################################

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
	H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
	H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
	locus=hla.id, assembly="hg19")

# divide HLA types randomly
set.seed(100)
hlatab <- hlaSplitAllele(hla, train.prop=0.5)
names(hlatab)
# "training"   "validation"
summary(hlatab$training)
summary(hlatab$validation)

# SNP predictors within the flanking region on each side
region <- 500   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
	hla.id, region*1000, assembly="hg19")
length(snpid)  # 275

# training and validation genotypes
train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
	snp.sel=match(snpid, HapMap_CEU_Geno$snp.id),
	samp.sel=match(hlatab$training$value$sample.id, HapMap_CEU_Geno$sample.id))
test.geno <- hlaGenoSubset(HapMap_CEU_Geno,
	samp.sel=match(hlatab$validation$value$sample.id, HapMap_CEU_Geno$sample.id))

# train a HIBAG model
set.seed(100)
# please use "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=4,
	verbose.detail=TRUE)
summary(model)

# validation
pred <- predict(model, test.geno)
summary(pred)

# compare
(comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model,
	call.threshold=0))
(comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model,
	call.threshold=0.5))


# save the parameter file
mobj <- hlaModelToObj(model)
save(mobj, file="HIBAG_model.RData")
save(test.geno, file="testgeno.RData")
save(hlatab, file="HLASplit.RData")

# Clear Workspace
hlaClose(model)  # release all resources of model
rm(list = ls())


######################################################################

# NOW, load a HIBAG model from the parameter file
mobj <- get(load("HIBAG_model.RData"))
model <- hlaModelFromObj(mobj)

# validation
test.geno <- get(load("testgeno.RData"))
hlatab <- get(load("HLASplit.RData"))

pred <- predict(model, test.geno)
# compare
(comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model,
	call.threshold=0.5))


#########################################################################
# import a PLINK BED file
#
bed.fn <- system.file("extdata", "HapMap_CEU.bed", package="HIBAG")
fam.fn <- system.file("extdata", "HapMap_CEU.fam", package="HIBAG")
bim.fn <- system.file("extdata", "HapMap_CEU.bim", package="HIBAG")
hapmap.ceu <- hlaBED2Geno(bed.fn, fam.fn, bim.fn, assembly="hg19")


#########################################################################
# predict
#
pred <- predict(model, hapmap.ceu, type="response")
head(pred$value)
#    sample.id allele1 allele2      prob
# 1    NA10859   01:01   03:01 0.9374860
# 2    NA11882   01:01   29:02 0.9999835
# ...


# delete the temporary files
unlink("HIBAG_model.RData", force=TRUE)
unlink("testgeno.RData", force=TRUE)
unlink("HLASplit.RData", force=TRUE)
}

\keyword{HLA}
\keyword{SNP}
\keyword{genetics}
