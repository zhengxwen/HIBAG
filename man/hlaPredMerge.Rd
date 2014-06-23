\name{hlaPredMerge}
\alias{hlaPredMerge}
\title{
	Merge prediction results from multiple HIBAG models
}
\description{
	Return an object of \code{\link{hlaAlleleClass}}, which contains predicted
HLA types.
}
\usage{
hlaPredMerge(..., weight=NULL, equivalence=NULL)
}
\arguments{
	\item{...}{The object(s) of \code{\link{hlaAlleleClass}}, having a field
		of 'postprob', and returned by
		\code{predict(..., type="response+prob", vote="majority")}}
	\item{weight}{the weight used for each prediction; if \code{NULL},
		equal weights}
	\item{equivalence}{a \code{data.frame} with two columns, the first column
		for new equivalent alleles, and the second for the alleles possibly
		existed in the object(s) passed to this function}
}
\details{
	Calculate a new probability matrix for each pair of HLA alleles, by
averaging (posterior) probabilities from all models with specified weights.
If \code{equivalence} is specified, multiple alleles might be collapsed into
one class.
}
\value{
	Return a \code{\link{hlaAlleleClass}} object.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{hlaAllele}},
	\code{\link{predict.hlaAttrBagClass}}
}

\examples{
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

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

# train HIBAG models
set.seed(100)

# please use "nclassifier=100" when you use HIBAG for real data
m1 <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=2,
	verbose.detail=TRUE)
m2 <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=2,
	verbose.detail=TRUE)


# validation
pd1 <- predict(m1, test.geno, type="response+prob", vote="majority")
pd2 <- predict(m2, test.geno, type="response+prob", vote="majority")

hlaCompareAllele(hlatab$validation, pd1)$overall
# acc.haplo = 0.8269231

hlaCompareAllele(hlatab$validation, pd2)$overall
# acc.haplo = 0.8076923


# merge predictions from multiple models, by voting from all classifiers
pd <- hlaPredMerge(pd1, pd2, weight=c(1,1))

hlaCompareAllele(hlatab$validation, pd)$overall
# acc.haplo = 0.9230769
}

\keyword{HLA}
\keyword{genetics}
