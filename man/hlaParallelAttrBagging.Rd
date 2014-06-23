\name{hlaParallelAttrBagging}
\alias{hlaParallelAttrBagging}
\title{
	Build a HIBAG model via parallel computation
}
\description{
	To build a HIBAG model for predicting HLA types via parallel computation.
}
\usage{
hlaParallelAttrBagging(cl, hla, snp, auto.save="",
	nclassifier=100, mtry=c("sqrt", "all", "one"), prune=TRUE, rm.na=TRUE,
	stop.cluster=FALSE, verbose=TRUE)
}
\arguments{
	\item{cl}{a cluster object, created by the package \link[parallel]{parallel}
		or \href{http://CRAN.R-project.org/package=snow}{snow}; if \code{NULL}
		is given, a uniprocessor implementation will be performed}
	\item{hla}{training HLA types, an object of \code{\link{hlaAlleleClass}}}
	\item{snp}{training SNP genotypes, an object of
		\code{\link{hlaSNPGenoClass}}}
	\item{auto.save}{specify a autosaved file, see details}
	\item{nclassifier}{the total number of individual classifiers}
	\item{mtry}{a character or a numeric value, the number of variables randomly
		sampled as candidates for each selection. See details}
	\item{prune}{if TRUE, to perform a parsimonious forward variable selection,
		otherwise, exhaustive forward variable selection. See details}
	\item{rm.na}{if TRUE, remove the samples with missing HLA types}
	\item{stop.cluster}{\code{TRUE}: stop cluster nodes after computing}
	\item{verbose}{if TRUE, show information}
}
\details{
	\code{mtry} (the number of variables randomly sampled as candidates for
		each selection):
	\code{"sqrt"}, using the square root of the total number of candidate SNPs;
	\code{"all"}, using all candidate SNPs;
	\code{"one"}, using one SNP;
	\code{an integer}, specifying the number of candidate SNPs;
	\code{0 < r < 1}, the number of candidate SNPs is
		"r * the total number of SNPs".

	\code{prune}: there is no significant difference on accuracy between
parsimonious and exhaustive forward variable selections. If \code{prune = TRUE},
the searching algorithm performs a parsimonious forward variable selection:
if a new SNP predictor reduces the current out-of-bag accuracy, then it is
removed from the candidate SNP set for future searching. Parsimonious selection
helps to improve the computational efficiency by reducing the searching times
of non-informative SNP markers.

	If \code{auto.save=""}, the function returns a HIBAG model (an object of
\code{\link{hlaAttrBagClass}}); otherwise, there is no return.
}
\value{
	Return an object of \code{\link{hlaAttrBagClass}} if \code{auto.save}
is specified.
}
\references{
	Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS;
	HIBAG -- HLA Genotype Imputation with Attribute Bagging;
	(Abstract 294, Platform/Oral Talk);
	Present at the 62nd Annual Meeting of the American Society of Human Genetics,
	November 9, 2012 in San Francisco, California.

	Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS;
	HIBAG -- HLA Genotype Imputation with Attribute Bagging.
	Pharmacogenomics Journal. doi: 10.1038/tpj.2013.18.
	\url{http://dx.doi.org/10.1038/tpj.2013.18}
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{hlaClose}}
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
	snp.sel = match(snpid, HapMap_CEU_Geno$snp.id),
	samp.sel = match(hlatab$training$value$sample.id, HapMap_CEU_Geno$sample.id))
test.geno <- hlaGenoSubset(HapMap_CEU_Geno,
	samp.sel=match(hlatab$validation$value$sample.id, HapMap_CEU_Geno$sample.id))


#############################################################################

library(parallel)

# use option cl.core to choose an appropriate cluster size.
cl <- makeCluster(getOption("cl.cores", 2))
set.seed(100)

# train a HIBAG model in parallel
# please use "nclassifier=100" when you use HIBAG for real data
hlaParallelAttrBagging(cl, hlatab$training, train.geno, nclassifier=4,
	auto.save="tmp_model.RData", stop.cluster=TRUE)

mobj <- get(load("tmp_model.RData"))
summary(mobj)
model <- hlaModelFromObj(mobj)

# validation
pred <- predict(model, test.geno)
summary(pred)

# compare
hlaCompareAllele(hlatab$validation, pred, allele.limit=model)$overall


# since 'stop.cluster=TRUE' used in 'hlaParallelAttrBagging'
# need a new cluster
cl <- makeCluster(getOption("cl.cores", 2))

pred <- predict(model, test.geno, cl=cl)
summary(pred)

# stop parallel nodes
stopCluster(cl)
}

\keyword{HLA}
\keyword{genetics}
