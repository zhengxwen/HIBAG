\name{hlaOutOfBag}
\alias{hlaOutOfBag}
\title{
	Out-of-bag estimation of overall accuracy, per-allele sensitivity, etc
}
\description{
	Out-of-bag estimation of overall accuracy, per-allele sensitivity,
specificity, positive predictive value, negative predictive value and
call rate.
}
\usage{
hlaOutOfBag(model, hla, snp, call.threshold=NaN, verbose=TRUE)
}
\arguments{
	\item{model}{an object of \code{\link{hlaAttrBagClass}} or
		\code{\link{hlaAttrBagObj}}}
	\item{hla}{the training HLA types, an object of
		\code{\link{hlaAlleleClass}}}
	\item{snp}{the training SNP genotypes, an object of
		\code{\link{hlaSNPGenoClass}}}
	\item{call.threshold}{the specified call threshold; if \code{NaN}, no
		threshold is used}
	\item{verbose}{if TRUE, show information}
}
\value{
	Return \code{\link{hlaAlleleClass}}.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaCompareAllele}}, \code{\link{hlaReport}}
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

# SNP predictors within the flanking region on each side
region <- 500   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
	hla.id, region*1000, assembly="hg19")
length(snpid)  # 275

# training and validation genotypes
geno <- hlaGenoSubset(HapMap_CEU_Geno,
	snp.sel = match(snpid, HapMap_CEU_Geno$snp.id),
	samp.sel = match(hla$value$sample.id, HapMap_CEU_Geno$sample.id))

# train a HIBAG model
set.seed(100)
# please use "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hla, geno, nclassifier=4)
summary(model)

# out-of-bag estimation
(comp <- hlaOutOfBag(model, hla, geno, call.threshold=NaN, verbose=TRUE))

# report
hlaReport(comp, type="txt")

hlaReport(comp, type="tex")

hlaReport(comp, type="html")
}

\keyword{HLA}
\keyword{genetics}
