\name{plot.hlaAttrBagObj}
\alias{plot.hlaAttrBagObj}
\title{
	Plot a "hlaAttrBagObj" object
}
\description{
	To show a scatterplot of the numbers of individual classifiers and SNP positions.
}
\usage{
\method{plot}{hlaAttrBagObj}(x, xlab=NULL, ylab=NULL,
	locus.color="red", locus.lty=2, locus.cex=1.25, assembly="auto", ...)
}
\arguments{
	\item{x}{an object of \code{\link{hlaAttrBagObj}}}
	\item{xlab}{the label of X-axis}
	\item{ylab}{the label of Y-axis}
	\item{locus.color}{the color of text and line for HLA locus}
	\item{locus.lty}{the type of line for HLA locus}
	\item{locus.cex}{the font size of HLA locus}
	\item{assembly}{genome assembly: "hg19" (default), "hg18", "NCBI37" or "NCBI36";
		"auto" refers to the assembly used in \code{x}}
	\item{...}{further arguments passed to or from other methods}
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{plot.hlaAttrBagClass}}
	\code{\link{print.hlaAttrBagClass}}, \code{\link{print.hlaAttrBagObj}},
	\code{\link{summary.hlaAttrBagClass}}, \code{\link{summary.hlaAttrBagObj}}
}

\examples{
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "C"
hla <- hlaAllele(HLA_Type_Table$sample.id, HLA_Type_Table[, paste(hla.id, ".1", sep="")],
	HLA_Type_Table[, paste(hla.id, ".2", sep="")], locus=hla.id, assembly="hg19")

# training genotypes
region <- 100   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
	hla.id, region*1000, assembly="hg19")
train.geno <- hlaGenoSubset(HapMap_CEU_Geno, snp.sel=match(snpid, HapMap_CEU_Geno$snp.id))

# train a HIBAG model
set.seed(1000)
# please use "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hla, train.geno, nclassifier=2, verbose.detail=TRUE)
plot(model)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
