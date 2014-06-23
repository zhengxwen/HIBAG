\name{plot.hlaAttrBagClass}
\alias{plot.hlaAttrBagClass}
\title{
	Plot a "hlaAttrBagClass" object
}
\description{
	To show a scatterplot of the numbers of individual classifiers and SNP positions
}
\usage{
\method{plot}{hlaAttrBagClass}(x, ...)
}
\arguments{
	\item{x}{an object of \code{\link{hlaAlleleClass}}}
	\item{...}{further arguments passed to or from other methods}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{plot.hlaAttrBagObj}},
	\code{\link{print.hlaAttrBagClass}}, \code{\link{print.hlaAttrBagObj}},
	\code{\link{summary.hlaAttrBagClass}}, \code{\link{summary.hlaAttrBagObj}}
}

\examples{
# load HLA types and SNP genotypes
data(hla.type.table, package="HIBAG")
data(hla.snp.geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "C"
hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
	hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)

# training genotypes
region <- 100   # kb
snpid <- hlaFlankingSNP(hla.snp.geno$snp.id, hla.snp.geno$snp.position, hla.id, region*1000)
train.geno <- hlaGenoSubset(hla.snp.geno, snp.sel=match(snpid, hla.snp.geno$snp.id))

# train an attribute bagging model
set.seed(1000)
# please try "nclassifier=25" or "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hla, train.geno, nclassifier=2, verbose.detail=TRUE)
plot(model)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
