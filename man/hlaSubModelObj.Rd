\name{hlaSubModelObj}
\alias{hlaSubModelObj}
\title{
	Get a subset of individual classifiers
}
\description{
	Get the first n individual classifiers.
}
\usage{
hlaSubModelObj(obj, n)
}
\arguments{
	\item{obj}{an object of \code{\link{hlaAttrBagObj}}}
	\item{n}{an integer, get the first n individual classifiers}
}
\value{
	Return an object of \code{\link{hlaAttrBagObj}}.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaAttrBagging}}
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
region <- 50   # kb
snpid <- hlaFlankingSNP(hla.snp.geno$snp.id, hla.snp.geno$snp.position, hla.id, region*1000)
train.geno <- hlaGenoSubset(hla.snp.geno, snp.sel=match(snpid, hla.snp.geno$snp.id))

# train an attribute bagging model
set.seed(1000)
# please try "nclassifier=25" or "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hla, train.geno, nclassifier=2, verbose.detail=TRUE)
mobj <- hlaModelToObj(model)
summary(mobj)

newmobj <- hlaSubModelObj(mobj, 1)
summary(newmobj)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
