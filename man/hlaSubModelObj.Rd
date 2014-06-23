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
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAttrBagging}}
}

\examples{
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "C"
hla <- hlaAllele(HLA_Type_Table$sample.id,
	H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
	H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
	locus=hla.id, assembly="hg19")

# training genotypes
region <- 50   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
	hla.id, region*1000, assembly="hg19")
train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
	snp.sel = match(snpid, HapMap_CEU_Geno$snp.id))

# train a HIBAG model
set.seed(1000)
# please use "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hla, train.geno, nclassifier=2, verbose.detail=TRUE)
mobj <- hlaModelToObj(model)
summary(mobj)

newmobj <- hlaSubModelObj(mobj, 1)
summary(newmobj)
}

\keyword{HLA}
\keyword{genetics}
