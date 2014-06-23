\name{hlaSplitAllele}
\alias{hlaSplitAllele}
\title{
	Divide the samples randomly
}
\description{
	Divide the samples to the training and validation sets randomly.
}
\usage{
hlaSplitAllele(HLA, train.prop=0.5)
}
\arguments{
	\item{HLA}{an object of \code{\link{hlaAlleleClass}}}
	\item{train.prop}{the proporion of training set}
}
\details{
	The algorithm tries to divide each HLA alleles into training and validation sets randomly
with a training proportion "\code{train.prop}".
}
\value{
	Return a list:
	\item{training}{an object of \code{\link{hlaAlleleClass}}}
	\item{validation}{an object of \code{\link{hlaAlleleClass}}}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaMakeSNPGeno}}, \code{\link{hlaMakeSNPHaplo}}, \code{\link{hlaHaploSubset}},
\code{\link{hlaGenoSubset}}
}

\examples{
# load HLA types and SNP genotypes
data(hla.type.table, package="HIBAG")
data(hla.snp.geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
	hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)

# divide HLA types randomly
set.seed(100)
hlatab <- hlaSplitAllele(hla, train.prop=0.5)
names(hlatab)
# "training"   "validation"
summary(hlatab$training)
summary(hlatab$validation)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
