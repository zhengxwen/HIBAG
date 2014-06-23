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
	The algorithm tries to divide each HLA alleles into training and validation
sets randomly with a training proportion \code{train.prop}.
}
\value{
	Return a list:
	\item{training}{an object of \code{\link{hlaAlleleClass}}}
	\item{validation}{an object of \code{\link{hlaAlleleClass}}}
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAllele}}
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
}

\keyword{HLA}
\keyword{genetics}
