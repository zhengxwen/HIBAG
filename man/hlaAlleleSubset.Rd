\name{hlaAlleleSubset}
\alias{hlaAlleleSubset}
\title{
	Get a subset of HLA types
}
\description{
	Get a subset of HLA types from an object of \code{\link{hlaAlleleClass}}.
}
\usage{
hlaAlleleSubset(hla, samp.sel=NULL)
}
\arguments{
	\item{hla}{an object of \code{\link{hlaAlleleClass}}}
	\item{samp.sel}{a logical vector, or an integer vector of indices}
}
\value{
	Return \code{\link{hlaAlleleClass}}.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAllele}}, \code{\link{hlaAlleleDigit}}
}

\examples{
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")

head(HLA_Type_Table)
dim(HLA_Type_Table)  # 60 13

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
	H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
	H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
	locus=hla.id, assembly="hg19")
summary(hla)

subhla <- hlaAlleleSubset(hla, 1:100)
summary(subhla)
}

\keyword{HLA}
\keyword{genetics}
