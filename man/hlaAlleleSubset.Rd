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
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaMakeSNPGeno}}, \code{\link{hlaMakeSNPHaplo}}, \code{\link{hlaHaploSubset}},
\code{\link{hlaGenoSubset}}
}

\examples{
# load HLA types and SNP genotypes
data(hla.type.table, package="HIBAG")

head(hla.type.table)
dim(hla.type.table)  # 60 13

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
	hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)
summary(hla)

subhla <- hlaAlleleSubset(hla, 1:100)
summary(subhla)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
