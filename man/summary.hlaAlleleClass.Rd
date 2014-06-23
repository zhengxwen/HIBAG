\name{summary.hlaAlleleClass}
\alias{summary.hlaAlleleClass}
\title{
	Get a subset of HLA types
}
\description{
	Get a subset of HLA types from an object of "hlaAlleleClass".
}
\usage{
\method{summary}{hlaAlleleClass}(object, show=TRUE, ...)
}
\arguments{
	\item{object}{an object of \code{\link{hlaAlleleClass}}}
	\item{show}{if TRUE, show information}
	\item{...}{further arguments passed to or from other methods}
}
\value{
	Return a \code{data.frame} of count and frequency for each HLA allele.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{summary.hlaSNPGenoClass}}
}

\examples{
# load HLA types and SNP genotypes
data(hla.type.table, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
	hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)
summary(hla)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
