\name{summary.hlaAlleleClass}
\alias{summary.hlaAlleleClass}
\title{
	Get a subset of HLA types
}
\description{
	Get a subset of HLA types from an object of \code{\link{hlaAlleleClass}}.
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
\author{Xiuwen Zheng}
\seealso{
	\code{\link{summary.hlaSNPGenoClass}}
}

\examples{
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id, HLA_Type_Table[, paste(hla.id, ".1", sep="")],
	HLA_Type_Table[, paste(hla.id, ".2", sep="")], locus=hla.id, assembly="hg19")
summary(hla)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
