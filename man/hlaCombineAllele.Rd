\name{hlaCombineAllele}
\alias{hlaCombineAllele}
\title{
	Combine two datasets of HLA types
}
\description{
	Get a subset of HLA types from an object of \code{\link{hlaAlleleClass}}.
}
\usage{
hlaCombineAllele(H1, H2)
}
\arguments{
	\item{H1}{the first \code{\link{hlaAlleleClass}} object}
	\item{H2}{the second \code{\link{hlaAlleleClass}} object}
}
\value{
	Return \code{\link{hlaAlleleClass}}.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAllele}}, \code{\link{hlaAlleleSubset}}
}

\examples{
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")

head(HLA_Type_Table)
dim(HLA_Type_Table)  # 60 13

# make a "hlaAlleleClass" object
hla.id <- "C"
hla <- hlaAllele(HLA_Type_Table$sample.id, HLA_Type_Table[, paste(hla.id, ".1", sep="")],
	HLA_Type_Table[, paste(hla.id, ".2", sep="")], locus=hla.id, assembly="hg19")
summary(hla)

subhla1 <- hlaAlleleSubset(hla,   1:100)
summary(subhla1)
subhla2 <- hlaAlleleSubset(hla, 201:300)
summary(subhla2)

H <- hlaCombineAllele(subhla1, subhla2)
summary(H)
}

\keyword{HLA}
\keyword{genetics}
