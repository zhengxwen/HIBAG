\name{hlaCombineAllele}
\alias{hlaCombineAllele}
\title{
	Combine two datasets of HLA types
}
\description{
	Get a subset of HLA types from an object of "hlaAlleleClass".
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
hla.id <- "C"
hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
	hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)
summary(hla)

subhla1 <- hlaAlleleSubset(hla,   1:100)
summary(subhla1)
subhla2 <- hlaAlleleSubset(hla, 201:300)
summary(subhla2)

H <- hlaCombineAllele(subhla1, subhla2)
summary(H)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
