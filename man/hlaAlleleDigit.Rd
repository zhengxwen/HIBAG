\name{hlaAlleleDigit}
\alias{hlaAlleleDigit}
\title{
	XX
}
\description{
	XX
}
\usage{
hlaAlleleDigit(obj, max.resolution="4-digit")
}
\arguments{
	\item{obj}{should be a "hlaAlleleClass" object or characters}
	\item{max.resolution}{"2-digit", "4-digit", "6-digit", "8-digit", "allele",
		"protein", "2", "4", "6", "8", "full" or "": "allele" = "2-digit",
		"protein" = "4-digit", "full" and "" indicating no limit on resolution}
}
\details{
	Either \code{HLAtypes} or \code{H1} \code{H2} should be specified.
	The format of \code{HLAtypes} is "allele 1 / allele 2", e.g., "0512/0102".
	If \code{max.resolution} is specified, the HLA alleles will be trimmed with the maximum resolution.
}
\value{
	Return a \code{\link{hlaAlleleClass}} object, and it is a list:
	\item{locus}{HLA locus}
	\item{pos.start}{the starting position in basepair}
	\item{pos.end}{the end position in basepair}
	\item{value}{a data frame}

	The component \code{value} includes:
	\item{sample.id}{sample ID}
	\item{types}{the format is an unordered pair of "allele 1 / allele 2", e.g., "0512/0102"}
	\item{allele1}{HLA allele}
	\item{allele2}{HLA allele}
	\item{prob}{the posterior probability}
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
hla <- hlaAllele(hla.type.table$sample.id, H1=hla.type.table[, paste(hla.id, ".1", sep="")],
	H2=hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)
summary(hla)

hla2 <- hlaAlleleDigit(hla, "2-digit")
summary(hla2)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
