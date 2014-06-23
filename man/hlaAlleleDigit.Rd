\name{hlaAlleleDigit}
\alias{hlaAlleleDigit}
\title{
	Trim HLA alleles
}
\description{
	Trim HLA alleles to specified width.
}
\usage{
hlaAlleleDigit(obj, max.resolution="4-digit", rm.suffix=FALSE)
}
\arguments{
	\item{obj}{should be a \code{\link{hlaAlleleClass}} object or characters}
	\item{max.resolution}{"2-digit", "4-digit", "6-digit", "8-digit", "allele",
		"protein", "2", "4", "6", "8", "full" or "": "allele" = "2-digit",
		"protein" = "4-digit", "full" and "" indicating no limit on resolution}
	\item{rm.suffix}{whether remove the suffix, e.g., for "01:22N", "N" is a
		suffix}
}
\details{
	Either \code{HLAtypes} or \code{H1} \code{H2} should be specified. The
format of \code{HLAtypes} is "allele 1 / allele 2", e.g., "0512/0102". If
\code{max.resolution} is specified, the HLA alleles will be trimmed with the
maximum resolution.
}
\value{
	Return a \code{\link{hlaAlleleClass}} object if \code{obj} is
\code{\link{hlaAlleleClass}}-type, or characters if \code{obj} is
character-type.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAllele}}
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
	locus = hla.id, assembly="hg19")
summary(hla)

hla2 <- hlaAlleleDigit(hla, "2-digit")
summary(hla2)
}

\keyword{HLA}
\keyword{genetics}
