\name{hlaAlleleDigit}
\alias{hlaAlleleDigit}
\title{
    Trim HLA alleles
}
\description{
    Trim HLA alleles to specified width.
}
\usage{
hlaAlleleDigit(obj, max.resolution=NA_character_, rm.suffix=FALSE)
}
\arguments{
    \item{obj}{should be a \code{\link{hlaAlleleClass}} object or characters}
    \item{max.resolution}{"2-digit", "1-field", "4-digit", "2-field", "6-digit",
        "3-field", "8-digit", "4-field", "allele", "protein", "full", "none",
        or "": "allele" = "2-digit"; "protein" = "4-digit"; "full", "none" or
        "" for no limit on resolution}
    \item{rm.suffix}{whether remove the non-digit suffix in the last field,
        e.g., for "01:22N", "N" is a non-digit suffix}
}
\details{
    If \code{max.resolution} is specified, the HLA alleles will be trimmed with
the maximum resolution.
See \url{https://hla.alleles.org/nomenclature/naming.html} for the HLA
nomenclature.
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
