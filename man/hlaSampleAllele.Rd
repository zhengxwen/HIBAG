\name{hlaSampleAllele}
\alias{hlaSampleAllele}
\title{
    Get sample IDs from HLA types with a filter
}
\description{
    Get sample IDs from HLA types limited to a set of HLA alleles.
}
\usage{
hlaSampleAllele(TrueHLA, allele.limit=NULL, max.resolution="")
}
\arguments{
    \item{TrueHLA}{an object of \code{\link{hlaAlleleClass}}}
    \item{allele.limit}{a list of HLA alleles, the validation samples are
        limited to those having HLA alleles in \code{allele.limit}, or
        \code{NULL} for no limit. \code{allele.limit} could be character-type,
        \code{\link{hlaAttrBagClass}} or \code{\link{hlaAttrBagObj}}}
    \item{max.resolution}{"2-digit", "4-digit", "6-digit", "8-digit", "allele",
        "protein", "2", "4", "6", "8", "full" or "": "allele" = "2-digit",
        "protein" = "4-digit", "full" and "" mean no limit on resolution}
}
\value{
    Return a list of sample IDs.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaCompareAllele}}
}

\examples{
# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
    HLA_Type_Table[, paste(hla.id, ".1", sep="")],
    HLA_Type_Table[, paste(hla.id, ".2", sep="")],
    locus=hla.id, assembly="hg19")
summary(hla)

hlaSampleAllele(hla)

hlaSampleAllele(hla, allele.limit=c(
    "01:01","02:01","02:06", "03:01", "11:01", "23:01"))
}

\keyword{HLA}
\keyword{genetics}
