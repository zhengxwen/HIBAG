\name{hlaAllele}
\alias{hlaAllele}
\title{
    A list of HLA/KIR types
}
\description{
    Return an object of \code{\link{hlaAlleleClass}}, which contains
HLA/KIR types.
}
\usage{
hlaAllele(sample.id, H1, H2, max.resolution="", locus="any", assembly="auto",
    locus.pos.start=NA_integer_, locus.pos.end=NA_integer_, prob=NULL,
    na.rm=TRUE)
}
\arguments{
    \item{sample.id}{sample IDs}
    \item{H1}{a vector of HLA/KIR alleles}
    \item{H2}{a vector of HLA/KIR alleles}
    \item{max.resolution}{"2-digit", "1-field", "4-digit", "2-field", "6-digit",
        "3-field", "8-digit", "4-field", "allele", "protein", "full", "none",
        or "": "allele" = "2-digit"; "protein" = "4-digit"; "full", "none" or
        "" for no limit on resolution}
    \item{locus}{the name of HLA locus: "A", "B", "C", "DRB1", "DRB5",
        "DQA1", "DQB1", "DPB1", KIR locus, or "any", where "any" indicates any other
        multiallelic locus; see \code{\link{hlaLociInfo}} for possible locus names}
    \item{assembly}{the human genome reference: "hg18", "hg19" (default),
        "hg38"; "auto" refers to "hg19"; "auto-silent" refers to "hg19" without
        any warning}
    \item{locus.pos.start}{the starting position in basepair}
    \item{locus.pos.end}{the end position in basepair}
    \item{prob}{the probabilities assigned to the samples}
    \item{na.rm}{if TRUE, remove the samples without valid HLA types}
}
\details{
    The format of \code{H1} and \code{H2} is "allele group : different protein :
synonymous mutations in exons : synonymous mutations in introns"L,
where the suffix L is express level (N, null; L, low; S, secreted; A, aberrant;
Q: questionable). For example, "44:02:01:02L".
    If \code{max.resolution} is specified, the HLA alleles will be trimmed with
a possible maximum resolution.
}
\value{
    Return a \code{\link{hlaAlleleClass}} object, and it is a list:
    \item{locus}{HLA locus}
    \item{pos.start}{the starting position in basepair}
    \item{pos.end}{the end position in basepair}
    \item{value}{a data frame}
    \item{assembly}{the human genome reference, such like "hg19"}

    The component \code{value} includes:
    \item{sample.id}{sample ID}
    \item{allele1}{HLA allele}
    \item{allele2}{HLA allele}
    \item{prob}{the posterior probability}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaAlleleDigit}}, \code{\link{hlaAlleleSubset}},
    \code{\link{hlaLociInfo}}, \code{\link{hlaAlleleToVCF}}
}

\examples{
head(HLA_Type_Table)
dim(HLA_Type_Table)  # 60 13

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
    H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
    H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
    locus=hla.id, assembly="hg19")
summary(hla)


# encode other loci
hlaAllele("HD0010", "1", "2", locus="NewLocus")
}

\keyword{HLA}
\keyword{genetics}
