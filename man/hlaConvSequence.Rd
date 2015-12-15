\name{hlaConvSequence}
\alias{hlaConvSequence}
\title{
    Conversion From HLA Alleles to Amino Acid Sequences
}
\description{
    Convert (P-coded or G-coded) HLA alleles to Amino Acid sequences.
}
\usage{
hlaConvSequence(hla, locus=NULL, method=c("AminoAcid_v3.22"), replace=NULL)
}
\arguments{
    \item{hla}{characters, or an object of \code{\link{hlaAlleleClass}}
        including at least 4-digit or 2-field (P-coded) HLA alleles}
    \item{locus}{"A", "B", "C", "DRB1", "DQA1", "DQB1", "DPB1" or "DPA1"}
    \item{method}{"AminoAcid_v3.22" -- output amino acids using
        IPD-IMGT/HLA 3.22.0 database (2015-10-07)}
    \item{replace}{NULL, or a character vector, like \code{c("09:02"="107:01")},
        any "09:02" will be replaced by "107:01"}
}
\details{
    The P or G codes for reporting of ambiguous allele typings can be found:
\url{http://hla.alleles.org/alleles/p_groups.html} or
\url{http://hla.alleles.org/alleles/g_groups.html}. The amino acids for each
HLA alleles could be found, e.g., 
\url{http://www.ebi.ac.uk/cgi-bin/ipd/imgt/hla/get_allele.cgi?A*01:02}.
Due to allelic ambiguity, multiple higher-resolution alleles are assigned to a
2-field P-coded allele. If the allele frequency is applicable from previews
study reports (\url{http://www.allelefrequencies.net}), the most possible
allele (amino acid sequences) is used for the output. Otherwise, the amino
acid sequences of the first higher-resolution allele is employed.
}
\value{
    Return an object of \code{\link{hlaSeqClass}} or characters.
\code{NA} indicates no matching.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaAlleleSubset}}
}

\examples{
# P and G codes
fn <- system.file("extdata", "Sequence", "HLA_Nom_v3.22.RData", package="HIBAG")
HLA_Nom <- get(load(fn))

HLA_Nom$Version  # "IPD-IMGT/HLA 3.22.0"
head(HLA_Nom$PCode, n=1L)  # Code / Alleles
head(HLA_Nom$GCode, n=1L)  # Code / Alleles


# amino acids
hlaConvSequence(c("01:01", "02:02", "01:01:01G", "01:01:01:01", "07"),
    locus="A")
hlaConvSequence(c("09:01", "09:02"), locus="DPB1", replace=c("09:02"="107:01"))

# make a "hlaAlleleClass" object
hla.id <- "C"
hla <- hlaAllele(HLA_Type_Table$sample.id,
    H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
    H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
    locus=hla.id, assembly="hg19")

summary(hlaConvSequence(hla))
}

\keyword{HLA}
\keyword{SNP}
\keyword{genetics}
