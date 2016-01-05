\name{hlaConvSequence}
\alias{hlaConvSequence}
\title{
    Conversion From HLA Alleles to Amino Acid Sequences
}
\description{
    Convert (P-coded or G-coded) HLA alleles to amino acid sequences.
}
\usage{
hlaConvSequence(hla=character(), locus=NULL,
    method=c("protein", "protein_reference"),
    code=c("exact", "P.code", "G.code", "P.code.merge", "G.code.merge"),
    region=c("auto", "all", "P.code", "G.code"), release=c("v3.22.0"),
    replace=NULL)
}
\arguments{
    \item{hla}{characters, or an object of \code{\link{hlaAlleleClass}},
        at least 4-digit or 2-field (P-coded) HLA alleles}
    \item{locus}{"A", "B", "C", "DRB1", "DQA1", "DQB1", "DPB1" or "DPA1"}
    \item{method}{"protein": returns protein sequence alignments,
        "protein_reference": returns the protein sequence alignment reference}
    \item{code}{"exact": requires full resolution; "P.code": allows ambiguous
        alleles according to P code; "G.code": allows ambiguous alleles
        according to G code; "P.code.merge" and "G.code.merge" merge multiple
        ambiguous allele sequences by masking unknown or ambiguous amino acid
        an asterisk}
    \item{region}{"all": returns all amino acid or nucleotide sequences;
        "P.code", "G.code": returns the exon 2 and 3 for HLA class I, and
        the exon 2 for HLA class II alleles; "auto": region="all" if
        code=="exact", region="P.code" if code=="P.code"|"P.code.merge",
        region="G.code" if code=="G.code"|"G.code.merge"}
    \item{release}{"v3.22.0" -- IPD-IMGT/HLA 3.22.0 database (2015-10-07)}
    \item{replace}{NULL, or a character vector, like \code{c("09:02"="107:01")},
        any "09:02" will be replaced by "107:01"}
}
\details{
    The P or G codes for reporting of ambiguous allele typings can be found:
\url{http://hla.alleles.org/alleles/p_groups.html} or
\url{http://hla.alleles.org/alleles/g_groups.html}. The protein sequences for
each HLA alleles could be found:
\url{http://hla.alleles.org/alleles/text_index.html}.

    Due to allelic ambiguity, multiple alleles are assigned to a 2-field
P-coded allele or 3-field G-coded allele. For HLA Class I alleles, identity
in the 'antigen binding domains' is based on identical protein sequences as
encoded by exons 2 and 3. For HLA Class II alleles this is based on identical
protein sequences as encoded by exon 2. P codes and G codes encode the same
protein sequence for the peptide binding domains (exon 2 and 3 for HLA class
I and exon 2 only for HLA class II alleles).

    1. the sequence is displayed as a hyphen "-" where it is identical to the
reference.

    2. an insertion or deletion is represented by a period ".".

    3. an unknown or ambiguous position in the alignment is represented by
an asterisk "*".

    4. a capital X is used for the 'stop' codons in protein alignments.

    \url{http://hla.alleles.org/alleles/formats.html}

    HLA class I and II sequence alignments (Text Index):
\url{http://hla.alleles.org/alleles/text_index.html}

    WARNING: if you are not familiar with HLA nomenclature, you might consult
with the package author or anyone who is familiar with HLA sequence alignments.
}
\value{
    Return an object of \code{\link{hlaAASeqClass}} or a list of characters.
\code{NULL} or \code{NA} in the list indicates no matching.
}
\references{
    The licence and disclaimer of distributed HLA data: Creative Commons
Attribution-NoDerivs Licence (\url{http://hla.alleles.org/terms.html}).

    Robinson J, Halliwell JA, Hayhurst JH, Flicek P, Parham P, Marsh SGE:
The IPD and IMGT/HLA database: allele variant databases.
Nucleic Acids Research. 2015 43:D423-431

    Robinson J, Malik A, Parham P, Bodmer JG, Marsh SGE:
IMGT/HLA - a sequence database for the human major histocompatibility complex.
Tissue Antigens. 2000 55:280-7
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaAlleleSubset}}
}

\examples{
hlaConvSequence(locus="A", method="protein_reference")
hlaConvSequence(c("01:01", "02:02", "01:01:01G", "01:01:01:01", "07"),
    locus="A")
hlaConvSequence(c("01:01", "02:02", "01:01:01G", "01:01:01:01", "07"),
    locus="A", code="P.code")
hlaConvSequence(c("01:01", "02:02", "01:01:01G", "01:01:01:01", "07"),
    locus="A", code="P.code.merge")


hlaConvSequence(locus="DPB1", method="protein_reference")
hlaConvSequence(c("09:01", "09:02"), locus="DPB1", replace=c("09:02"="107:01"))
hlaConvSequence(c("09:01", "09:02"), locus="DPB1", code="P.code",
    replace=c("09:02"="107:01"))
hlaConvSequence(c("09:01", "09:02"), locus="DPB1", code="P.code.merge",
    replace=c("09:02"="107:01"))


hlaConvSequence(locus="DQB1", method="protein_reference")
hlaConvSequence(c("05:01:01:01", "06:01:01"), locus="DQB1")
hlaConvSequence(c("05:01", "06:01"), locus="DQB1", code="P.code")
hlaConvSequence(c("05:01", "06:01"), locus="DQB1", code="P.code.merge")



hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
    H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
    H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
    locus=hla.id, assembly="hg19")

(v <- hlaConvSequence(hla, code="P.code.merge"))
summary(v)

v <- hlaConvSequence(hla, code="P.code.merge", region="all")
summary(v)



hla.id <- "DQB1"
hla <- hlaAllele(HLA_Type_Table$sample.id,
    H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
    H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
    locus=hla.id, assembly="hg19")

(v <- hlaConvSequence(hla, code="P.code.merge"))
summary(v)

v <- hlaConvSequence(hla, code="P.code.merge", region="all")
summary(v)
}

\keyword{HLA}
\keyword{SNP}
\keyword{genetics}
