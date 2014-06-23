\name{hlaUniqueAllele}
\alias{hlaUniqueAllele}
\title{
	Get unique HLA alleles
}
\description{
	Get unique HLA alleles, which are in ascending order.
}
\usage{
hlaUniqueAllele(hla)
}
\arguments{
	\item{hla}{character-type HLA alleles, or a
		\code{\link{hlaAlleleClass}} object}
}
\details{
	Each HLA allele name has a unique number corresponding to up to four sets
of digits separated by colons. The name designation depends on the sequence of
the allele and that of its nearest relative. The digits before the first colon
describe the type, which often corresponds to the serological antigen carried
by an allotype. The next set of digits are used to list the subtypes, numbers
being assigned in the order in which DNA sequences have been determined.
Alleles whose numbers differ in the two sets of digits must differ in one or
more nucleotide substitutions that change the amino acid sequence of the
encoded protein. Alleles that differ only by synonymous nucleotide substitutions
(also called silent or non-coding substitutions) within the coding sequence
are distinguished by the use of the third set of digits. Alleles that only
differ by sequence polymorphisms in the introns or in the 5' or 3' untranslated
regions that flank the exons and introns are distinguished by the use of the
fourth set of digits.

	In addition to the unique allele number there are additional optional
suffixes that may be added to an allele to indicate its expression status.
Alleles that have been shown not to be expressed, 'Null' alleles have been
given the suffix 'N'. Those alleles which have been shown to be alternatively
expressed may have the suffix 'L', 'S', 'C', 'A' or 'Q'.

	\url{http://hla.alleles.org/nomenclature/index.html}
}
\value{
	Return a vector of HLA alleles
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAllele}}, \code{\link{hlaAlleleDigit}}
}

\examples{
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
	H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
	H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
	locus=hla.id, assembly="hg19")
summary(hla)
hlaUniqueAllele(hla)

hlaUniqueAllele(c("01", "01:03", "01:01", "03:05", "03:01g", "03:01", "104:01"))
}

\keyword{HLA}
\keyword{genetics}
