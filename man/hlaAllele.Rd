\name{hlaAllele}
\alias{hlaAllele}
\title{
	A list of HLA types
}
\description{
	Return an object of "hlaAlleleClass", which contains HLA types.
}
\usage{
hlaAllele(sample.id, H1, H2, max.resolution="",
	locus="", locus.pos.start=NA, locus.pos.end=NA, prob=NULL, na.rm=TRUE)
}
\arguments{
	\item{sample.id}{sample IDs}
	\item{H1}{a vector of HLA alleles}
	\item{H2}{a vector of HLA alleles}
	\item{max.resolution}{"2-digit", "4-digit", "6-digit", "8-digit", "allele",
		"protein", "2", "4", "6", "8", "full" or "": "allele" = "2-digit",
		"protein" = "4-digit", "full" and "" indicating no limit on resolution}
	\item{locus}{the name of HLA locus: "A", "B", "C", "DRB1", "DRB5", "DQA1",
		"DQB1", "DPB1", or "any", where "any" indicates any other multiallelic locus}
	\item{locus.pos.start}{the starting position in basepair}
	\item{locus.pos.end}{the end position in basepair}
	\item{prob}{the probabilities assigned to the samples}
	\item{na.rm}{if TRUE, remove the samples without valid HLA types}
}
\details{
	The format of \code{H1} and \code{H2} is "allele group : different protein :
synonymous mutations in exons : synonymous mutations in introns"L,
where the suffix L is express level (N, null; L, low; S, secreted; A, aberrant; Q: questionable).
For example, "44:02:01:02L".
	If \code{max.resolution} is specified, the HLA alleles will be trimmed with
a possible maximum resolution.
}
\value{
	Return a \code{\link{hlaAlleleClass}} object, and it is a list:
	\item{locus}{HLA locus}
	\item{pos.start}{the starting position in basepair}
	\item{pos.end}{the end position in basepair}
	\item{value}{a data frame}

	The component \code{value} includes:
	\item{sample.id}{sample ID}
	\item{allele1}{HLA allele}
	\item{allele2}{HLA allele}
	\item{prob}{the posterior probability}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaMakeSNPGeno}}, \code{\link{hlaMakeSNPHaplo}},
	\code{\link{hlaHaploSubset}}, \code{\link{hlaGenoSubset}}
}

\examples{
# load HLA types and SNP genotypes
data(hla.type.table, package="HIBAG")

head(hla.type.table)
dim(hla.type.table)  # 60 13

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
	hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)
summary(hla)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
