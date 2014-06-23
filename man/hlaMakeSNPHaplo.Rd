\name{hlaMakeSNPHaplo}
\alias{hlaMakeSNPHaplo}
\title{
	Make a SNP haplotype object
}
\description{
	To create a \code{\link{hlaSNPHaploClass}} object (SNP haplotypic object).
}
\usage{
hlaMakeSNPHaplo(haplotype, sample.id, snp.id, snp.position,
	A.allele, B.allele, assembly=c("auto", "hg18", "hg19", "unknown"))
}
\arguments{
	\item{haplotype}{a haplotype matrix, ``# of SNPs'' - by - ``2 x # of individuals''}
	\item{sample.id}{a vector of sample ids}
	\item{snp.id}{a vector of SNP ids}
	\item{snp.position}{a vector of SNP positions}
	\item{A.allele}{a vector of A alleles in the SNP list}
	\item{B.allele}{a vector of B alleles in the SNP list}
	\item{assembly}{the human genome reference: "hg19" (default), "hg18",
		"auto" refers to "hg19"}
}
\details{
	\code{haplotype} is a numeric matrix, with an entry value 0 standing for
B (ZERO A allele), 1 for A (ONE A allele) and others for missing values
(missing genotypes are usually set to be NA).
}
\value{
	Return a \code{\link{hlaSNPHaploClass}} object, and it is a list:
	\item{haplotype}{a haplotype  matrix, ``# of SNPs'' - by - ``2 x # of individuals''}
	\item{sample.id}{a vector of sample IDs}
	\item{snp.id}{a vector of SNP IDs}
	\item{snp.position}{a vector of SNP positions in basepair}
	\item{snp.allele}{a vector of characters with the format of
		``A allele/B allele''}
	\item{assembly}{the human genome reference}
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaMakeSNPGeno}}, \code{\link{hlaGenoSubset}},
	\code{\link{hlaHaploSubset}}
}

\keyword{SNP}
\keyword{genetics}
