\name{hlaSNPGenoClass}
\alias{hlaSNPGenoClass}
\title{The class of SNP genotypes}
\description{
	The class of SNP genotypes, and its instance is returned from
\code{\link{hlaMakeSNPGeno}}.
}
\value{
	There are five components:
	\item{genotype}{a genotype matrix, ``# of SNPs'' - by - ``# of individuals''}
	\item{sample.id}{a vector of sample IDs}
	\item{snp.id}{a vector of SNP IDs}
	\item{snp.position}{a vector of SNP positions in basepair}
	\item{snp.allele}{a vector of characters with a format of
		``A allele/B allele''}
	\item{assembly}{the human genome reference, such like "hg19"}
}

\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaMakeSNPGeno}}
}

\keyword{SNP}
\keyword{genetics}
