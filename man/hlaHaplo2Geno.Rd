\name{hlaHaplo2Geno}
\alias{hlaHaplo2Geno}
\title{
	Get a genotype object from a specified haplotype object
}
\description{
	To get a genotype object from a specified haplotype object by combining
SNP alleles.
}
\usage{
hlaHaplo2Geno(hapobj)
}
\arguments{
	\item{hapobj}{a haplotype object of \code{\link{hlaSNPHaploClass}}}
}
\details{
	\code{hapobj$haplotype} is a numeric matrix, with an entry value 0 standing
for B (ZERO A allele), 1 for A (ONE A allele) and others for missing values
(missing SNP alleles are usually set to be NA).
}
\value{
	Return a \code{\link{hlaSNPHaploClass}} object, and it is a list:
	\item{haplotype}{a haplotype  matrix, ``# of SNPs'' - by - ``2 x # of individuals''}
	\item{sample.id}{a vector of sample IDs}
	\item{snp.id}{a vector of SNP IDs}
	\item{snp.position}{a vector of SNP positions in basepair}
	\item{snp.allele}{a vector of characters with the format of
		``A allele/B allele''}
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaMakeSNPGeno}}, \code{\link{hlaMakeSNPHaplo}},
	\code{\link{hlaGenoSubset}}, \code{\link{hlaHaploSubset}}
}

\keyword{SNP}
\keyword{genetics}
