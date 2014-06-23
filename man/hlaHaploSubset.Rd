\name{hlaHaploSubset}
\alias{hlaHaploSubset}
\title{
	Get a subset of haplotypes
}
\description{
	To get a subset of haplotypes from a \code{\link{hlaSNPHaploClass}} object.
}
\usage{
hlaHaploSubset(haploobj, samp.sel=NULL, snp.sel=NULL)
}
\arguments{
	\item{haploobj}{a haplotype object of \code{\link{hlaSNPHaploClass}}}
	\item{samp.sel}{a logical vector, or an integer vector of indices}
	\item{snp.sel}{a logical vector, or an integer vector of indices}
}
\details{
	\code{haploobj$haplotype} is a numeric matrix, with an entry value 0
standing for B (ZERO A allele), 1 for A (ONE A allele) and others for missing
values (missing genotypes are usually set to be NA).
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
	\code{\link{hlaGenoSubset}}, \code{\link{hlaGenoCombine}}
}

\keyword{SNP}
\keyword{genetics}
