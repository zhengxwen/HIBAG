\name{hlaSNPHaploClass}
\alias{hlaSNPHaploClass}
\title{The class of SNP haplotypes}
\description{
	The class of SNP haplotypes, and its instance is returned from
\code{\link{hlaMakeSNPHaplo}}.
}
\value{
	There are five components:
	\item{haplotype}{a haplotype  matrix, ``# of SNPs'' - by - ``2 x # of individuals''}
	\item{sample.id}{a vector of sample IDs}
	\item{snp.id}{a vector of SNP IDs}
	\item{snp.position}{a vector of SNP positions in basepair}
	\item{snp.allele}{a vector of characters with the format of ``A allele/B allele''}
}

\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaMakeSNPHaplo}}
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
