\name{hlaSNPID}
\alias{hlaSNPID}
\title{
	Get SNP IDs and positions
}
\description{
	Get the information of SNP ID and position
}
\usage{
hlaSNPID(obj)
}
\arguments{
	\item{obj}{a genotype object of \code{\link{hlaSNPGenoClass}},
		a haplotype object of \code{\link{hlaSNPHaploClass}},
		an object of \code{\link{hlaAttrBagClass}}, or
		an object of \code{\link{hlaAttrBagObj}}}
}
\value{
	return \code{paste(obj$snp.id, obj$snp.position, sep="-")}.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaGenoSwitchStrand}}, \code{\link{hlaGenoCombine}}
}

\examples{
# load SNP genotypes
data(hla.snp.geno, package="HIBAG")

x <- hlaSNPID(hla.snp.geno)
head(x)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
