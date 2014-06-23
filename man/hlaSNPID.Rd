\name{hlaSNPID}
\alias{hlaSNPID}
\title{
	Get SNP IDs and positions
}
\description{
	Get the information of SNP ID with or without position
}
\usage{
hlaSNPID(obj, with.pos=TRUE)
}
\arguments{
	\item{obj}{a genotype object of \code{\link{hlaSNPGenoClass}},
		a haplotype object of \code{\link{hlaSNPHaploClass}},
		an object of \code{\link{hlaAttrBagClass}}, or
		an object of \code{\link{hlaAttrBagObj}}}
	\item{with.pos}{if TRUE, return position information; otherwise, no position information}
}
\value{
	if \code{with.pos = TRUE},
		return \code{paste(obj$snp.id, obj$snp.position, sep="-")},
	otherwise,
		return \code{obj$snp.id} only.
}
\details{
	This function can be used to distinguish SNPs with both of SNP ID and position.
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
