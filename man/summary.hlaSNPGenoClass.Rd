\name{summary.hlaSNPGenoClass}
\alias{summary.hlaSNPGenoClass}
\title{
	Summarize the genotypic dataset
}
\description{
	Summarize the genotypic dataset.
}
\usage{
\method{summary}{hlaSNPGenoClass}(object, show=TRUE, ...)
}
\arguments{
	\item{object}{a genotype object of \code{\link{hlaSNPGenoClass}}}
	\item{show}{if TRUE, print information}
	\item{...}{further arguments passed to or from other methods}
}
\value{
	None.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaMakeSNPGeno}}, \code{\link{hlaMakeSNPHaplo}},
	\code{\link{hlaGenoSubset}}, \code{\link{hlaHaploSubset}},
	\code{\link{summary.hlaSNPHaploClass}}
}

\examples{
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

summary(HapMap_CEU_Geno)
}

\keyword{SNP}
\keyword{genetics}
