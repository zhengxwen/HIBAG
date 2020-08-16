\name{summary.hlaSNPGenoClass}
\alias{summary.hlaSNPGenoClass}
\alias{print.hlaSNPGenoClass}
\title{
    Summarize a SNP dataset
}
\description{
    Summarize the genotypic dataset.
}
\usage{
\method{summary}{hlaSNPGenoClass}(object, show=TRUE, ...)
\method{print}{hlaSNPGenoClass}(x, ...)
}
\arguments{
    \item{object}{a genotype object of \code{\link{hlaSNPGenoClass}}}
    \item{x}{a genotype object of \code{\link{hlaSNPGenoClass}}}
    \item{show}{if TRUE, print information}
    \item{...}{further arguments passed to or from other methods}
}
\value{
    None.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaMakeSNPGeno}}, \code{\link{hlaGenoSubset}}
}

\examples{
summary(HapMap_CEU_Geno)
}

\keyword{SNP}
\keyword{genetics}
