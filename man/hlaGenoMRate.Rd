\name{hlaGenoMRate}
\alias{hlaGenoMRate}
\title{
    Missing Rates Per SNP
}
\description{
    To calculate the missing rates from genotypes or haplotypes per SNP.
}
\usage{
hlaGenoMRate(obj)
}
\arguments{
    \item{obj}{an object of \code{\link{hlaSNPGenoClass}}}
}
\value{
    Return missing rates per SNP.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaGenoAFreq}}, \code{\link{hlaGenoMFreq}},
    \code{\link{hlaGenoMRate}}, \code{\link{hlaGenoMRate_Samp}}
}

\examples{
summary(HapMap_CEU_Geno)

summary(hlaGenoMRate(HapMap_CEU_Geno))
}

\keyword{SNP}
\keyword{genetics}
