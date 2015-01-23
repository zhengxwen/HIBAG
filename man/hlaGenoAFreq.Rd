\name{hlaGenoAFreq}
\alias{hlaGenoAFreq}
\title{
    Allele Frequency
}
\description{
    To calculate the allele frequencies from genotypes or haplotypes.
}
\usage{
hlaGenoAFreq(obj)
}
\arguments{
    \item{obj}{an object of \code{\link{hlaSNPGenoClass}}}
}
\value{
    Return allele frequecies.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaGenoAFreq}}, \code{\link{hlaGenoMFreq}},
    \code{\link{hlaGenoMRate}}, \code{\link{hlaGenoMRate_Samp}}
}

\examples{
summary(HapMap_CEU_Geno)

summary(hlaGenoAFreq(HapMap_CEU_Geno))
}

\keyword{SNP}
\keyword{genetics}
