\name{hlaGenoMRate_Samp}
\alias{hlaGenoMRate_Samp}
\title{
	Missing Rates Per Sample
}
\description{
	To calculate the missing rates from genotypes or haplotypes per sample.
}
\usage{
hlaGenoMRate_Samp(obj)
}
\arguments{
	\item{obj}{an object of \code{\link{hlaSNPGenoClass}} or
		\code{\link{hlaSNPHaploClass}}}
}
\value{
	Return missing rates per sample.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaGenoAFreq}}, \code{\link{hlaGenoMFreq}},
	\code{\link{hlaGenoMRate}}, \code{\link{hlaGenoMRate_Samp}}
}

\examples{
# load SNP genotypes
data(HapMap_CEU_Geno, package="HIBAG")

summary(HapMap_CEU_Geno)

summary(hlaGenoMRate_Samp(HapMap_CEU_Geno))
}

\keyword{SNP}
\keyword{genetics}
