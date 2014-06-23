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
	\item{obj}{an object of \code{\link{hlaSNPGenoClass}} or \code{\link{hlaSNPHaploClass}}}
}
\value{
	Return missing rates per sample.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaGenoAFreq}}, \code{\link{hlaGenoMFreq}}, \code{\link{hlaGenoMRate}},
\code{\link{hlaGenoMRate_Samp}}
}

\examples{
# load SNP genotypes
data(hla.snp.geno, package="HIBAG")

summary(hla.snp.geno)

summary(hlaGenoMRate_Samp(hla.snp.geno))
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
