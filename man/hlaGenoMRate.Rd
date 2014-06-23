\name{hlaGenoMRate}
\alias{hlaGenoMRate}
\title{
	Missing Rates Per SNP
}
\description{
	To calculate the missing rates from genotypes or haplotypes per SNP
}
\usage{
hlaGenoMRate(obj)
}
\arguments{
	\item{obj}{an object of \code{\link{hlaSNPGenoClass}} or \code{\link{hlaSNPHaploClass}}}
}
\value{
	Return missing rates per SNP.
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

summary(hlaGenoMRate(hla.snp.geno))
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
