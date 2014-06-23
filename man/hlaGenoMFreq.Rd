\name{hlaGenoMFreq}
\alias{hlaGenoMFreq}
\title{
	Minor Allele Frequency
}
\description{
	To calculate the minor allele frequencies from genotypes or haplotypes
}
\usage{
hlaGenoMFreq(obj)
}
\arguments{
	\item{obj}{an object of \code{\link{hlaSNPGenoClass}} or \code{\link{hlaSNPHaploClass}}}
}
\value{
	Return minor allele frequecies.
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

summary(hlaGenoMFreq(hla.snp.geno))
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
