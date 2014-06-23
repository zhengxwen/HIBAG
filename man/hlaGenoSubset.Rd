\name{hlaGenoSubset}
\alias{hlaGenoSubset}
\title{
	Get a subset of genotypes
}
\description{
	To get a subset of genotypes from a \code{\link{hlaSNPGenoClass}} object.
}
\usage{
hlaGenoSubset(genoobj, samp.sel=NULL, snp.sel=NULL)
}
\arguments{
	\item{genoobj}{a genotype object of \code{\link{hlaSNPGenoClass}}}
	\item{samp.sel}{a logical vector, or an integer vector of indices}
	\item{snp.sel}{a logical vector, or an integer vector of indices}
}
\details{
	\code{genoobj$genotype} is a numeric matrix, with an entry value 0 standing
for BB (ZERO A allele), 1 for AB (ONE A allele), 2 for AA (TWO A alleles) and
others for missing values (missing genotypes are usually set to be NA).
}
\value{
	Return a \code{\link{hlaSNPGenoClass}} object, and it is a list:
	\item{genotype}{a genotype matrix, ``# of SNPs'' - by - ``# of individuals''}
	\item{sample.id}{a vector of sample IDs}
	\item{snp.id}{a vector of SNP IDs}
	\item{snp.position}{a vector of SNP positions in basepair}
	\item{snp.allele}{a vector of characters with the format of
		``A allele/B allele''}
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaMakeSNPGeno}}, \code{\link{hlaMakeSNPHaplo}},
	\code{\link{hlaHaploSubset}}, \code{\link{hlaGenoCombine}}
}

\examples{
# load SNP genotypes
data(HapMap_CEU_Geno, package="HIBAG")
summary(HapMap_CEU_Geno)

geno <- hlaGenoSubset(HapMap_CEU_Geno,
	snp.sel = (hlaGenoMFreq(HapMap_CEU_Geno)>0.10))
summary(geno)
}

\keyword{SNP}
\keyword{genetics}
