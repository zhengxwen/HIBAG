\name{hlaGenoCombine}
\alias{hlaGenoCombine}
\title{
	Combine two genotypic data sets into one
}
\description{
	To combine two genotypic data sets into one dataset.
}
\usage{
hlaGenoCombine(geno1, geno2, match.pos=TRUE, allele.check=TRUE)
}
\arguments{
	\item{geno1}{the first genotype object of \code{\link{hlaSNPGenoClass}}}
	\item{geno2}{the second genotype object of \code{\link{hlaSNPGenoClass}}}
	\item{match.pos}{if \code{TRUE}, it uses both SNP ID and position to match SNPs;
		otherwise, it uses SNP ID only; calling \code{hlaSNPID(, with.pos=match.pos)}}
	\item{allele.check}{if \code{TRUE}, call \code{\link{hlaGenoSwitchStrand}} to check
		and then switch allele pairs if needed}
}
\details{
	The function merges two SNP dataset \code{geno1} and \code{geno2}, and returns a SNP
data consisting of the SNP intersect between \code{geno1} and \code{geno2}, and having
the same SNP information (allele and position) as \code{geno1}.
}
\value{
	A \code{\link{hlaSNPGenoClass}} object.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaMakeSNPGeno}}, \code{\link{hlaMakeSNPHaplo}},
	\code{\link{hlaHaploSubset}}, \code{\link{hlaGenoSubset}}
}

\examples{
# load SNP genotypes
data(HapMap_CEU_Geno, package="HIBAG")

# import a PLINK BED file
bed.fn <- system.file("extdata", "HapMap_CEU.bed", package="HIBAG")
fam.fn <- system.file("extdata", "HapMap_CEU.fam", package="HIBAG")
bim.fn <- system.file("extdata", "HapMap_CEU.bim", package="HIBAG")
hapmap.ceu <- hlaBED2Geno(bed.fn, fam.fn, bim.fn, assembly="hg19")

# combine two datasets together
geno <- hlaGenoCombine(HapMap_CEU_Geno, hapmap.ceu)
summary(geno)
#	390 samples X 1564 SNPs
#	SNPs range from 25877003bp to 33529555bp
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
