\name{hlaGenoSwitchStrand}
\alias{hlaGenoSwitchStrand}
\title{
	Allele switching
}
\description{
	Determine the ordered pair of A and B alleles, using the allele information
provided by \code{template}.
}
\usage{
hlaGenoSwitchStrand(target, template, match.pos=TRUE, verbose=TRUE)
}
\arguments{
	\item{target}{a genotype object of \code{\link{hlaSNPGenoClass}}}
	\item{template}{a genotype object of \code{\link{hlaSNPGenoClass}}}
	\item{match.pos}{if TRUE, it uses both SNP ID and position to match SNPs; otherwise,
		it uses SNP ID only; calling \code{hlaSNPID(, with.pos=match.pos)}}
	\item{verbose}{show information, if TRUE}
}
\details{
	The A/B pairs of \code{target} are determined using the information from \code{template}.
}
\value{
	Return a \code{\link{hlaSNPGenoClass}} object consisting of the SNP intersect
between \code{target} and \code{template}.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaMakeSNPGeno}}, \code{\link{hlaMakeSNPHaplo}}, \code{\link{hlaHaploSubset}},
\code{\link{hlaGenoSubset}}
}

\examples{
# load SNP genotypes
data(HapMap_CEU_Geno, package="HIBAG")
summary(HapMap_CEU_Geno)
# A/C A/G C/T G/T 
# 136 655 632 141 

# import a PLINK BED file
bed.fn <- system.file("extdata", "HapMap_CEU.bed", package="HIBAG")
fam.fn <- system.file("extdata", "HapMap_CEU.fam", package="HIBAG")
bim.fn <- system.file("extdata", "HapMap_CEU.bim", package="HIBAG")
hapmap.ceu <- hlaBED2Geno(bed.fn, fam.fn, bim.fn, assembly="hg19")
summary(hapmap.ceu)
#  A/C  A/G  A/T  C/G  C/T  G/T 
#  471 2102   83  134 2046  480 

# combine two datasets together
geno <- hlaGenoSwitchStrand(HapMap_CEU_Geno, hapmap.ceu)
summary(geno)
# There are 1564 SNPs in common.
# The allele pairs of 763 SNPs need to be switched.
# A/C A/G C/T G/T 
# 136 655 632 141 
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
