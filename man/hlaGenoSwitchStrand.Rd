\name{hlaGenoSwitchStrand}
\alias{hlaGenoSwitchStrand}
\title{
	Determine the ordered pair of A and B alleles
}
\description{
	To get a genotype object, using the allele information provided by \code{template}.
}
\usage{
hlaGenoSwitchStrand(target, template, verbose=TRUE)
}
\arguments{
	\item{target}{a genotype object of \code{\link{hlaSNPGenoClass}}}
	\item{template}{a genotype object of \code{\link{hlaSNPGenoClass}}}
	\item{verbose}{show information, if TRUE}
}
\details{
	The A/B pairs of \code{target} are determined using the information from \code{template}.
}
\value{
	Return a \code{\link{hlaSNPGenoClass}} object.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaMakeSNPGeno}}, \code{\link{hlaMakeSNPHaplo}}, \code{\link{hlaHaploSubset}},
\code{\link{hlaGenoSubset}}
}

\examples{
# load SNP genotypes
data(hla.snp.geno, package="HIBAG")
summary(hla.snp.geno)
# A/C A/G C/A C/T G/A G/T T/C T/G 
#  64 369  64 305 276  69 337  80 

# import a PLINK BED file
bed.fn <- system.file("extdata", "hapmap.ceu.bed", package="HIBAG")
fam.fn <- system.file("extdata", "hapmap.ceu.fam", package="HIBAG")
bim.fn <- system.file("extdata", "hapmap.ceu.bim", package="HIBAG")
hapmap.ceu <- hlaBED2Geno(bed.fn, fam.fn, bim.fn)
summary(hapmap.ceu)
#  A/C  A/G  A/T  C/G  C/T  G/T 
#  471 2102   83  134 2046  480 

# combine two datasets together
geno <- hlaGenoSwitchStrand(hla.snp.geno, hapmap.ceu)
summary(geno)
# There are 1564 SNPs in common.
# The allele pairs of 763 SNPs need to be switched.
# A/C A/G C/T G/T 
# 136 655 632 141 
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
