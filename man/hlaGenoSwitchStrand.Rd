\name{hlaGenoSwitchStrand}
\alias{hlaGenoSwitchStrand}
\title{
	Allele switching
}
\description{
	Determine the ordered pair of A and B alleles, using the allele information provided by \code{template}.
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
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaMakeSNPGeno}}, \code{\link{hlaMakeSNPHaplo}}, \code{\link{hlaHaploSubset}},
\code{\link{hlaGenoSubset}}
}

\examples{
# load SNP genotypes
data(hla.snp.geno, package="HIBAG")
summary(hla.snp.geno)
# A/C A/G C/T G/T 
# 136 655 632 141 

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
