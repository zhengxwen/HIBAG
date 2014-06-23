\name{hlaGenoCombine}
\alias{hlaGenoCombine}
\title{
	Combine two genotype dataset into one
}
\description{
	To combine two genotype datasets into one dataset.
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
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaMakeSNPGeno}}, \code{\link{hlaMakeSNPHaplo}},
	\code{\link{hlaHaploSubset}}, \code{\link{hlaGenoSubset}}
}

\examples{
# load SNP genotypes
data(hla.snp.geno, package="HIBAG")

# import a PLINK BED file
bed.fn <- system.file("extdata", "hapmap.ceu.bed", package="HIBAG")
fam.fn <- system.file("extdata", "hapmap.ceu.fam", package="HIBAG")
bim.fn <- system.file("extdata", "hapmap.ceu.bim", package="HIBAG")
hapmap.ceu <- hlaBED2Geno(bed.fn, fam.fn, bim.fn)

# combine two datasets together
geno <- hlaGenoCombine(hla.snp.geno, hapmap.ceu)
summary(geno)
#	390 samples X 1564 SNPs
#	SNPs range from 25877003bp to 33529555bp
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
