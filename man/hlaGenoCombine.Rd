\name{hlaGenoCombine}
\alias{hlaGenoCombine}
\title{
	Combine two genotypic data sets into one
}
\description{
	To combine two genotypic data sets into one dataset.
}
\usage{
hlaGenoCombine(geno1, geno2,
	match.type=c("RefSNP+Position", "RefSNP", "Position"),
	allele.check=TRUE, same.strand=FALSE, verbose=TRUE)
}
\arguments{
	\item{geno1}{the first genotype object of \code{\link{hlaSNPGenoClass}}}
	\item{geno2}{the second genotype object of \code{\link{hlaSNPGenoClass}}}
	\item{match.type}{\code{"RefSNP+Position"} (by default) -- using both of
		RefSNP IDs and positions; \code{"RefSNP"} -- using RefSNP IDs only;
		\code{"Position"} -- using positions only}
	\item{allele.check}{if \code{TRUE}, call \code{\link{hlaGenoSwitchStrand}}
		to check and then switch allele pairs if needed}
	\item{same.strand}{\code{TRUE} assuming alleles are on the same strand
		(e.g., forward strand); otherwise, \code{FALSE} not assuming whether
		on the same strand or not}
	\item{verbose}{show information, if TRUE}
}
\details{
	The function merges two SNP dataset \code{geno1} and \code{geno2}, and
returns a SNP dataset consisting of the SNP intersect between \code{geno1} and
\code{geno2}, and having the same SNP information (allele and position) as
\code{geno1}.
}
\value{
	An object of \code{\link{hlaSNPGenoClass}}.
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
}

\keyword{SNP}
\keyword{genetics}
