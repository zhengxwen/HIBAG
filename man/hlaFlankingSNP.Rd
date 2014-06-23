\name{hlaFlankingSNP}
\alias{hlaFlankingSNP}
\title{
	SNP IDs in Flanking Region
}
\description{
	To select SNPs in the flanking region of a specified HLA locus.
}
\usage{
hlaFlankingSNP(snp.id, position, hla.id, flank.bp=500*1000)
}
\arguments{
	\item{snp.id}{a vector of SNP IDs}
	\item{position}{a vector of positions}
	\item{hla.id}{the name of HLA locus}
	\item{flank.bp}{the size of flanking region on each side in basepair}
}
\details{
	\code{hla.id} is "A", "B", "C", "DRB1", "DRB5", "DQA1", "DQB1", "DPB1" or "any".
}
\value{
	Return allele frequecies.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{hlaClose}}
}

\examples{
# load HLA types and SNP genotypes
data(hla.type.table, package="HIBAG")
data(hla.snp.geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
	hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)

# training genotypes
region <- 500   # kb
snpid <- hlaFlankingSNP(hla.snp.geno$snp.id, hla.snp.geno$snp.position, hla.id, region*1000)

train.geno <- hlaGenoSubset(hla.snp.geno, snp.sel=match(snpid, hla.snp.geno$snp.id),
	samp.sel=match(hla$training$value$sample.id, hla.snp.geno$sample.id))
summary(train.geno)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
