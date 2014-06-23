\name{hlaFlankingSNP}
\alias{hlaFlankingSNP}
\title{
	SNP IDs in Flanking Region
}
\description{
	To select SNPs in the flanking region of a specified HLA locus.
}
\usage{
hlaFlankingSNP(snp.id, position, hla.id, flank.bp=500*1000, assembly="auto")
}
\arguments{
	\item{snp.id}{a vector of SNP IDs}
	\item{position}{a vector of positions}
	\item{hla.id}{the name of HLA locus}
	\item{flank.bp}{the size of flanking region on each side in basepair}
	\item{assembly}{genome assembly: "hg19" (default), "hg18", "NCBI37" or "NCBI36";
		"auto" refers to "hg19"}
}
\details{
	\code{hla.id} is "A", "B", "C", "DRB1", "DRB5", "DQA1", "DQB1", "DPB1" or "any".
}
\value{
	Return allele frequecies.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{hlaClose}}
}

\examples{
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id, HLA_Type_Table[, paste(hla.id, ".1", sep="")],
	HLA_Type_Table[, paste(hla.id, ".2", sep="")], locus=hla.id, assembly="hg19")

# training genotypes
region <- 500   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
	hla.id, region*1000, assembly="hg19")

train.geno <- hlaGenoSubset(HapMap_CEU_Geno, snp.sel=match(snpid, HapMap_CEU_Geno$snp.id),
	samp.sel=match(hla$training$value$sample.id, HapMap_CEU_Geno$sample.id))
summary(train.geno)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
