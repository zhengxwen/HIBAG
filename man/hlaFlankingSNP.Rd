\name{hlaFlankingSNP}
\alias{hlaFlankingSNP}
\title{
    SNP IDs in Flanking Region
}
\description{
    To select SNPs in the flanking region of a specified HLA locus.
}
\usage{
hlaFlankingSNP(snp.id, position, hla.id, flank.bp=500*1000, assembly="auto",
    pos.start=NA_integer_, pos.end=NA_integer_)
}
\arguments{
    \item{snp.id}{a vector of SNP IDs}
    \item{position}{a vector of positions}
    \item{hla.id}{the name of HLA locus, or "any" for other genes}
    \item{flank.bp}{the size of flanking region on each side in basepair}
    \item{assembly}{the human genome reference: "hg18", "hg19" (default),
        "hg38"; "auto" refers to "hg19"; "auto-silent" refers to "hg19" without
        any warning}
    \item{pos.start}{the start position of the gene}
    \item{pos.end}{the end position of the gene}
}
\details{
    \code{hla.id} is "A", "B", "C", "DRB1", "DRB5", "DQA1", "DQB1", "DPB1" or
        "any".
}
\value{
    Return selected SNP IDs from \code{snp.id}.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaGenoSubset}}, \code{\link{hlaLociInfo}}
}

\examples{
# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
    H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
    H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
    locus=hla.id, assembly="hg19")

# training genotypes
region <- 500   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
    hla.id, region*1000, assembly="hg19")

train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
    snp.sel  = match(snpid, HapMap_CEU_Geno$snp.id),
    samp.sel = match(hla$training$value$sample.id, HapMap_CEU_Geno$sample.id))
summary(train.geno)


## customize positions
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
    "any", 500*1000, pos.start=29457816, pos.end=30450205)
}

\keyword{SNP}
\keyword{genetics}
