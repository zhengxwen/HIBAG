\name{plot.hlaAttrBagObj}
\alias{plot.hlaAttrBagObj}
\alias{plot.hlaAttrBagClass}
\title{
    Plot a HIBAG model
}
\description{
    To show a scatterplot of the numbers of individual classifiers and SNP
positions.
}
\usage{
\method{plot}{hlaAttrBagObj}(x, snp.col="gray33", snp.pch=1, snp.sz=1,
    locus.col="blue", locus.lty=1L, locus.lty2=2L, addplot=NULL,
    assembly="auto", ...)
\method{plot}{hlaAttrBagClass}(x, ...)
}
\arguments{
    \item{x}{an object of \code{\link{hlaAttrBagObj}}}
    \item{snp.col}{the color of SNP uses}
    \item{snp.pch}{the point type of SNP uses}
    \item{snp.sz}{the point size of SNP uses}
    \item{locus.col}{the color of text and line for HLA locus}
    \item{locus.lty}{the type of line for the bounds of HLA locus}
    \item{locus.lty2}{the type of line for HLA locus}
    \item{addplot}{NULL for creating a plot, or a ggplot object to be appended}
    \item{assembly}{the human genome reference: "hg18", "hg19" (default),
        "hg38"; "auto" refers to "hg19"; "auto-silent" refers to "hg19" without
        any warning}
    \item{...}{further arguments passed to or from other methods}
}
\value{None}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{print.hlaAttrBagObj}}, \code{\link{summary.hlaAttrBagObj}}
}

\examples{
# make a "hlaAlleleClass" object
hla.id <- "C"
hla <- hlaAllele(HLA_Type_Table$sample.id,
    H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
    H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
    locus=hla.id, assembly="hg19")

# training genotypes
region <- 100   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
    hla.id, region*1000, assembly="hg19")
train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
    snp.sel = match(snpid, HapMap_CEU_Geno$snp.id))

# train a HIBAG model
set.seed(1000)
# please use "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hla, train.geno, nclassifier=2, verbose.detail=TRUE)
plot(model)
}

\keyword{HLA}
\keyword{genetics}
