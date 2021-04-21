\name{hlaAlleleToVCF}
\alias{hlaAlleleToVCF}
\title{
    Convert HLA alleles to VCF
}
\description{
    To convert the HLA allele data to a VCF file.
}
\usage{
hlaAlleleToVCF(hla, outfn, DS=TRUE, verbose=TRUE)
}
\arguments{
    \item{hla}{the training HLA types, an object of \code{\link{hlaAlleleClass}}}
    \item{outfn}{a VCF file name or a \code{connection}}
    \item{DS}{if TRUE, output dosages in the DS field}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    None.
}
\references{
    Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS;
    HIBAG -- HLA Genotype Imputation with Attribute Bagging.
    Pharmacogenomics Journal. doi: 10.1038/tpj.2013.18.
    \url{https://www.nature.com/articles/tpj201318}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaAttrBagging}}, \code{\link{hlaAllele}}
}

\examples{
# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
    H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
    H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
    locus=hla.id, assembly="hg19")

# SNP predictors within the flanking region on each side
region <- 500   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
    hla.id, region*1000, assembly="hg19")
length(snpid)  # 275

# train a HIBAG model
set.seed(100)
# please use "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hla, HapMap_CEU_Geno, nclassifier=2)
summary(model)

# validation
pred <- hlaPredict(model, HapMap_CEU_Geno)
summary(pred)


# output to standard output with dosages
hlaAlleleToVCF(hlaAlleleSubset(pred, 1:4), stdout())
}

\keyword{HLA}
\keyword{SNP}
\keyword{genetics}
