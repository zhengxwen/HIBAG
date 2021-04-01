\name{hlaCheckSNPs}
\alias{hlaCheckSNPs}
\title{
    Check the SNP predictors in a HIBAG model
}
\description{
    Check the SNP predictors in a HIBAG model, by calculating the overlapping
between the model and SNP genotypes.
}
\usage{
hlaCheckSNPs(model, object,
    match.type=c("Position", "Pos+Allele", "RefSNP+Position", "RefSNP"), verbose=TRUE)
}
\arguments{
    \item{model}{an object of \code{\link{hlaAttrBagClass}}, or an object of
        \code{\link{hlaAttrBagObj}}}
    \item{object}{a genotype object of \code{\link{hlaSNPGenoClass}}, or a
        character vector like c("rs2523442", "rs9257863", ...)}
    \item{match.type}{\code{"RefSNP+Position"} (by default) -- using both of
        RefSNP IDs and positions; \code{"RefSNP"} -- using RefSNP IDs only;
        \code{"Position"} -- using positions only}
    \item{verbose}{if TRUE, show information}
}
\value{
    Return a \code{data.frame} for individual classifiers:
    \item{NumOfValidSNP}{the number of non-missing SNPs in an individual
        classifier}
    \item{NumOfSNP}{the number of SNP predictors in an individual classifier}
    \item{fraction}{NumOfValidSNP / NumOfSNP}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaAttrBagging}}, \code{\link{predict.hlaAttrBagClass}}
}

\examples{
# make a "hlaAlleleClass" object
hla.id <- "DQB1"
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
model <- hlaAttrBagging(hla, train.geno, nclassifier=2)
print(model)


hlaCheckSNPs(model, train.geno)

# close the HIBAG model explicitly
hlaClose(model)
}

\keyword{SNP}
\keyword{genetics}
