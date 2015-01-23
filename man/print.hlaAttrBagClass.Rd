\name{print.hlaAttrBagClass}
\alias{print.hlaAttrBagClass}
\alias{print.hlaAttrBagObj}
\alias{summary.hlaAttrBagClass}
\alias{summary.hlaAttrBagObj}
\title{
    Summarize a ``hlaAttrBagClass'' or ``hlaAttrBagObj'' object.
}
\description{
    Summarize an object of \code{\link{hlaAttrBagClass}} or
\code{\link{hlaAttrBagObj}}.
}
\usage{
\method{print}{hlaAttrBagClass}(x, ...)
\method{print}{hlaAttrBagObj}(x, ...)
\method{summary}{hlaAttrBagClass}(object, show=TRUE, ...)
\method{summary}{hlaAttrBagObj}(object, show=TRUE, ...)
}
\arguments{
    \item{x}{an object of \code{\link{hlaAttrBagClass}} or
        \code{\link{hlaAttrBagObj}}}
    \item{object}{an object of \code{\link{hlaAttrBagClass}} or
        \code{\link{hlaAttrBagObj}}}
    \item{show}{if TRUE, show information}
    \item{...}{further arguments passed to or from other methods}
}
\value{
    \code{print} returns \code{NULL}.

    \code{summary.hlaAttrBagClass} and \code{summary.hlaAttrBagObj} return
a \code{list}:
    \item{num.classifier}{the total number of classifiers}
    \item{num.snp}{the total number of SNPs}
    \item{snp.id}{SNP IDs}
    \item{snp.position}{SNP position in basepair}
    \item{snp.hist}{the number of classifier for each SNP, and it could be
        used for SNP importance}
    \item{info}{a \code{data.frame} for the average number of SNPs
        (\code{num.snp}), haplotypes (\code{num.haplo}), out-of-bag accuracies
        (\code{accuracy}) among all classifiers:
        mean, standard deviation, min, max}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{plot.hlaAttrBagClass}}, \code{\link{plot.hlaAttrBagObj}}
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
print(model)
}

\keyword{HLA}
\keyword{genetics}
