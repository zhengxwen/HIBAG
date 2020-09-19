\name{hlaReportPlot}
\alias{hlaReportPlot}
\title{
    Format a report with figures
}
\description{
    Create figures for evaluating prediction accuracies.
}
\usage{
hlaReportPlot(PredHLA=NULL, TrueHLA=NULL, model=NULL,
    fig=c("matching", "call.rate", "call.threshold"), match.threshold=NaN,
    log_scale=TRUE)
}
\arguments{
    \item{PredHLA}{NULL, an object of \code{\link{hlaAlleleClass}}, the
        predicted HLA types}
    \item{TrueHLA}{NULL, an object of \code{\link{hlaAlleleClass}}, the true
        HLA types}
    \item{model}{NULL, or a model of \code{\link{hlaAttrBagClass}}}
    \item{fig}{"matching": violin plot for matching measurements;
        "call.rate": relationship between accuracy and call rate;
        "call.threshold": relationship between accuracy and call threshold}
    \item{match.threshold}{the threshold for matching proportion}
    \item{log_scale}{if TRUE, use log scale for matching violin plot}
}
\value{
    Return a ggplot2 object.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaReport}}
}

\examples{
# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
    H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
    H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
    locus=hla.id, assembly="hg19")

# divide HLA types randomly
set.seed(100)
hlatab <- hlaSplitAllele(hla, train.prop=0.5)
names(hlatab)
# "training"   "validation"
summary(hlatab$training)
summary(hlatab$validation)

# SNP predictors within the flanking region on each side
region <- 500   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
    hla.id, region*1000, assembly="hg19")
length(snpid)  # 275

# training and validation genotypes
train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
    snp.sel = match(snpid, HapMap_CEU_Geno$snp.id),
    samp.sel = match(hlatab$training$value$sample.id,
    HapMap_CEU_Geno$sample.id))
test.geno <- hlaGenoSubset(HapMap_CEU_Geno,
    samp.sel=match(hlatab$validation$value$sample.id,
    HapMap_CEU_Geno$sample.id))

# train a HIBAG model
set.seed(100)
# please use "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=4,
    verbose.detail=TRUE)
summary(model)

# validation
pred <- hlaPredict(model, test.geno)


# visualize
hlaReportPlot(pred, fig="matching")

hlaReportPlot(model=model, fig="matching")

hlaReportPlot(pred, model=model, fig="matching")

hlaReportPlot(pred, hlatab$validation, fig="call.rate")

hlaReportPlot(pred, hlatab$validation, fig="call.threshold")
}

\keyword{genetics}
