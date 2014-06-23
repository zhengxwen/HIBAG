\name{hlaReport}
\alias{hlaReport}
\title{
	Format a report
}
\description{
	Create a report for evaluating prediction accuracies.
}
\usage{
hlaReport(object, export.fn="", type=c("txt", "tex", "html"),
	header=TRUE)
}
\arguments{
	\item{object}{an object returned by \code{\link{hlaCompareAllele}}}
	\item{export.fn}{a file name for output, or "" for \code{stdout}}
	\item{type}{\code{"txt"} -- tab-delimited text format;
		\code{"tex"} -- tex format using the 'longtable' package;
		\code{"html"} -- html file}
	\item{header}{if \code{TRUE}, output the header of text file associated
		corresponding format}
}
\value{
	None.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaCompareAllele}}
}

\examples{
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

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
	samp.sel = match(hlatab$training$value$sample.id, HapMap_CEU_Geno$sample.id))
test.geno <- hlaGenoSubset(HapMap_CEU_Geno,
	samp.sel=match(hlatab$validation$value$sample.id, HapMap_CEU_Geno$sample.id))

# train a HIBAG model
set.seed(100)
# please use "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=4,
	verbose.detail=TRUE)
summary(model)

# validation
pred <- predict(model, test.geno)
# compare
(comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model,
	call.threshold=0))


# report
hlaReport(comp, type="txt")

hlaReport(comp, type="tex")

hlaReport(comp, type="html")
}

\keyword{genetics}
