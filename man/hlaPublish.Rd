\name{hlaPublish}
\alias{hlaPublish}
\title{
	Finalize a HIBAG model
}
\description{
	Finalize a HIBAG model by removing unused SNP predictors and adding
appendix information (platform, training set, authors, warning, etc)
}
\usage{
hlaPublish(mobj, platform=NULL, information=NULL, warning=NULL,
	rm.unused.snp=TRUE, anonymize=TRUE, verbose=TRUE)
}
\arguments{
	\item{mobj}{an object of \code{\link{hlaAttrBagObj}} or
		\code{\link{hlaAttrBagClass}}}
	\item{platform}{the text of platform information}
	\item{information}{the other information, like authors}
	\item{warning}{any warning message}
	\item{rm.unused.snp}{if \code{TRUE}, remove unused SNPs from the model}
	\item{anonymize}{if \code{TRUE}, remove sample IDs}
	\item{verbose}{if TRUE, show information}
}
\value{
	Returns a new object of \code{\link{hlaAttrBagObj}}.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaModelFromObj}}, \code{\link{hlaModelToObj}}
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

# training genotypes
region <- 250   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
	hla.id, region*1000, assembly="hg19")
train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
	snp.sel = match(snpid, HapMap_CEU_Geno$snp.id),
	samp.sel = match(hla$value$sample.id, HapMap_CEU_Geno$sample.id))


#
# train a HIBAG model
#
set.seed(1000)

# please use "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hla, train.geno, nclassifier=2, verbose.detail=TRUE)
summary(model)
length(model$snp.id)

mobj <- hlaPublish(model,
	platform = "Illumina 1M Duo",
	information = "Training set -- HapMap Phase II")
model2 <- hlaModelFromObj(mobj)
length(mobj$snp.id)
mobj$appendix
summary(mobj)

p1 <- predict(model, train.geno)
p2 <- predict(model2, train.geno)

# check
cbind(p1$value, p2$value)
}

\keyword{genetics}
