\name{hlaModelFiles}
\alias{hlaModelFiles}
\title{
	Load a model object from files
}
\description{
	To load HIBAG models from a list of files, and merge all together.
}
\usage{
hlaModelFiles(fn.list, action.missingfile=c("ignore", "stop"), verbose=TRUE)
}
\arguments{
	\item{fn.list}{a vector of file names}
	\item{action.missingfile}{"ignore", ignore the missing files, by default;
		"stop", stop if missing}
	\item{verbose}{if TRUE, show information}
}
\value{
	Return \code{\link{hlaAttrBagObj}}.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{hlaAllele}}, \code{\link{hlaFlankingSNP}},
\code{\link{hlaModelToObj}}
}

\examples{
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "C"
hla <- hlaAllele(HLA_Type_Table$sample.id, HLA_Type_Table[, paste(hla.id, ".1", sep="")],
	HLA_Type_Table[, paste(hla.id, ".2", sep="")], locus=hla.id, assembly="hg19")

# training genotypes
region <- 100   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
	hla.id, region*1000, assembly="hg19")
train.geno <- hlaGenoSubset(HapMap_CEU_Geno, snp.sel=match(snpid, HapMap_CEU_Geno$snp.id),
	samp.sel=match(hla$value$sample.id, HapMap_CEU_Geno$sample.id))

#
# train a HIBAG model
#
set.seed(1000)

# please use "nclassifier=100" when you use HIBAG for real data
model1 <- hlaAttrBagging(hla, train.geno, nclassifier=1, verbose.detail=TRUE)
mobj1 <- hlaModelToObj(model1)
save(mobj1, file="test-model1.RData")

model2 <- hlaAttrBagging(hla, train.geno, nclassifier=1, verbose.detail=TRUE)
mobj2 <- hlaModelToObj(model2)
save(mobj2, file="test-model2.RData")

model3 <- hlaAttrBagging(hla, train.geno, nclassifier=1, verbose.detail=TRUE)
mobj3 <- hlaModelToObj(model3)
save(mobj3, file="test-model3.RData")

# load all of mobj1, mobj2 and mobj3
mobj <- hlaModelFiles(c("test-model1.RData", "test-model2.RData", "test-model3.RData"))
summary(mobj)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
