\name{hlaModelFiles}
\alias{hlaModelFiles}
\title{
	Load a model object from files
}
\description{
	To load HIBAG models from a list of files, and merge all together.
}
\usage{
hlaModelFiles(fn.list, action.missingfile=c("ignore", "stop"))
}
\arguments{
	\item{fn.list}{a vector of file names}
	\item{action.missingfile}{"ignore", ignore the missing files, by default;
		"stop", stop if missing}
}
\value{
	Return \code{\link{hlaAttrBagObj}}.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{hlaAllele}}, \code{\link{hlaFlankingSNP}},
\code{\link{hlaModelToObj}}
}

\examples{
# load HLA types and SNP genotypes
data(hla.type.table, package="HIBAG")
data(hla.snp.geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "C"
hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
	hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)

# training genotypes
region <- 100   # kb
snpid <- hlaFlankingSNP(hla.snp.geno$snp.id, hla.snp.geno$snp.position, hla.id, region*1000)
train.geno <- hlaGenoSubset(hla.snp.geno, snp.sel=match(snpid, hla.snp.geno$snp.id),
	samp.sel=match(hla$value$sample.id, hla.snp.geno$sample.id))

#
# train an attribute bagging model
#
set.seed(1000)

# please try "nclassifier=25" or "nclassifier=100" when you use HIBAG for real data
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
