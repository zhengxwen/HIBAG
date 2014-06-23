\name{hlaCombineModelObj}
\alias{hlaCombineModelObj}
\title{
	Combine two HIBAG models together
}
\description{
	Merge two objects of \code{\link{hlaAttrBagObj}} together, which is useful
for building an ensemble model in parallel.
}
\usage{
hlaCombineModelObj(obj1, obj2)
}
\arguments{
	\item{obj1}{an object of \code{\link{hlaAttrBagObj}}}
	\item{obj2}{an object of \code{\link{hlaAttrBagObj}}}
}
\value{
	Return an object of \code{\link{hlaAttrBagObj}}.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{hlaModelFiles}}
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

# SNP predictors within the flanking region on each side
region <- 500   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
	hla.id, region*1000, assembly="hg19")
length(snpid)  # 275

# training genotypes
train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
	snp.sel = match(snpid, HapMap_CEU_Geno$snp.id))

# train a HIBAG model
set.seed(100)
m1 <- hlaAttrBagging(hla, train.geno, nclassifier=1)
m2 <- hlaAttrBagging(hla, train.geno, nclassifier=1)

m1.obj <- hlaModelToObj(m1)
m2.obj <- hlaModelToObj(m2)

m.obj <- hlaCombineModelObj(m1.obj, m2.obj)
summary(m.obj)
}

\keyword{HLA}
\keyword{genetics}
