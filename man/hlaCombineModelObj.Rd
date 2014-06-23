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
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{hlaParallelAttrBagging}}
}

\examples{
# load HLA types and SNP genotypes
data(hla.type.table, package="HIBAG")
data(hla.snp.geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
	hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)

# SNP predictors within the flanking region on each side
region <- 500   # kb
snpid <- hlaFlankingSNP(hla.snp.geno$snp.id, hla.snp.geno$snp.position, hla.id, region*1000)
length(snpid)  # 275

# training genotypes
train.geno <- hlaGenoSubset(hla.snp.geno, snp.sel=match(snpid, hla.snp.geno$snp.id))

# train an attribute bagging model
set.seed(100)
m1 <- hlaAttrBagging(hla, train.geno, nclassifier=1)
m2 <- hlaAttrBagging(hla, train.geno, nclassifier=1)

m1.obj <- hlaModelToObj(m1)
m2.obj <- hlaModelToObj(m2)

m.obj <- hlaCombineModelObj(m1.obj, m2.obj)
summary(m.obj)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
