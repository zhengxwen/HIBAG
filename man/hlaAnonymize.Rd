\name{hlaAnonymize}
\alias{hlaAnonymize}
\title{
	Make the IDs of training samples anonymous
}
\description{
	Changes the \code{sample.id} of a model object to ``unmeaningful'' index.
}
\usage{
hlaAnonymize(modelobj)
}
\arguments{
	\item{modelobj}{an object of \code{\link{hlaAttrBagObj}}}
}
\value{
	Returns a new object of \code{\link{hlaAttrBagObj}}.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaModelFiles}}, \code{\link{hlaModelfromObj}},
	\code{\link{hlaModelToObj}}
}

\examples{
# load HLA types and SNP genotypes
data(hla.type.table, package="HIBAG")
data(hla.snp.geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
	hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)

# training genotypes
region <- 250   # kb
snpid <- hlaFlankingSNP(hla.snp.geno$snp.id, hla.snp.geno$snp.position, hla.id, region*1000)
train.geno <- hlaGenoSubset(hla.snp.geno, snp.sel=match(snpid, hla.snp.geno$snp.id),
	samp.sel=match(hla$value$sample.id, hla.snp.geno$sample.id))

#
# train an attribute bagging model
#
set.seed(1000)

# please try "nclassifier=25" or "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hla, train.geno, nclassifier=2, verbose.detail=TRUE)
mobj <- hlaModelToObj(model)
summary(mobj)


# Make the IDs of training samples anonymous
head(mobj$sample.id)
# [1] "NA11882" "NA11881" "NA11993" "NA11992" "NA11995" "NA11994"

head(hlaAnonymize(mobj)$sample.id)
# [1] 1 2 3 4 5 6
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
