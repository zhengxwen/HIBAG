\name{HIBAG-package}
\alias{HIBAG-package}
\alias{HIBAG}
\docType{package}
\title{HIBAG -- HLA Genotype Imputation with Attribute Bagging}
\description{
	To impute HLA types from unphased SNP data using an attribute bagging method.
}
\details{
	\tabular{ll}{
		Package: \tab HIBAG\cr
		Type:    \tab Package\cr
		Version: \tab 0.9.1\cr
		Date:    \tab 2012-06-20\cr
		License: \tab GPL version 3\cr
	}
}
\references{
	\url{http://students.washington.edu/zhengx/HIBAG/}

	X. Zheng, J. Shen, C. Cox, J. Wakefield, M. Ehm, M. Nelson, B. Weir.
	HIBAG -- HLA Genotype Imputation with Attribute Bagging.
	submitted.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}

\examples{
# load HLA types and SNP genotypes
data(hla.type.table, package="HIBAG")
data(hla.snp.geno, package="HIBAG")

head(hla.type.table)
dim(hla.type.table)  # 60 13

summary(hla.snp.geno)


# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
	hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)
summary(hla)


# plot linkage disequilibrium
ymax <- 0.16
plot(NaN, NaN, xlab="SNP Position (in KB)", ylab="Linkage Disequilibrium",
	xlim=range(hla.snp.geno$snp.position)/1000, ylim=c(0, ymax),
	main="Major Histocompatibility Complex")
# for-loop
for (i in 1:3)
{
	hla.id <- c("A", "C", "DQA1")[i]
	# make a "hlaAlleleClass" object
	hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
		hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)
	# linkage disequilibrium between HLA locus and SNP markers
	ld <- hlaGenoLD(hla, hla.snp.geno)
	# draw
	points(hla.snp.geno$snp.position/1000, ld, pch="*", col=i)
	x <- (hla$pos.start/1000 + hla$pos.end/1000)/2
	abline(v=x, col=i, lty=3)
	text(x, ymax-(i-1)*0.005, hla.id, adj=1, col=i)
}



######################################################################

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
	hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)

# divide HLA types randomly
set.seed(100)
hlatab <- hlaSplitAllele(hla, train.prop=0.5)
names(hlatab)
# "training"   "validation"
summary(hlatab$training)
summary(hlatab$validation)

# SNP predictors within the flanking region on each side
region <- 500   # kb
snpid <- hlaFlankingSNP(hla.snp.geno$snp.id, hla.snp.geno$snp.position, hla.id, region*1000)
length(snpid)  # 275

# training and validation genotypes
train.geno <- hlaGenoSubset(hla.snp.geno, snp.sel=match(snpid, hla.snp.geno$snp.id),
	samp.sel=match(hlatab$training$value$sample.id, hla.snp.geno$sample.id))
test.geno <- hlaGenoSubset(hla.snp.geno,
	samp.sel=match(hlatab$validation$value$sample.id, hla.snp.geno$sample.id))

# train an attribute bagging model
set.seed(100)
# please try "nclassifier=25" or "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=2, verbose.detail=TRUE)
summary(model)

# validation
pred <- predict(model, test.geno)
# compare
comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model$hla.allele,
	call.threshold=0.5)
comp
# total.num.ind crt.num.ind crt.num.haplo   acc.ind acc.haplo call.threshold n.call call.rate
#            26          20            42 0.9090909 0.9545455            0.5     22 0.8461538


# save the parameter file
mobj <- hlaModelToObj(model)
save(mobj, file="HIBAG_model.RData")
save(test.geno, file="testgeno.RData")
save(hlatab, file="HLASplit.RData")

# Clear Workspace
hlaClose(model)  # release all resources of model
rm(list = ls())


######################################################################

# NOW, load a HIBAG model from the parameter file
mobj <- get(load("HIBAG_model.RData"))
model <- hlaModelfromObj(mobj)

# validation
test.geno <- get(load("testgeno.RData"))
hlatab <- get(load("HLASplit.RData"))

pred <- predict(model, test.geno)
# compare
comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model$hla.allele,
	call.threshold=0.5)
comp
# total.num.ind crt.num.ind crt.num.haplo   acc.ind acc.haplo call.threshold n.call call.rate
#            26          20            42 0.9090909 0.9545455            0.5     22 0.8461538


#########################################################################
# import a PLINK BED file
#
bed.fn <- system.file("extdata", "hapmap.ceu.bed", package="HIBAG")
fam.fn <- system.file("extdata", "hapmap.ceu.fam", package="HIBAG")
bim.fn <- system.file("extdata", "hapmap.ceu.bim", package="HIBAG")
hapmap.ceu <- hlaBED2Geno(bed.fn, fam.fn, bim.fn)


#########################################################################
# predict
#
pred <- predict(model, hapmap.ceu, type="response")
pred$value
#    sample.id allele1 allele2      prob
# 1    NA10859   01:01   03:01 0.5000786
# 2    NA11882   01:01   29:02 0.4999865
# ...
# 89   NA12892   02:01   11:01 0.9999992
# 90   NA12891   01:01   24:02 0.8888880
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
