\name{hlaBED2Geno}
\alias{hlaBED2Geno}
\title{
	Convert from PLINK BED format
}
\description{
	To convert a PLINK BED file to an object of \code{\link{hlaSNPGenoClass}}.
}
\usage{
hlaBED2Geno(bed.fn, fam.fn, bim.fn, rm.invalid.allele=FALSE, import.chr="xMHC", verbose=TRUE)
}
\arguments{
	\item{bed.fn}{binary file, genotype information}
	\item{fam.fn}{family, individual information, etc}
	\item{bim.fn}{extended MAP file: two extra cols = allele names}
	\item{rm.invalid.allele}{if TRUE, remove SNPs with invalid alleles}
	\item{import.chr}{the chromosome, "1" .. "22", "X", "Y", "XY", "MT" or "xMHC",
		where "xMHC" implies the xMHC on chromosome 6}
	\item{verbose}{if TRUE, show information}
}
\value{
	Return an object of \code{\link{hlaSNPGenoClass}}.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaGeno2PED}}
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


#########################################################################
# train an attribute bagging model
#
set.seed(1000)
# please try "nclassifier=25" or "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hla, train.geno, nclassifier=2, verbose.detail=TRUE)
summary(model)


#########################################################################
# import a PLINK BED file
#
bed.fn <- system.file("extdata", "hapmap.ceu.bed", package="HIBAG")
fam.fn <- system.file("extdata", "hapmap.ceu.fam", package="HIBAG")
bim.fn <- system.file("extdata", "hapmap.ceu.bim", package="HIBAG")
hapmap.ceu <- hlaBED2Geno(bed.fn, fam.fn, bim.fn)


#########################################################################
# validation
#
pred <- predict(model, hapmap.ceu)
pred$value
#    sample.id     types allele1 allele2      prob
# 1    NA10859 0101/0301    0101    0301 0.9999993
# 2    NA11882 0101/2902    0101    2902 0.9999993
# 3    NA11881 0301/2601    0301    2601 0.9999791
# 4    NA10860 0201/0201    0201    0201 0.1596430
# 5    NA11993 2501/2902    2501    2902 0.8465468
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
