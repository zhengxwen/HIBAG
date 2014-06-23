\name{hlaGeno2PED}
\alias{hlaGeno2PED}
\title{
	Convert to PLINK PED format
}
\description{
	Convert an object of \code{\link{hlaSNPGenoClass}} to a file of PLINK PED format.
}
\usage{
hlaGeno2PED(geno, out.fn)
}
\arguments{
	\item{geno}{a genotype object of \code{\link{hlaSNPGenoClass}}}
	\item{out.fn}{the file name of output ped file}
}
\details{
	Two files ".map" and ".ped" are created.
}
\value{
	None.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaBED2Geno}}
}

\examples{
# load SNP genotypes
data(hla.type.table, package="HIBAG")
data(hla.snp.geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
	hla.type.table[, paste(hla.id, ".2", sep="")], max.resolution=4, locus=hla.id)

# training genotypes
region <- 500   # kb
snpid <- hlaFlankingSNP(hla.snp.geno$snp.id, hla.snp.geno$snp.position, hla.id, region*1000)

train.geno <- hlaGenoSubset(hla.snp.geno, snp.sel=match(snpid, hla.snp.geno$snp.id))

hlaGeno2PED(train.geno, "test")
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
