\name{hlaGeno2PED}
\alias{hlaGeno2PED}
\title{
	Convert to PLINK PED format
}
\description{
	Convert an object of \code{\link{hlaSNPGenoClass}} to a file of
PLINK PED format.
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
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaBED2Geno}}
}

\examples{
# load SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
	H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
	H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
	max.resolution=4, locus=hla.id, assembly="hg19")

# training genotypes
region <- 500   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
	hla.id, region*1000, assembly="hg19")

train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
	snp.sel = match(snpid, HapMap_CEU_Geno$snp.id))

hlaGeno2PED(train.geno, "test")
}

\keyword{SNP}
\keyword{genetics}
