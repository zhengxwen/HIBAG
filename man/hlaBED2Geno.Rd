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
	\code{\link{hlaGeno2PED}}, \code{\link{hlaAttrBagging}}
}

\examples{
#########################################################################
# Import a PLINK BED file
#
bed.fn <- system.file("extdata", "hapmap.ceu.bed", package="HIBAG")
fam.fn <- system.file("extdata", "hapmap.ceu.fam", package="HIBAG")
bim.fn <- system.file("extdata", "hapmap.ceu.bim", package="HIBAG")
hapmap.ceu <- hlaBED2Geno(bed.fn, fam.fn, bim.fn)

summary(hapmap.ceu)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
