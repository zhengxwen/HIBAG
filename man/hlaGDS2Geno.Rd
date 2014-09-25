\name{hlaGDS2Geno}
\alias{hlaGDS2Geno}
\title{
	Convert from SNP GDS format
}
\description{
	To convert a SNP GDS file to an object of \code{\link{hlaSNPGenoClass}}.
}
\usage{
hlaGDS2Geno(gds.fn, rm.invalid.allele=FALSE,
	import.chr="xMHC", assembly=c("auto", "hg18", "hg19", "unknown"),
	verbose=TRUE)
}
\arguments{
	\item{gds.fn}{the SNP GDS file used by the \code{SNPRelate} package}
	\item{rm.invalid.allele}{if TRUE, remove SNPs with invalid alleles}
	\item{import.chr}{the chromosome, "1" .. "22", "X", "Y", "XY", "MT",
		"xMHC", or "", where "xMHC" implies the extended MHC on chromosome 6,
		and "" for all SNPs}
	\item{assembly}{the human genome reference: "hg19" (default), "hg18",
		"auto" refers to "hg19"}
	\item{verbose}{if TRUE, show information}
}
\value{
	Return an object of \code{\link{hlaSNPGenoClass}}.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaGeno2PED}}, \code{\link{hlaBED2Geno}}
}

\keyword{SNP}
\keyword{genetics}
