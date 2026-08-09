\name{hlaGDS2Geno}
\alias{hlaGDS2Geno}
\title{
    Import genotypes from a GDS file
}
\description{
    To convert a SNPRelate or SeqArray GDS file to an object of \code{hlaSNPGenoClass}.
}
\usage{
hlaGDS2Geno(gds.fn, rm.invalid.allele=FALSE, import.chr="xMHC", position=NULL,
    assembly="auto", verbose=TRUE)
}
\arguments{
    \item{gds.fn}{a file name for the GDS file defined in the SNPRelate or SeqArray
        package}
    \item{rm.invalid.allele}{if TRUE, remove SNPs with non-standard alleles
        (except A,G,C,T); FALSE is suggested (by default), since the A allele
        is set to be "0" at a multi-allelic site in a SeqArray GDS file (e.g.,
        "0/G"), and such a SNP would be removed if TRUE}
    \item{import.chr}{the chromosome, "1" .. "22", "X", "Y", "XY", "MT",
        "xMHC", or "", where "xMHC" implies the extended MHC on chromosome 6,
        and "" for all SNPs}
    \item{position}{NULL for no restriction on the SNP positions (by default),
        or an integer vector of the base positions of the SNPs to be imported;
        it is used together with \code{import.chr}, i.e., only the SNPs in
        \code{import.chr} and at the specified positions are imported}
    \item{assembly}{the human genome reference: "hg18", "hg19" (default),
        "hg38"; "auto" refers to "hg19"; "auto-silent" refers to "hg19" without
        any warning}
    \item{verbose}{if TRUE, show information}
}
\value{
    Return an object of \code{\link{hlaSNPGenoClass}}.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaGeno2PED}}, \code{\link{hlaBED2Geno}}
}

\examples{
# Import a SNP GDS file
fn <- system.file("extdata", "HapMap_CEU_Chr6.gds", package="HIBAG")

geno <- hlaGDS2Geno(fn, assembly="hg18")

summary(geno)

# import the SNPs at the specified positions only
geno2 <- hlaGDS2Geno(fn, position=geno$snp.position[1:100], assembly="hg18")

summary(geno2)
}

\keyword{SNP}
\keyword{genetics}
