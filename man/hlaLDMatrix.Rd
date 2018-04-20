\name{hlaLDMatrix}
\alias{hlaLDMatrix}
\title{
    Composite Linkage Disequilibrium in a Region
}
\description{
    To calculate composite linkage disequilibrium (r2) among SNPs within a
region.
}
\usage{
hlaLDMatrix(geno, loci=NULL, maf=0.01, assembly="auto", draw=TRUE,
    verbose=TRUE)
}
\arguments{
    \item{geno}{an object of \code{\link{hlaSNPGenoClass}}}
    \item{maf}{MAF filter \code{>= maf}}
    \item{loci}{NULL or a character vector, e.g., "A", "B"}
    \item{assembly}{the human genome reference: "hg18", "hg19" (default),
        "hg38"; "auto" refers to "hg19"; "auto-silent" refers to "hg19" without
        any warning}
    \item{draw}{if TRUE, return a ggplot2 object}
    \item{verbose}{if TRUE, show information}
}
\value{
    Return a ggplot2 object if \code{draw=TRUE} or a matrix correlation.
}
\references{
    Weir BS, Cockerham CC:
Complete characterization of disequilibrium at two loci;
in Feldman MW (ed): Mathematical Evolutionary Theory.
Princeton, NJ: Princeton University Press, 1989.
}
\author{Xiuwen Zheng}

\examples{
region <- 500*1000   # basepair
geno <- hlaGenoSubsetFlank(HapMap_CEU_Geno, "A", region)
summary(geno)

hlaLDMatrix(geno, "A")
}

\keyword{SNP}
\keyword{genetics}
