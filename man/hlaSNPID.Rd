\name{hlaSNPID}
\alias{hlaSNPID}
\title{
    Get SNP IDs and positions
}
\description{
    Get the information of SNP ID with or without position.
}
\usage{
hlaSNPID(obj, type=c("Position", "Pos+Allele", "RefSNP+Position", "RefSNP"))
}
\arguments{
    \item{obj}{a genotypic object of \code{\link{hlaSNPGenoClass}},
        a model object of \code{\link{hlaAttrBagClass}} or
        a model object of \code{\link{hlaAttrBagObj}}}
    \item{type}{\code{"RefSNP+Position"} (by default), \code{"RefSNP"} or
        \code{"Position"}}
}
\value{
    If \code{type = "RefSNP+Position"},
        return \code{paste(obj$snp.id, obj$snp.position, sep="-")};
    if \code{type = "RefSNP"},
        return \code{obj$snp.id};
    if \code{type = "Position"},
        return \code{obj$snp.position};
    if \code{type = "Pos+Allele"},
        return \code{paste(obj$snp.position, obj$snp.allele, sep="-")}.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaGenoSwitchStrand}}, \code{\link{hlaGenoCombine}}
}

\examples{
x <- hlaSNPID(HapMap_CEU_Geno)
head(x)

x <- hlaSNPID(HapMap_CEU_Geno, "RefSNP")
head(x)

x <- hlaSNPID(HapMap_CEU_Geno, "Position")
head(x)
}

\keyword{SNP}
\keyword{genetics}
