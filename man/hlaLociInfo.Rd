\name{hlaLociInfo}
\alias{hlaLociInfo}
\title{
    HLA/KIR Locus Information
}
\description{
    To get the starting and ending positions in basepair of HLA/KIR loci.
}
\usage{
hlaLociInfo(assembly=c("auto", "auto-silent", "hg18", "hg19", "hg38",
    "unknown"))
}
\arguments{
    \item{assembly}{the human genome reference: "hg18", "hg19" (default),
        "hg38"; "auto" refers to "hg19"; "auto-silent" refers to "hg19" without
        any warning}
}
\value{
    Return a data frame include the genomic locations.
}

\references{
    NCBI Resources: \url{https://www.ncbi.nlm.nih.gov/gene},
    HLA Nomenclature: \url{http://hla.alleles.org/genes/index.html}
}
\author{Xiuwen Zheng}

\examples{
hlaLociInfo()
}

\keyword{HLA}
\keyword{genetics}
