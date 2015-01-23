\name{hlaLociInfo}
\alias{hlaLociInfo}
\title{
    HLA Locus Information
}
\description{
    To get the starting and ending positions in basepair of HLA loci.
}
\usage{
hlaLociInfo(assembly=c("auto", "auto-silent", "hg18", "hg19", "hg20",
    "unknown"))
}
\arguments{
    \item{assembly}{the human genome reference: "hg18", "hg19" (default),
        "hg20"; "auto" refers to "hg19"; "auto-silent" refers to "hg19" without
        any warning}
}
\value{
    Return a \code{\link[GenomicRanges]{GRanges}} object (a container for
the genomic locations and their associated annotations).
}

\references{
    NCBI Resources: \url{http://www.ncbi.nlm.nih.gov/gene}
}
\author{Xiuwen Zheng}

\examples{
hlaLociInfo()
}

\keyword{HLA}
\keyword{genetics}
