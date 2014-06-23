\name{hlaLociInfo}
\alias{hlaLociInfo}
\title{
	HLA Locus Information
}
\description{
	To get the starting and ending positions in basepair of HLA loci.
}
\usage{
hlaLociInfo(assembly=c("auto", "hg19", "hg18", "NCBI37", "NCBI36"))
}
\arguments{
	\item{assembly}{genome assembly: "hg19" (default), "hg18", "NCBI37" or "NCBI36";
		"auto" refers to "hg19"}
}
\details{
	HLA loci are "A", "B", "C", "DRB1", "DRB5", "DQA1", "DQB1" and "DPB1".
}
\value{
	Return a list:
	\item{loci}{the names of HLA genes}
	\item{pos.HLA.start}{the starting position in basepair}
	\item{pos.HLA.end}{the ending position in basepair}
	\item{length.HLA}{the length of HLA genes in basepair}
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAllele}}, \code{\link{hlaAlleleSubset}}, \code{\link{hlaAlleleClass}}
}

\examples{
hlaLociInfo()
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
