\name{hlaLociInfo}
\alias{hlaLociInfo}
\title{
	HLA Locus Information
}
\description{
	To get the starting and ending positions in basepair of HLA loci.
}
\usage{
hlaLociInfo()
}
\details{
	HLA loci are "A", "B", "C", "DRB1", "DRB5", "DQA1", "DQB1" and "DPB1".
}
\value{
	Return a list:
	\item{loci}{the names of HLA genes}
	\item{pos.HLA.start}{the starting positions in basepair}
	\item{pos.HLA.end}{the ending positions in basepair}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaAllele}}, \code{\link{hlaAlleleSubset}}, \code{\link{hlaAlleleClass}}
}

\examples{
hlaLociInfo()
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
