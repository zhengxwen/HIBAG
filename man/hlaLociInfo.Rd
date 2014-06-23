\name{hlaLociInfo}
\alias{hlaLociInfo}
\title{
	HLA Locus Information
}
\description{
	To get the starting and ending positions in basepair of HLA loci.
}
\usage{
hlaLociInfo(assembly=c("auto", "hg18", "hg19", "unknown"))
}
\arguments{
	\item{assembly}{the human genome reference: "hg19" (default), "hg18",
		"auto" refers to "hg19"}
}
\value{
	Return a list:
	\item{loci}{the names of HLA genes: A, B, C, DRB1, DRB5, DQA1, DQB1
		and DPB1}
	\item{pos.HLA.start}{the starting position in basepair}
	\item{pos.HLA.end}{the ending position in basepair}
	\item{length.HLA}{the length of HLA genes in basepair}
}
\author{Xiuwen Zheng}

\examples{
hlaLociInfo()
}

\keyword{HLA}
\keyword{genetics}
