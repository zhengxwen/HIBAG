\name{hlaAlleleClass}
\alias{hlaAlleleClass}
\title{Class of HLA Type}
\description{
	The definition of a class for HLA types, returned from
\code{\link{hlaAllele}}.
}
\value{
	There are following components:
	\item{locus}{HLA locus}
	\item{pos.start}{the starting position in basepair}
	\item{pos.end}{the end position in basepair}
	\item{value}{a data frame}
	\item{assembly}{the human genome reference, such like "hg19"}
	\item{postprob}{a matrix of all posterior probabilities}

	The component \code{value} includes:
	\item{sample.id}{sample ID}
	\item{allele1}{HLA allele}
	\item{allele2}{HLA allele}
	\item{prob}{the posterior probability}
}

\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAllele}}
}

\keyword{HLA}
\keyword{genetics}
