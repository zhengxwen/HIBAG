\name{hlaAlleleClass}
\alias{hlaAlleleClass}
\title{Class of HLA Type}
\description{
	The definition of a class for HLA types.
}
\value{
	There are following components:
	\item{locus}{HLA locus}
	\item{pos.start}{the starting position in basepair}
	\item{pos.end}{the end position in basepair}
	\item{value}{a data frame}
	\item{assembly}{genome assembly, such like "hg19"}

	The component \code{value} includes:
	\item{sample.id}{sample ID}
	\item{allele1}{HLA allele}
	\item{allele2}{HLA allele}
	\item{prob}{the posterior probability}
}

\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaMakeSNPGeno}}
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
