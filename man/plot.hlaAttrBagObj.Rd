\name{plot.hlaAttrBagObj}
\alias{plot.hlaAttrBagObj}
\title{
	Plot a "hlaAttrBagObj" object
}
\description{
	To show a scatterplot of the numbers of individual classifiers and SNP positions.
}
\usage{
\method{plot}{hlaAttrBagObj}(x, xlab=NULL, ylab=NULL, locus.color="red",
	locus.lty=2, locus.cex=1.25, ...)
}
\arguments{
	\item{x}{an object of \code{\link{hlaAttrBagObj}}}
	\item{xlab}{the label of X-axis}
	\item{ylab}{the label of Y-axis}
	\item{locus.color}{the color of text and line for HLA locus}
	\item{locus.lty}{the type of line for HLA locus}
	\item{locus.cex}{the font size of HLA locus}
	\item{...}{further arguments passed to or from other methods}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{plot.hlaAttrBagClass}},
	\code{\link{print.hlaAttrBagClass}}, \code{\link{print.hlaAttrBagObj}},
	\code{\link{summary.hlaAttrBagClass}}, \code{\link{summary.hlaAttrBagObj}}
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
