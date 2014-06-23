\name{summary.hlaAttrBagObj}
\alias{summary.hlaAttrBagObj}
\title{
	Summarize a ``hlaAttrBagObj'' object
}
\description{
	Show the information of a \code{\link{hlaAttrBagObj}} object
}
\usage{
\method{summary}{hlaAttrBagObj}(object, show=TRUE, ...)
}
\arguments{
	\item{object}{an object of \code{\link{hlaAttrBagObj}}}
	\item{show}{if TRUE, show information}
	\item{...}{further arguments passed to or from other methods}
}
\value{
	Return a \code{list}:
	\item{num.classifier}{the total number of classifiers}
	\item{num.snp}{the total number of SNPs}
	\item{snp.id}{SNP IDs}
	\item{snp.position}{SNP position in basepair}
	\item{snp.hist}{the number of classifier for each SNP, and it could be
		used for SNP importance}
	\item{info}{a \code{data.frame} for the average number of SNPs
		(\code{num.snp}), haplotypes (\code{num.haplo}), out-of-bag accuracies
		(\code{accuracy}) among all classifiers:
		mean, standard deviation, min, max}
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{plot.hlaAttrBagObj}}, \code{\link{print.hlaAttrBagObj}}
}

\keyword{HLA}
\keyword{genetics}
