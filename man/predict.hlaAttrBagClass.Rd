\name{predict.hlaAttrBagClass}
\alias{predict.hlaAttrBagClass}
\title{
	HIBAG model predictions
}
\description{
	To predict HLA type based on an attribute bagging model.
}
\usage{
\method{predict}{hlaAttrBagClass}(object, genotypes, type=c("response", "prob"),
	allele.check=TRUE, verbose=TRUE, ...)
}
\arguments{
	\item{object}{a model of \code{\link{hlaAttrBagClass}}}
	\item{genotypes}{a genotype object of \code{\link{hlaSNPGenoClass}}}
	\item{type}{"response": return the best-guess type plus posterior probability;
		"prob": return all posterior probabilites}
	\item{allele.check}{if TRUE, check and switch allele pairs}
	\item{verbose}{if TRUE, show information}
	\item{...}{further arguments passed to or from other methods}
}
\value{
	Return a \code{\link{hlaAlleleClass}} object with posterior probabilities of predicted
HLA types, or a matrix with all posterior probabilities. If a probability matrix is returned,
\code{colnames} is \code{sample.id} and \code{rownames} is an unorderd pair of HLA alleles.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{hlaAllele}}
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
