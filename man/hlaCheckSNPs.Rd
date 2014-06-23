\name{hlaCheckSNPs}
\alias{hlaCheckSNPs}
\title{
	Check the SNP predictors in a HIBAG model
}
\description{
	Check the SNP predictors in a HIBAG model, by calculating the overlapping between
the model and SNP genotypes.
}
\usage{
hlaCheckSNPs(model, object,
	match.type=c("RefSNP+Position", "RefSNP", "Position"), verbose=TRUE)
}
\arguments{
	\item{model}{an object of \code{\link{hlaAttrBagClass}}, or an object of
		\code{\link{hlaAttrBagObj}}}
	\item{object}{a genotype object of \code{\link{hlaSNPGenoClass}}, or a character
		vector like c("rs2523442", "rs9257863", ...)}
	\item{match.type}{\code{"RefSNP+Position"} (by default) -- using both of
		RefSNP IDs and positions; \code{"RefSNP"} -- using RefSNP IDs only;
		\code{"Position"} -- using positions only}
	\item{verbose}{if TRUE, show information}
}
\value{
	return a \code{data.frame} for individual classifiers:
	\item{NumOfValidSNP}{the number of non-missing SNPs in an individual classifier}
	\item{NumOfSNP}{the number of SNP predictors in an individual classifier}
	\item{fraction}{NumOfValidSNP / NumOfSNP}
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{predict.hlaAttrBagClass}}
}

\keyword{SNP}
\keyword{genetics}
