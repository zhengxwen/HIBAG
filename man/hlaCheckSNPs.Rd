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
hlaCheckSNPs(model, object, match.pos=TRUE, verbose=TRUE)
}
\arguments{
	\item{model}{an object of \code{\link{hlaAttrBagClass}}, or an object of
		\code{\link{hlaAttrBagObj}}}
	\item{object}{a genotype object of \code{\link{hlaSNPGenoClass}}, or a character
		vector like c("rs2523442", "rs9257863", ...)}
	\item{match.pos}{if \code{TRUE}, it uses both SNP ID and position to match SNPs;
		otherwise, it uses SNP ID only; calling \code{hlaSNPID(, with.pos=match.pos)}}
	\item{verbose}{if TRUE, show information}
}
\value{
	return a \code{data.frame} for individual classifiers:
	\item{NumOfValidSNP}{the number of non-missing SNPs in an individual classifier}
	\item{NumOfSNP}{the number of SNP predictors in an individual classifier}
	\item{fraction}{NumOfValidSNP / NumOfSNP}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{predict.hlaAttrBagClass}}
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
