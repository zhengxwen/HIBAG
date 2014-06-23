\name{predict.hlaAttrBagClass}
\alias{predict.hlaAttrBagClass}
\title{
	HIBAG model prediction
}
\description{
	To predict HLA type based on an attribute bagging model.
}
\usage{
\method{predict}{hlaAttrBagClass}(object, genotypes, type=c("response", "prob"),
	allele.check=TRUE, match.pos=TRUE, verbose=TRUE, ...)
}
\arguments{
	\item{object}{a model of \code{\link{hlaAttrBagClass}}}
	\item{genotypes}{a genotype object of \code{\link{hlaSNPGenoClass}}}
	\item{type}{"response": return the best-guess type plus posterior probability;
		"prob": return all posterior probabilites}
	\item{allele.check}{if \code{TRUE}, check and then switch allele pairs if needed}
	\item{match.pos}{if \code{TRUE}, it uses both SNP ID and position to match SNPs;
		otherwise, it uses SNP ID only; calling \code{hlaSNPID(, with.pos=match.pos)}}
	\item{verbose}{if TRUE, show information}
	\item{...}{further arguments passed to or from other methods}
}
\value{
	Return a \code{\link{hlaAlleleClass}} object with posterior probabilities of
predicted HLA types, or a matrix with all posterior probabilities.
	If a probability matrix is returned, \code{colnames} is \code{sample.id} and
\code{rownames} is an unordered pair of HLA alleles.
}
\details{
	If more than 50\% of SNP predictors are missing, a warning will be given.

	When \code{match.pos=TRUE}, the matching of SNPs requires both SNP ID and position,
and a lower missing fraction maybe gained by matching reference SNP ID only. Call
\code{predict(, match.pos=FALSE)} for this purpose. However, the SNPs with identical
reference SNP ID but different position may raise a problem of SNP mismatching.
Any concern about SNP mismatching should be emailed to the genotyping platform provider.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{hlaAllele}}
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
