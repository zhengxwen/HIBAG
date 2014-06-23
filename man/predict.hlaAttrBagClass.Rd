\name{predict.hlaAttrBagClass}
\alias{predict.hlaAttrBagClass}
\title{
	HIBAG model prediction (in parallel)
}
\description{
	To predict HLA type based on a HIBAG model (in parallel).
}
\usage{
\method{predict}{hlaAttrBagClass}(object, snp, cl,
	type=c("response", "prob", "response+prob"), vote=c("prob", "majority"),
	allele.check=TRUE, match.type=c("RefSNP+Position", "RefSNP", "Position"),
	same.strand=FALSE, verbose=TRUE, ...)
}
\arguments{
	\item{object}{a model of \code{\link{hlaAttrBagClass}}}
	\item{snp}{a genotypic object of \code{\link{hlaSNPGenoClass}}}
	\item{cl}{a cluster object, created by the package \link[parallel]{parallel}
		or \href{http://CRAN.R-project.org/package=snow}{snow}; if \code{NULL}
		is given, a uniprocessor implementation will be performed}
	\item{type}{"response": return the best-guess type plus its posterior
		probability; "prob": return all posterior probabilities;
		"response+prob": return the best-guess and all posterior probabilities}
	\item{vote}{\code{"prob"} (default behavior) -- make a prediction based on
		the averaged posterior probabilities from all individual classifiers;
		\code{"majority"} -- majority voting from all individual
		classifiers, where each classifier votes for an HLA type}
	\item{allele.check}{if \code{TRUE}, check and then switch allele pairs
		if needed}
	\item{match.type}{\code{"RefSNP+Position"} (by default) -- using both of
		RefSNP IDs and positions; \code{"RefSNP"} -- using RefSNP IDs only;
		\code{"Position"} -- using positions only}
	\item{same.strand}{\code{TRUE} assuming alleles are on the same strand
		(e.g., forward strand); otherwise, \code{FALSE} not assuming whether
		on the same strand or not}
	\item{verbose}{if TRUE, show information}
	\item{...}{further arguments passed to or from other methods}
}
\value{
	Return a \code{\link{hlaAlleleClass}} object with posterior probabilities
of predicted HLA types, or a matrix of pairwise possible HLA types with all
posterior probabilities. If \code{type = "response+prob"}, return a
\code{\link{hlaAlleleClass}} object with a matrix of \code{postprob} for the
probabilities of all pairs of alleles.
	If a probability matrix is returned, \code{colnames} is \code{sample.id}
and \code{rownames} is an unordered pair of HLA alleles.
}
\details{
	If more than 50\% of SNP predictors are missing, a warning will be given.

	When \code{match.type="RefSNP+Position"}, the matching of SNPs requires
both RefSNP IDs and positions. A lower missing fraction maybe gained by
matching RefSNP IDs or positions only. Call
	\code{predict(..., match.type="RefSNP")} or
	\code{predict(..., match.type="Position")} for this purpose.
It might be safe to assume that the SNPs with the same positions on the same
genome reference (e.g., hg19) are the same variant albeit the different RefSNP
IDs. Any concern about SNP mismatching should be emailed to the genotyping
platform provider.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{hlaAllele}},
	\code{\link{hlaCompareAllele}}, \code{\link{hlaParallelAttrBagging}}
}

\keyword{HLA}
\keyword{SNP}
\keyword{genetics}
