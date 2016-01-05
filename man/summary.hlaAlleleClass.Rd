\name{summary.hlaAlleleClass}
\alias{summary.hlaAlleleClass}
\alias{summary.hlaAASeqClass}
\title{
    Summarize a ``hlaAlleleClass'' or ``hlaAASeqClass'' object
}
\description{
    Show the information of a \code{\link{hlaAlleleClass}} or
\code{\link{hlaAASeqClass}} object.
}
\usage{
\method{summary}{hlaAlleleClass}(object, verbose=TRUE, ...)
\method{summary}{hlaAASeqClass}(object, poly.only=TRUE, head=0L,
    verbose=TRUE, ...)
}
\arguments{
    \item{object}{an object of \code{\link{hlaAlleleClass}} or
        \code{\link{hlaAASeqClass}}}
    \item{poly.only}{if TRUE, only show the amino acid positions with
        polymorphism; otherwise, show all sequences}
    \item{head}{show the first \code{head} rows of cross tabulation, or
        \code{0L} for all rows}
    \item{verbose}{if TRUE, show information}
    \item{...}{further arguments passed to or from other methods}
}
\value{
    Return a \code{data.frame} of count and frequency for each HLA allele, if
\code{object} is \code{hlaAlleleClass}; a matrix of cross tabulation of amino
acids at each position, if \code{object} is \code{hlaAASeqClass}.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaAllele}}, \code{\link{hlaConvSequence}}
}

\keyword{HLA}
\keyword{genetics}
