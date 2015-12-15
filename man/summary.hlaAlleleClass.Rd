\name{summary.hlaAlleleClass}
\alias{summary.hlaAlleleClass}
\alias{summary.hlaSeqClass}
\title{
    Summarize a ``hlaAlleleClass'' or ``hlaSeqClass'' object
}
\description{
    Show the information of a \code{\link{hlaAlleleClass}} or
\code{\link{hlaSeqClass}} object.
}
\usage{
\method{summary}{hlaAlleleClass}(object, verbose=TRUE, ...)
\method{summary}{hlaSeqClass}(object, head=100L, verbose=TRUE, ...)
}
\arguments{
    \item{object}{an object of \code{\link{hlaAlleleClass}} or
        \code{\link{hlaSeqClass}}}
    \item{head}{show the first \code{head} rows of cross tabulation}
    \item{verbose}{if TRUE, show information}
    \item{...}{further arguments passed to or from other methods}
}
\value{
    Return a \code{data.frame} of count and frequency for each HLA allele, if
\code{object} is \code{hlaAlleleClass}; a matrix of cross tabulation of amino
acids at each position, if \code{object} is \code{hlaSeqClass}.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaAllele}}, \code{\link{hlaConvSequence}}
}

\keyword{HLA}
\keyword{genetics}
