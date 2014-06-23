\name{hlaClose}
\alias{hlaClose}
\title{
	Dispose a model object
}
\description{
	Release all resources stored in the \code{\link{hlaAttrBagClass}} object.
The HIBAG package allows up to 256 \code{\link{hlaAttrBagClass}} objects
stored in memory.
}
\usage{
hlaClose(model)
}
\arguments{
	\item{model}{an object of \code{\link{hlaAttrBagClass}}}
}
\value{
	None.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{summary.hlaAttrBagClass}}
}

\keyword{HLA}
\keyword{genetics}
