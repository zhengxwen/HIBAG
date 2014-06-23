\name{hlaClose}
\alias{hlaClose}
\title{
	Dispose a model object
}
\description{
	Release all resources stored in the "hlaAttrBagClass" object. HIBAG allows up to 256
"hlaAttrBagClass" objects.
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
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{summary.hlaAttrBagClass}}
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
