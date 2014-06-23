\name{hlaModelToObj}
\alias{hlaModelToObj}
\title{
	Convert ``hlaAttrBagClass'' to ``hlaAttrBagObj''
}
\description{
	Convert \code{\link{hlaAttrBagClass}} to \code{\link{hlaAttrBagObj}}, and
\code{\link{hlaAttrBagObj}} can be written to the specified file with an
external representation of R objects.
}
\usage{
hlaModelToObj(model)
}
\arguments{
	\item{model}{an object of \code{\link{hlaAttrBagClass}}}
}
\value{
	Return an object of \code{\link{hlaAttrBagObj}}.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaModelFromObj}}, \code{\link{hlaAttrBagging}}
}

\keyword{HLA}
\keyword{genetics}
