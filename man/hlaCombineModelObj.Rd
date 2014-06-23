\name{hlaCombineModelObj}
\alias{hlaCombineModelObj}
\title{
	Combine two objects of "hlaAttrBagObj"
}
\description{
	Merge two objects of "hlaAttrBagObj" together, which is useful for building
an ensemble model in parallel.
}
\usage{
hlaCombineModelObj(obj1, obj2)
}
\arguments{
	\item{obj1}{an object of \code{\link{hlaAttrBagObj}}}
	\item{obj2}{an object of \code{\link{hlaAttrBagObj}}}
}
\value{
	Return an object of \code{\link{hlaAttrBagObj}}.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{hlaClusterAttrBagging}}
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
