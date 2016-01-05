\name{hlaCheckAllele}
\alias{hlaCheckAllele}
\title{
    Check SNP alleles
}
\description{
    Check SNP reference and non-reference alleles.
}
\usage{
hlaCheckAllele(allele1, allele2)
}
\arguments{
    \item{allele1}{two alleles for the first individual, like
        \code{c("A/G", "C/G")}}
    \item{allele2}{two alleles for the second individual, like
        \code{c("A/G", "C/G")}}
}
\value{
    Return a logical vector, where \code{TRUE} indicates the alleles are
matching at that locus.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaCheckSNPs}}
}

\examples{
hlaCheckAllele(c("A/G", "T/G", "0/A"), c("G/A", "C/A", "G/0"))
}

\keyword{SNP}
\keyword{genetics}
