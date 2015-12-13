\name{hlaSeqClass}
\alias{hlaSeqClass}
\title{Class of HLA Sequence Type}
\description{
    The definition of a class for HLA sequence genotype.
}
\value{
    There are following components:
    \item{locus}{HLA locus}
    \item{pos.start}{the starting position in basepair}
    \item{pos.end}{the end position in basepair}
    \item{value}{a data frame}
    \item{assembly}{the human genome reference, such like "hg19"}

    The component \code{value} includes:
    \item{sample.id}{sample ID}
    \item{allele1}{amino acid or nucleotide sequence}
    \item{allele2}{amino acid or nucleotide sequence}
    \item{P1, ..., Pn}{if applicable, a matrix of posterior probability,
        row -- sample, column -- position of amino acid}
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaConvSequence}}
}

\keyword{HLA}
\keyword{genetics}
