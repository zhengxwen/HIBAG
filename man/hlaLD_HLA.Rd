\name{hlaLD_HLA}
\alias{hlaLD_HLA}
\title{
	Composite Linkage disequilibrium
}
\description{
	To calculate composite linkage disequilibrium between HLA locus and a SNP marker
}
\usage{
hlaLD_HLA(HLA, snp)
}
\arguments{
	\item{HLA}{an object of \code{\link{hlaAlleleClass}}}
	\item{snp}{a vector of SNP genotypes, 0 -- BB, 1 -- AB, 2 -- AA, NA -- missing}
}
\value{
	The value of linkage disequilibrium.
}
\references{
	Weir BS, Cockerham CC: Complete characterization of disequilibrium at two loci;
		in Feldman MW (ed): Mathematical Evolutionary Theory.
		Princeton, NJ: Princeton University Press, 1989.

	Zaykin, D. V., Pudovkin, A., and Weir, B. S. (2008).
		Correlation-based inference for linkage disequilibrium with multiple alleles.
		Genetics 180, 533-545.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaGenoLD}}
}

\examples{
# load HLA types and SNP genotypes
data(hla.type.table, package="HIBAG")
data(hla.snp.geno, package="HIBAG")

hla.id <- "A"
hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
	hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)

hlaLD_HLA(hla, hla.snp.geno$genotype[5, ])
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
