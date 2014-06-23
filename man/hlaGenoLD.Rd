\name{hlaGenoLD}
\alias{hlaGenoLD}
\title{
	Composite Linkage disequilibrium
}
\description{
	To calculate composite linkage disequilibrium between HLA locus and SNP markers
}
\usage{
hlaGenoLD(hla, geno)
}
\arguments{
	\item{hla}{an object of \code{\link{hlaAlleleClass}}}
	\item{geno}{an object of \code{\link{hlaSNPGenoClass}}}
}
\value{
	Return a vector of linkage disequilibrium for each SNP marker.
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
	\code{\link{hlaLD_HLA}}
}

\examples{
# load HLA types and SNP genotypes
data(hla.type.table, package="HIBAG")
data(hla.snp.geno, package="HIBAG")

# plot linkage disequilibrium
ymax <- 0.16
plot(NaN, NaN, xlab="SNP Position (in KB)", ylab="Linkage Disequilibrium",
	xlim=range(hla.snp.geno$snp.position)/1000, ylim=c(0, ymax),
	main="Major Histocompatibility Complex")
# for-loop
for (i in 1:3)
{
	hla.id <- c("A", "C", "DQA1")[i]
	# make a "hlaAlleleClass" object
	hla <- hlaAllele(hla.type.table$sample.id, hla.type.table[, paste(hla.id, ".1", sep="")],
		hla.type.table[, paste(hla.id, ".2", sep="")], locus=hla.id)
	# linkage disequilibrium between HLA locus and SNP markers
	ld <- hlaGenoLD(hla, hla.snp.geno)
	# draw
	points(hla.snp.geno$snp.position/1000, ld, pch="*", col=i)
	x <- (hla$pos.start/1000 + hla$pos.end/1000)/2
	abline(v=x, col=i, lty=3)
	text(x, ymax-(i-1)*0.005, hla.id, adj=1, col=i)
}
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
