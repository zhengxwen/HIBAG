\name{hlaGenoLD}
\alias{hlaGenoLD}
\title{
	Composite Linkage Disequilibrium
}
\description{
	To calculate composite linkage disequilibrium (r2) between HLA locus
and SNP markers.
}
\usage{
hlaGenoLD(hla, geno)
}
\arguments{
	\item{hla}{an object of \code{\link{hlaAlleleClass}}}
	\item{geno}{an object of \code{\link{hlaSNPGenoClass}}, or a vector or
		matrix for SNP data}
}
\value{
	Return a vector of linkage disequilibrium (r2) for each SNP marker.
}
\references{
	Weir BS, Cockerham CC:
	Complete characterization of disequilibrium at two loci;
	in Feldman MW (ed): Mathematical Evolutionary Theory.
	Princeton, NJ: Princeton University Press, 1989.

	Zaykin, D. V., Pudovkin, A., and Weir, B. S. (2008).
	Correlation-based inference for linkage disequilibrium with multiple alleles.
	Genetics 180, 533-545.
}
\author{Xiuwen Zheng}

\examples{
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

# plot linkage disequilibrium
ymax <- 0.16
plot(NaN, NaN, xlab="SNP Position (in KB)",
	ylab="Composite Linkage Disequilibrium (r2)",
	xlim=range(HapMap_CEU_Geno$snp.position)/1000, ylim=c(0, ymax),
	main="Major Histocompatibility Complex")

hla.list <- c("A", "C", "DQA1")
col.list <- 1:3

# for-loop
for (i in 1:3)
{
	hla.id <- hla.list[i]

	# make a "hlaAlleleClass" object
	hla <- hlaAllele(HLA_Type_Table$sample.id,
		H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
		H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
		locus=hla.id, assembly="hg19")

	# linkage disequilibrium between HLA locus and SNP markers
	ld <- hlaGenoLD(hla, HapMap_CEU_Geno)

	# draw
	points(HapMap_CEU_Geno$snp.position/1000, ld, pch="*", col=i)
	x <- (hla$pos.start/1000 + hla$pos.end/1000)/2
	abline(v=x, col=col.list[i], lty=3, lwd=2.5)
	points(x, ymax, pch=25, col=7, bg=col.list[i], cex=1.5)
}
legend("topleft", col=col.list, pt.bg=col.list, text.col=col.list, pch=25,
	legend=paste("HLA -", hla.list))
}

\keyword{SNP}
\keyword{genetics}
