\name{hlaMakeSNPGeno}
\alias{hlaMakeSNPGeno}
\title{
	Make a SNP genotype object
}
\description{
	To create a \code{\link{hlaSNPGenoClass}} object (SNP genotype object).
}
\usage{
hlaMakeSNPGeno(genotype, sample.id, snp.id, snp.position, A.allele, B.allele)
}
\arguments{
	\item{genotype}{a genotype matrix, ``# of SNPs'' - by - ``# of individuals''}
	\item{sample.id}{a vector of sample IDs}
	\item{snp.id}{a vector of SNP IDs}
	\item{snp.position}{a vector of SNP positions}
	\item{A.allele}{a vector of A alleles in the SNP list}
	\item{B.allele}{a vector of B alleles in the SNP list}
}
\details{
	\code{genotype} is a numeric matrix, with an entry value 0 standing for BB (ZERO A allele),
1 for AB (ONE A allele), 2 for AA (TWO A alleles) and others for missing values (missing genotypes
are usually set to be NA).
}
\value{
	Return a \code{\link{hlaSNPGenoClass}} object, and it is a list:
	\item{genotype}{a genotype matrix, ``# of SNPs'' - by - ``# of individuals''}
	\item{sample.id}{a vector of sample IDs}
	\item{snp.id}{a vector of SNP IDs}
	\item{snp.position}{a vector of SNP positions in basepair}
	\item{snp.allele}{a vector of characters with the format of ``A allele/B allele''}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{hlaMakeSNPHaplo}}, \code{\link{hlaGenoSubset}}, \code{\link{hlaHaploSubset}},
\code{\link{hlaGenoCombine}}
}

\examples{
# load SNP genotypes
data(hla.snp.geno, package="HIBAG")
summary(hla.snp.geno)

allele <- strsplit(hla.snp.geno$snp.allele, "/")
A.allele <- sapply(allele, function(x) { x[1] })
B.allele <- sapply(allele, function(x) { x[2] })

geno <- hlaMakeSNPGeno(hla.snp.geno$genotype, hla.snp.geno$sample.id, hla.snp.geno$snp.id,
	hla.snp.geno$snp.position, A.allele, B.allele)
summary(geno)
}

\keyword{HLA}
\keyword{GWAS}
\keyword{genetics}
