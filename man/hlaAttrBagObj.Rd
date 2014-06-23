\name{hlaAttrBagObj}
\alias{hlaAttrBagObj}
\title{The class of HIBAG object}
\description{
	The class of a HIBAG model.
}
\value{
	A list of:
	\item{n.samp}{the total number of training samples}
	\item{n.snp}{the total number of candidate SNP predictors}
	\item{sample.id}{the sample IDs}
	\item{snp.id}{the SNP IDs}
	\item{snp.position}{SNP position in basepair}
	\item{snp.allele}{a vector of characters with the format of ``A allele/B allele''}
	\item{snp.allele.freq}{the allele frequencies}
	\item{hla.locus}{the name of HLA locus}
	\item{hla.allele}{the HLA alleles used in the model}
	\item{hla.freq}{the HLA allele frequencies}
	\item{assembly}{the human genome reference, such like "hg19"}
	\item{classifiers}{a list of all classifiers (described as follows)}
	\item{appendix}{
		\code{platform} -- supported platform(s);
		\code{information} -- other information, like training sets, authors;
		\code{warning} -- any warning message}

	\code{classifiers} has the following components:
	\item{samp.num}{the number of copies of samples in a bootstrap sample}
	\item{haplos}{a data.frame of haplotype frequencies}
	\item{ }{\code{freq}  -- haplotype frequency}
	\item{ }{\code{hla}   -- a HLA allele}
	\item{ }{\code{haplo} -- a SNP haplotype, with an entry value 0 standing
		for B (ZERO A allele), 1 for A (ONE A allele)}
	\item{snpidx}{the SNP indices used in this classifier}
	\item{outofbag.acc}{the out-of-bag accuracy of this classifier}
}

\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{hlaParallelAttrBagging}},
	\code{\link{hlaModelToObj}}, \code{\link{hlaModelFiles}},
	\code{\link{hlaAttrBagClass}}
}

\keyword{HLA}
\keyword{genetics}
