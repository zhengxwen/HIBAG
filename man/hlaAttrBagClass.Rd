\name{hlaAttrBagClass}
\alias{hlaAttrBagClass}
\title{The class of HIBAG model}
\description{
	The class of a HIBAG model, and its instance is returned from
\code{\link{hlaAttrBagging}}.
}
\value{
	Return a list of:
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
	\item{model}{internal use}
	\item{appendix}{an optional list:
		\code{platform} -- supported platform(s);
		\code{information} -- other information, like training sets, authors;
		\code{warning} -- any warning message}
}

\references{
	Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS;
	HIBAG -- HLA Genotype Imputation with Attribute Bagging;
	(Abstract 294, Platform/Oral Talk);
	Present at the 62nd Annual Meeting of the American Society of Human Genetics,
	November 9, 2012 in San Francisco, California.

	Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS;
	HIBAG -- HLA Genotype Imputation with Attribute Bagging.
	To appear in the Pharmacogenomics Journal.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{hlaAttrBagging}}, \code{\link{hlaParallelAttrBagging}},
	\code{\link{hlaAttrBagObj}}
}

\keyword{HLA}
\keyword{genetics}
