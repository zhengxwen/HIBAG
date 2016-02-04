\name{hlaAssocTest}
\alias{hlaAssocTest}
\alias{hlaAssocTest.hlaAlleleClass}
\alias{hlaAssocTest.hlaAASeqClass}
\title{
    Statistical Association Tests
}
\description{
    Perform statistical association tests via Pearson's Chi-squared test,
Fisher's exact test and logistic regressions.
}
\usage{
\method{hlaAssocTest}{hlaAlleleClass}(hla, formula, data,
    model=c("dominant", "additive", "recessive", "genotype"),
    model.fit=c("glm"), prob.threshold=NaN, use.prob=FALSE, showOR=FALSE,
    verbose=TRUE, ...)
\method{hlaAssocTest}{hlaAASeqClass}(hla, formula, data,
    model=c("dominant", "additive", "recessive", "genotype"),
    model.fit=c("glm"), prob.threshold=NaN, use.prob=FALSE, showOR=FALSE,
    show.all=FALSE, verbose=TRUE, ...)
}
\arguments{
    \item{hla}{an object of \code{\link{hlaAlleleClass}}}
    \item{formula}{an object of class \code{"formula"} (or one that can be
        coerced to that class): a symbolic description of the model to be
        fitted, e.g., \code{y ~ 1}, \code{y ~ h + a}}
    \item{data}{an optional data frame, list or environment containing the
        variables in the model. If not found in \code{data}, the variables
        are taken from \code{environment(formula)}}
    \item{model}{dominant, additive, recessive or genotype models:
        \code{"dominant"} is default}
    \item{model.fit}{"glm" -- generalized linear regression}
    \item{prob.threshold}{the probability threshold to exclude individuals
        with low confidence scores}
    \item{use.prob}{if \code{TRUE}, use the posterior probabilities as weights
        in \code{glm} models}
    \item{showOR}{show odd ratio (OR) instead of log OR if \code{TRUE}}
    \item{show.all}{if \code{TRUE}, show both significant and non-significant
        results; if \code{FALSE}, only show significant results}
    \item{verbose}{if TRUE, show information}
    \item{...}{optional arguments to \code{\link{glm}} or \code{nlme} call}
}
\details{
\tabular{ll}{
    \bold{model} \tab \bold{description (given a specific HLA allele h)} \cr
    dominant  \tab [-/-] vs. [-/h,h/h] (0 vs. 1 in design matrix) \cr
    additive  \tab [-] vs. [h] in Chi-squared and Fisher's exact test,
        the allele dosage in regressions (0: -/-, 1: -/h, 2: h/h) \cr
    recessive \tab [-/-,-/h] vs. [h/h] (0 vs. 1 in design matrix) \cr
    genotype  \tab [-/-], [-/h], [h/h] (0 vs. 1 in design matrix) \cr
}

    In allelic associations, Chi-squared and Fisher exact tests are preformed
on the cross tabulation, which is constructed according to the specified
model (dominant, additive, recessive and gneotype).

    In amino acid associations, Fisher exact test is performed on a cross
tabulation with the numbers of each amino acid stratified by response variable
(e.g., disease status).

    In linear and logistic regressions, 95\% confidence intervals are
calculated based on asymptotic normality. The option \code{use.prob=TRUE} might
be useful in the sensitivity analysis.
}
\value{
    Return a \code{data.frame} with
    \item{[-]}{the number of haplotypes not carrying the specified HLA allele}
    \item{[h]}{the number of haplotype carrying the specified HLA allele}
    \item{\%.[-], ...}{case/disease proportion in the group [-], ...}
    \item{[-/-]}{the number of individuals or haplotypes not carrying the
        specified HLA allele}
    \item{[-/h]}{the number of individuals or haplotypes carrying one
        specified HLA allele}
    \item{[-/h]}{the number of individuals or haplotypes carrying two
        specified HLA alleles}
    \item{[-/h,h/h]}{the number of individuals or haplotypes carrying one or
        two specified HLA alleles}
    \item{[-/-,-/h]}{the number of individuals or haplotypes carrying at most
        one specified HLA allele}
    \item{\%.[-/-], ...}{case/disease proportion in the group [-/-], ...}
    \item{avg.[-/-], ...}{outcome average in the group [-/-], ...}
    \item{chisq.st}{the value the chi-squared test statistic}
    \item{chisq.p}{the p-value for the Chi-squared test}
    \item{fisher.p}{the p-value for the Fisher's exact test} 
    \item{h.est}{the coefficient estimate of HLA allele}
    \item{h.25\%, h.75\%}{the 95\% confidence interval for HLA allele}
    \item{h.pval}{p value for HLA allele}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaConvSequence}}, \code{\link{summary.hlaAASeqClass}}
}

\examples{
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
    H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
    H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
    locus=hla.id, assembly="hg19")

set.seed(1000)
n <- nrow(hla$value)
dat <- data.frame(case = c(rep(0, n/2), rep(1, n/2)), y = rnorm(n),
    pc1 = rnorm(n))

hlaAssocTest(hla, case ~ 1, data=dat)
hlaAssocTest(hla, case ~ 1, data=dat, model="additive")
hlaAssocTest(hla, case ~ 1, data=dat, model="recessive")
hlaAssocTest(hla, case ~ 1, data=dat, model="genotype")

hlaAssocTest(hla, y ~ 1, data=dat)
hlaAssocTest(hla, y ~ 1, data=dat, model="genotype")

hlaAssocTest(hla, case ~ h, data=dat)
hlaAssocTest(hla, case ~ h + pc1, data=dat)
hlaAssocTest(hla, case ~ h + pc1, data=dat, showOR=TRUE)

hlaAssocTest(hla, y ~ h, data=dat)
hlaAssocTest(hla, y ~ h + pc1, data=dat)
hlaAssocTest(hla, y ~ h + pc1, data=dat, showOR=TRUE)

hlaAssocTest(hla, case ~ h, data=dat, model="additive")
hlaAssocTest(hla, case ~ h, data=dat, model="recessive")
hlaAssocTest(hla, case ~ h, data=dat, model="genotype")
}

\keyword{HLA}
\keyword{SNP}
\keyword{genetics}
