\name{hlaAssocTest}
\alias{hlaAssocTest}
\title{
    Statistical Association Tests
}
\description{
    To perform statistical association tests via Pearson's Chi-squared test,
Fisher's exact test, logistic regression, random effect models.
}
\usage{
hlaAssocTest(hla, formula, data,
    model=c("dominant", "additive", "recessive", "genotype"),
    model.fit=c("glm"), showOR=FALSE, ...)
}
\arguments{
    \item{hla}{an object of \code{\link{hlaAlleleClass}}}
    \item{formula}{an object of class \code{"formula"} (or one that can be
        coerced to that class): a symbolic description of the model to be
        fitted, e.g., \code{~group}}
    \item{data}{an optional data frame, list or environment containing the
        variables in the model. If not found in \code{data}, the variables
        are taken from \code{environment(formula)}}
    \item{model}{dominant, additive, recessive or genotype models:
        \code{"dominant"} is default}
    \item{model.fit}{"glm" -- generalized linear regression}
    \item{showOR}{show odd ratio (OR) instead of log OR if \code{TRUE}}
    \item{...}{optional arguments to \code{\link{glm}} or \code{nlme} call}
}
\details{
\tabular{ll}{
    \bold{model} \tab \bold{description (given a specific HLA allele h)} \cr
    dominant  \tab [-/-] vs. [-/h,h/h] \cr
    additive  \tab the allele dosage, 0: -/-, 1: -/h, 2: h/h \cr
    recessive \tab [-/-,-/h] vs. [h/h] \cr
    genotype  \tab [-/-], [-/h], [h/h] \cr
}
}
\value{
    Return a \code{data.frame} with
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
    \item{chisq.st}{the value the chi-squared test statistic}
    \item{chisq.p}{the p-value for the Chi-squared test}
    \item{fisher.p}{the p-value for the Fisher's exact test}
  
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{summary.hlaAttrBagClass}}
}

\examples{
# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
    H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
    H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
    locus=hla.id, assembly="hg19")

set.seed(1000)
n <- nrow(hla$value)
dat <- data.frame(case = as.factor(c(rep(0, n/2), rep(1, n/2))),
    y = rnorm(n), pc1 = rnorm(n))

hlaAssocTest(hla, case ~ 1, data=dat)
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
