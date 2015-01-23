\name{hlaModelFromObj}
\alias{hlaModelFromObj}
\alias{hlaModelToObj}
\title{
    Conversion between the in-memory model and the object that can be
saved in a file
}
\description{
    Build a model \code{\link{hlaAttrBagClass}} from an object of
\code{\link{hlaAttrBagObj}} which is stored in an R object file, or convert
\code{\link{hlaAttrBagClass}} to \code{\link{hlaAttrBagObj}}.
}
\usage{
hlaModelFromObj(obj)
hlaModelToObj(model)
}
\arguments{
    \item{obj}{an object of \code{\link{hlaAttrBagObj}}}
    \item{model}{an object of \code{\link{hlaAttrBagClass}}}
}
\value{
    \code{hlaModelFromObj} returns \code{hlaAttrBagClass}, and
    \code{hlaModelToObj} returns \code{hlaAttrBagObj}.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaModelToObj}}, \code{\link{hlaAttrBagging}}
}

\examples{
# make a "hlaAlleleClass" object
hla.id <- "DQB1"
hla <- hlaAllele(HLA_Type_Table$sample.id,
    H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
    H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
    locus=hla.id, assembly="hg19")

# training genotypes
region <- 100   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
    hla.id, region*1000, assembly="hg19")
train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
    snp.sel = match(snpid, HapMap_CEU_Geno$snp.id))

# train a HIBAG model
set.seed(1000)
model <- hlaAttrBagging(hla, train.geno, nclassifier=2)
print(model)

mobj <- hlaModelToObj(model)

is(model)
is(mobj)


# close the HIBAG model explicitly
hlaClose(model)
}

\keyword{HLA}
\keyword{genetics}
