\name{hlaDistance}
\alias{hlaDistance}
\title{
    Distance matrix of HLA alleles
}
\description{
    To calculate the distance matrix of HLA alleles from a HIBAG model.
}
\usage{
hlaDistance(model)
}
\arguments{
    \item{model}{a model of \code{\link{hlaAttrBagClass}} or
        \code{\link{hlaAttrBagObj}}}
}
\value{
    Return a distance matrix with row and column names for HLA alleles.
}
\author{Xiuwen Zheng}

\examples{
# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
    H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
    H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
    locus=hla.id, assembly="hg19")

# flanking genotypes
train.geno <- hlaGenoSubsetFlank(HapMap_CEU_Geno, hla.id, 500000)
summary(train.geno)

# train a HIBAG model
set.seed(100)
model <- hlaAttrBagging(hla, train.geno, nclassifier=10)
summary(model)

# distance matrix
d <- hlaDistance(model)

# draw
p <- hclust(as.dist(d))
plot(p, xlab="HLA alleles")
}

\keyword{HLA}
\keyword{genetics}
