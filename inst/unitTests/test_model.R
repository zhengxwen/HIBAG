# ===========================================================================
#
# runTests.R: unit testing
#

CreateModel <- function()
{
    # load the package HIBAG
    library(HIBAG)

    # load HLA types and SNP genotypes
    data(HLA_Type_Table, package="HIBAG")
    data(HapMap_CEU_Geno, package="HIBAG")

    # set seed
    RNGkind("L'Ecuyer-CMRG")
    set.seed(1000)

    ans <- list()

    # make HIBAG models
    for (hla.id in c("A", "B", "C", "DQA1", "DQB1", "DRB1"))
    {
        train.hla <- hlaAllele(HLA_Type_Table$sample.id,
            H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
            H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
            locus=hla.id, assembly="hg19")

        # SNP predictors within the flanking region on each side
        region <- 500   # kb
        snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id,
            HapMap_CEU_Geno$snp.position, hla.id, region*1000, assembly="hg19")

        # training genotypes
        train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
            snp.sel=match(snpid, HapMap_CEU_Geno$snp.id))

        # train a HIBAG model
        model <- hlaAttrBagging(train.hla, train.geno, nclassifier=10)
        summary(model)
        mobj <- hlaModelToObj(model)

        ans[[hla.id]] <- mobj
    }

    ans
}


test.HIBAG <- function()
{
    ##  unit testing

    library(HIBAG)

    SSE <- .Call("HIBAG_SSE_Flag", PACKAGE="HIBAG")
    if (SSE[1] > 0)
    {
        # scientific reproducibility with random number
        # when SSE2 is used since the model file was created with SSE2
        # more details: x87 might use 80-bit precision interally but
        #   SSE2 uses 64-bit precision only

        # the pre-fit models
        fn <- system.file("extdata", "ModelListTest.RData", package="HIBAG")
        modellist <- get(load(fn))

        # fit the models
        mlist <- CreateModel()

        # check
        for (hla.id in c("A", "B", "C", "DQA1", "DQB1", "DRB1"))
        {
            checkEquals(modellist[[hla.id]], modellist[[hla.id]],
                sprintf("Model HLA-%s", hla.id))
        }
    }

    invisible()
}
