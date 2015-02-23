# load the HIBAG package
library(HIBAG)


# a list of HLA genes
hla.list <- c("A", "B", "C", "DQA1", "DQB1", "DRB1")

# pre-defined lower bound of prediction accuracy
hla.acc <- c(0.9, 0.8, 0.8, 0.85, 0.85, 0.75)


for (hla.idx in seq_len(length(hla.list)))
{
	hla.id <- hla.list[hla.idx]

	# make a "hlaAlleleClass" object
	hla <- hlaAllele(HLA_Type_Table$sample.id,
		H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
		H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
		locus=hla.id, assembly="hg19")

	# divide HLA types randomly
	set.seed(100)
	hlatab <- hlaSplitAllele(hla, train.prop=0.5)

	# SNP predictors within the flanking region on each side
	region <- 500	# kb
	snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id,
		HapMap_CEU_Geno$snp.position,
		hla.id, region*1000, assembly="hg19")

	# training and validation genotypes
	train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
		snp.sel=match(snpid, HapMap_CEU_Geno$snp.id),
		samp.sel=match(hlatab$training$value$sample.id,
		HapMap_CEU_Geno$sample.id))
	test.geno <- hlaGenoSubset(HapMap_CEU_Geno,
		samp.sel=match(hlatab$validation$value$sample.id,
		HapMap_CEU_Geno$sample.id))


	# train a HIBAG model
	set.seed(100)
	model <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=10)
	summary(model)

	# validation
	pred <- predict(model, test.geno)
	summary(pred)

	# compare
	comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model,
		call.threshold=0)
	print(comp$overall)

	# check
	if (comp$overall$acc.haplo < hla.acc[hla.idx])
	{
		stop("HLA - ", hla.id, ", 'acc.haplo' should be >= ",
			hla.acc[hla.idx], ".")
	}

	cat("\n\n")
}
