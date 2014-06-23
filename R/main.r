#######################################################################
#
# Package Name: HIBAG v1.0
#
# Description:
#   HIBAG -- HLA Genotype Imputation with Attribute Bagging
#
# Author: Xiuwen Zheng
# Email: zhengx@u.washington.edu
#



#######################################################################
#
# the functions for SNP genotypes and haplotypes
#
#######################################################################

#
#
# hlaSNPGenoClass is a class of SNP genotypes
# list:
#     genotype -- a genotype matrix, ``# of SNPs'' X ``# of samples''
#     sample.id -- sample id
#     snp.id -- snp id
#     snp.position -- snp positions in basepair
#     snp.allele -- snp alleles, ``A allele/B allele''
#
#
# hlaSNPHaploClass is a class of SNP haplotypes
# list:
#     haplotype -- a haplotype matrix, ``# of SNPs'' X ``2 x # of samples''
#     sample.id -- sample id
#     snp.id -- snp id
#     snp.position -- snp positions in basepair
#     snp.allele -- snp alleles, ``A allele/B allele''
#
#


#######################################################################
# To create a "hlaSNPGenoClass" object (SNP genotype object)
#

hlaMakeSNPGeno <- function(genotype, sample.id, snp.id, snp.position, A.allele, B.allele)
{
	# check
	stopifnot(is.matrix(genotype))
	stopifnot(length(snp.id) == nrow(genotype))
	stopifnot(length(sample.id) == ncol(genotype))
	stopifnot(length(snp.id) == length(snp.position))
	stopifnot(length(snp.id) == length(A.allele))
	stopifnot(length(snp.id) == length(B.allele))
	stopifnot(is.character(A.allele))
	stopifnot(is.character(B.allele))

	rv <- list(genotype = genotype, sample.id = sample.id, snp.id = snp.id,
		snp.position = snp.position,
		snp.allele = paste(A.allele, B.allele, sep="/"))
	class(rv) <- "hlaSNPGenoClass"

    # valid SNP alleles
    flag <- !((A.allele %in% c("A", "T", "G", "C")) & (B.allele %in% c("A", "T", "G", "C")))
    if (any(flag))
    {
        warning(sprintf("There are %d SNPs with invalid alleles, and they have been removed.",
        	sum(flag)))
        rv <- hlaGenoSubset(rv, snp.sel=!flag)
    }
    # valid snp.id
    flag <- is.na(rv$snp.id)
    if (any(flag))
    {
        warning(sprintf("There are %d SNPs with missing SNP id, and they have been removed.",
        	sum(flag)))
        rv <- hlaGenoSubset(rv, snp.sel=!flag)
    }
    # valid snp.position
    flag <- is.na(rv$snp.position)
    if (any(flag))
    {
        warning(sprintf("There are %d SNPs with missing SNP positions, and they have been removed.",
        	sum(flag)))
        rv <- hlaGenoSubset(rv, snp.sel=!flag)
    }

	return(rv)
}


#######################################################################
# To create a "hlaSNPGenoClass" object (SNP genotype object)
#

hlaMakeSNPHaplo <- function(haplotype, sample.id, snp.id, snp.position, A.allele, B.allele)
{
	# check
	stopifnot(is.matrix(haplotype))
	stopifnot(length(snp.id) == nrow(haplotype))
	stopifnot(2*length(sample.id) == ncol(haplotype))
	stopifnot(length(snp.id) == length(snp.position))
	stopifnot(length(snp.id) == length(A.allele))
	stopifnot(length(snp.id) == length(B.allele))
	stopifnot(is.character(A.allele))
	stopifnot(is.character(B.allele))

	rv <- list(haplotype = haplotype, sample.id = sample.id, snp.id = snp.id,
		snp.position = snp.position,
		snp.allele = paste(A.allele, B.allele, sep="/"))
	class(rv) <- "hlaSNPHaploClass"

    # valid SNP alleles
    flag <- !((A.allele %in% c("A", "T", "G", "C")) & (B.allele %in% c("A", "T", "G", "C")))
    if (any(flag))
    {
        warning(sprintf("There are %d SNPs with invalid alleles, and they have been removed.",
        	sum(flag)))
        rv <- hlaHaploSubset(rv, snp.sel=!flag)
    }
    # valid snp.id
    flag <- is.na(rv$snp.id)
    if (any(flag))
    {
        warning(sprintf("There are %d SNPs with missing SNP id, and they have been removed.",
        	sum(flag)))
        rv <- hlaHaploSubset(rv, snp.sel=!flag)
    }
    # valid snp.position
    flag <- is.na(rv$snp.position)
    if (any(flag))
    {
        warning(sprintf("There are %d SNPs with missing SNP positions, and they have been removed.",
        	sum(flag)))
        rv <- hlaHaploSubset(rv, snp.sel=!flag)
    }

    return(rv)
}


#######################################################################
# To select a subset of SNP genotypes
#

hlaGenoSubset <- function(genoobj, samp.sel=NULL, snp.sel=NULL)
{
	# check
	stopifnot(class(genoobj) == "hlaSNPGenoClass")
	stopifnot(is.null(samp.sel) | is.logical(samp.sel) | is.integer(samp.sel))
	if (is.logical(samp.sel))
		stopifnot(length(samp.sel) == length(genoobj$sample.id))
	stopifnot(is.null(snp.sel) | is.logical(snp.sel) | is.integer(snp.sel))
	if (is.logical(snp.sel))
		stopifnot(length(snp.sel) == length(genoobj$snp.id))
	if (is.integer(samp.sel))
		stopifnot(length(unique(samp.sel)) == length(samp.sel))
	if (is.integer(snp.sel))
		stopifnot(length(unique(snp.sel)) == length(snp.sel))

	# subset
	if (is.null(samp.sel))
		samp.sel <- rep(TRUE, length(genoobj$sample.id))
	if (is.null(snp.sel))
		snp.sel <- rep(TRUE, length(genoobj$snp.id))
	rv <- list(genotype = genoobj$genotype[snp.sel, samp.sel],
		sample.id = genoobj$sample.id[samp.sel],
		snp.id = genoobj$snp.id[snp.sel],
		snp.position = genoobj$snp.position[snp.sel],
		snp.allele = genoobj$snp.allele[snp.sel]
	)
	class(rv) <- "hlaSNPGenoClass"
	return(rv)
}


#######################################################################
# To select a subset of SNP haplotypes
#

hlaHaploSubset <- function(haploobj, samp.sel=NULL, snp.sel=NULL)
{
	# check
	stopifnot(class(haploobj) == "hlaSNPHaploClass")
	stopifnot(is.null(samp.sel) | is.logical(samp.sel) | is.integer(samp.sel))
	if (is.logical(samp.sel))
		stopifnot(length(samp.sel) == length(haploobj$sample.id))
	stopifnot(is.null(snp.sel) | is.logical(snp.sel) | is.integer(snp.sel))
	if (is.logical(snp.sel))
		stopifnot(length(snp.sel) == length(haploobj$snp.id))
	if (is.integer(samp.sel))
		stopifnot(length(unique(samp.sel)) == length(samp.sel))
	if (is.integer(snp.sel))
		stopifnot(length(unique(snp.sel)) == length(snp.sel))

	# subset
	if (is.null(samp.sel))
		samp.sel <- rep(TRUE, length(haploobj$sample.id))
	if (is.numeric(samp.sel))
	{
		v <- samp.sel
		samp.sel <- rep(FALSE, length(haploobj$sample.id))
		samp.sel[v] <- TRUE
	}
	samp.sel <- rep(samp.sel, each=2)
	if (is.null(snp.sel))
		snp.sel <- rep(TRUE, length(haploobj$snp.id))

	rv <- list(haplotype = haploobj$haplotype[snp.sel, samp.sel],
		sample.id = haploobj$sample.id[samp.sel],
		snp.id = haploobj$snp.id[snp.sel],
		snp.position = haploobj$snp.position[snp.sel],
		snp.allele = haploobj$snp.allele[snp.sel]
	)
	class(rv) <- "hlaSNPHaploClass"
	return(rv)
}


#######################################################################
# To select a subset of SNP genotypes
#

hlaHaplo2Geno <- function(hapobj)
{
	stopifnot(class(hapobj) == "hlaSNPHaploClass")
	n <- dim(hapobj$haplotype)[2]
	rv <- list(
		genotype = hapobj$haplotype[, seq(1,n,2)] + hapobj$haplotype[, seq(2,n,2)],
		sample.id = hapobj$sample.id,
		snp.id = hapobj$snp.id,
		snp.position = hapobj$snp.position,
		snp.allele = hapobj$snp.allele
	)
	class(rv) <- "hlaSNPGenoClass"
	return(rv)
}


#######################################################################
# To get the overlapping SNPs between target and template with
#   corrected strand.
#

hlaGenoSwitchStrand <- function(target, template, verbose=TRUE)
{
	# check
	stopifnot(class(target) %in% c("hlaSNPGenoClass", "hlaSNPHaploClass"))
	stopifnot(class(template) %in% c("hlaSNPGenoClass", "hlaSNPHaploClass", "hlaAttrBagClass"))

	# initialize
	s1 <- hlaSNPID(template)
	s2 <- hlaSNPID(target)
	s <- intersect(s1, s2)
	if (length(s) <= 0) stop("There is no common SNP.")
	I1 <- match(s, s1); I2 <- match(s, s2)

	# compute allele frequencies
	if (class(template) == "hlaSNPGenoClass")
	{
		template.afreq <- rowMeans(template$genotype, na.rm=TRUE) * 0.5
	} else if (class(template) == "hlaSNPHaploClass") {
		template.afreq <- rowMeans(template$haplotype, na.rm=TRUE)
	} else {
		template.afreq <- template$snp.allele.freq
	}
	if (class(target) == "hlaSNPGenoClass")
	{
		target.afreq <- rowMeans(target$genotype, na.rm=TRUE) * 0.5
	} else {
		target.afreq <- rowMeans(target$haplotype, na.rm=TRUE)
	}

	# call
	gz <- .C("hlaAlleleStrand", template$snp.allele, template.afreq, I1,
		target$snp.allele, target.afreq, I2, length(s), out = logical(length(s)),
		err=integer(1), NAOK=TRUE, PACKAGE="HIBAG")
	if (gz$err != 0) stop(hlaErrMsg())
	if (verbose)
		cat(sprintf("The allele pairs of %d SNPs need to be switched.\n", sum(gz$out)))

	# result
	if (class(target) == "hlaSNPGenoClass")
	{
		geno <- target$genotype[I2, ]
		for (i in which(gz$out)) geno[i, ] <- 2 - geno[i, ]
		rv <- list(genotype = geno)
		rv$sample.id <- target$sample.id
		rv$snp.id <- target$snp.id[I2]
		rv$snp.position <- target$snp.position[I2]
		rv$snp.allele <- template$snp.allele[I1]
		class(rv) <- "hlaSNPGenoClass"
	} else {
		haplo <- target$haplotype[I2, ]
		for (i in which(gz$out)) haplo[i, ] <- 1 - haplo[i, ]
		rv <- list(haplotype = haplo)
		rv$sample.id <- target$sample.id
		rv$snp.id <- target$snp.id[I2]
		rv$snp.position <- target$snp.position[I2]
		rv$snp.allele <- template$snp.allele[I1]
		class(rv) <- "hlaSNPHaploClass"
	}

	return(rv)
}


#######################################################################
# To get the information of SNP ID and position
#

hlaSNPID <- function(obj)
{
	stopifnot(class(obj) %in% c("hlaSNPGenoClass", "hlaSNPHaploClass",
		"hlaAttrBagClass", "hlaAttrBagObj"))
	paste(obj$snp.id, obj$snp.position, sep="-")
}


#######################################################################
# To combine two SNP genotype dataset
#

hlaGenoCombine <- function(geno1, geno2, allele.check=TRUE)
{
    # check
    stopifnot(class(geno1) == "hlaSNPGenoClass")
    stopifnot(class(geno2) == "hlaSNPGenoClass")

	if (allele.check)
	{
	    tmp2 <- hlaGenoSwitchStrand(geno2, geno1)
    	tmp1 <- hlaGenoSubset(geno1, snp.sel=match(hlaSNPID(tmp2), hlaSNPID(geno1)))
    } else {
		s1 <- hlaSNPID(geno1); s2 <- hlaSNPID(geno2)
		set <- unique(intersect(s1, s2))
		tmp1 <- hlaGenoSubset(geno1, snp.sel=match(set, s1))
		tmp2 <- hlaGenoSubset(geno2, snp.sel=match(set, s2))
    }

    rv <- list(genotype = cbind(tmp1$genotype, tmp2$genotype),
        sample.id = c(tmp1$sample.id, tmp2$sample.id),
        snp.id = tmp1$snp.id, snp.position = tmp1$snp.position,
        snp.allele = tmp1$snp.allele)
    class(rv) <- "hlaSNPGenoClass"

    return(rv)
}


#######################################################################
# Convert to PLINK PED format
#

hlaGeno2PED <- function(geno, out.fn)
{
	# check
	stopifnot(class(geno) == "hlaSNPGenoClass")
	stopifnot(is.character(out.fn))

	# MAP file
	rv <- data.frame(chr=rep(6, length(geno$snp.id)), rs=geno$snp.id,
		morgan=rep(0, length(geno$snp.id)), bp=geno$snp.position, stringsAsFactors=FALSE)
	write.table(rv, file=paste(out.fn, ".map", sep=""),
		row.names=FALSE, col.names=FALSE, quote=FALSE)

	# PED file
	n <- length(geno$sample.id)
	m <- matrix("", nrow=n, ncol=2*length(geno$snp.id))
	for (i in 1:length(geno$snp.id))
	{
		allele <- unlist(strsplit(geno$snp.allele[i], "/"))
		g <- geno$genotype[i, ] + 1
		m[, 2*(i-1)+1] <- allele[c(2, 1, 1)[g]]
		m[, 2*(i-1)+2] <- allele[c(2, 2, 1)[g]]
	}
	rv <- cbind(Family=geno$sample.id, Ind=geno$sample.id,
		Paternal=rep(0, n), Maternal=rep(0, n), Sex=rep(0, n), Pheno=rep(-9, n),
		m)
	write.table(rv, file=paste(out.fn, ".ped", sep=""),
		row.names=FALSE, col.names=FALSE, quote=FALSE)

	# return
	return(invisible(NULL))
}


#######################################################################
# Convert from PLINK BED format
#

hlaBED2Geno <- function(bed.fn, fam.fn, bim.fn, rm.invalid.allele=TRUE,
	import.chr="xMHC", verbose=TRUE)
{
	# check
	stopifnot(is.character(bed.fn) & is.character(fam.fn) & is.character(bim.fn))
	stopifnot(is.logical(rm.invalid.allele))
	stopifnot(is.logical(verbose))

	# detect bed.fn
	bed <- .C("hlaBEDFlag", bed.fn, snporder=integer(1), err=integer(1),
		NAOK=TRUE, PACKAGE="HIBAG")
	if (bed$err != 0) stop(hlaErrMsg())
	if (verbose)
	{
		cat("\topen \"", bed.fn, sep="")
		if (bed$snporder == 0)
			cat("\" in the individual-major mode.\n")
		else
			cat("\" in the SNP-major mode.\n")
	}

	# read fam.fn
	famD <- read.table(fam.fn, header=FALSE, stringsAsFactors=FALSE)
	names(famD) <- c("FamilyID", "InvID", "PatID", "MatID", "Sex", "Pheno")
	if (length(unique(famD$InvID)) == dim(famD)[1])
	{
		sample.id <- famD$InvID
	} else {
		sample.id <- paste(famD$FamilyID, famD$InvID, sep="-")
		if (length(unique(sample.id)) != dim(famD)[1])
			stop("IDs in PLINK bed are not unique!")
	}
	if (verbose)
		cat("\topen \"", fam.fn, "\" DONE.\n", sep="")

	# read bim.fn
	bimD <- read.table(bim.fn, header=FALSE, stringsAsFactors=FALSE)
	names(bimD) <- c("chr", "snp.id", "map", "pos", "allele1", "allele2")

	# chromosome
	chr <- bimD$chr; chr[is.na(chr)] <- ""
	# position
	snp.pos <- bimD$pos
	snp.pos[!is.finite(snp.pos)] <- 0
	# snp.id
	snp.id <- bimD$snp.id
	# snp allele
	snp.allele <- paste(bimD$allele1, bimD$allele2, sep="/")
	if (verbose)
		cat("\topen \"", bim.fn, "\" DONE.\n", sep="")

	# create a hlaSNPGenoClass object
	if (import.chr == "xMHC")
	{
		snp.flag <- (chr==6) & (25759242<=snp.pos) & (snp.pos<=33534827)
		n.snp <- as.integer(sum(snp.flag))
		if (verbose)
			cat(sprintf("\tImport %d SNPs within the MHC region on chromosome 6.\n", n.snp))
	} else {
		snp.flag <- (chr %in% import.chr) & (snp.pos>0)
		n.snp <- as.integer(sum(snp.flag))
		if (verbose)
		{
			cat(sprintf("\tImport %d SNPs from chromosome %s.\n", n.snp,
				as.character(import.chr)))
		}
	}
	if (n.snp <= 0) stop("There is no SNP imported.")

	# call the C function
	rv <- .C("hlaConvBED", bed.fn, length(sample.id), length(snp.id), n.snp,
		(bed$snporder==0), snp.flag, verbose,
		geno = matrix(as.integer(0), nrow=n.snp, ncol=length(sample.id)),
		err=integer(1), NAOK=TRUE, PACKAGE="HIBAG")
	if (rv$err != 0) stop(hlaErrMsg())

	# result
	v <- list(genotype = rv$geno, sample.id = sample.id, snp.id = snp.id[snp.flag],
		snp.position = snp.pos[snp.flag], snp.allele = snp.allele[snp.flag])
	class(v) <- "hlaSNPGenoClass"

	# remove invalid snps
	if (rm.invalid.allele)
	{
		snp.allele <- v$snp.allele
		snp.allele[is.na(snp.allele)] <- "?/?"
		flag <- sapply(strsplit(snp.allele, "/"),
			function(x)
			{
				if (length(x) == 2)
				{
					all(x %in% c("A", "G", "C", "T"))
				} else {
					FALSE
				}
			}
		)
		if (any(!flag) & verbose)
			cat(sprintf("\t%d SNPs with invalid alleles have been removed.\n", sum(!flag)))

		# get a subset
		v <- hlaGenoSubset(v, snp.sel=flag)
	}

	return(v)
}






#######################################################################
# Summarize a "hlaSNPGenoClass" object
#

summary.hlaSNPGenoClass <- function(object, show=TRUE, ...)
{
	# check
	stopifnot(class(object) == "hlaSNPGenoClass")
	geno <- object

	fn <- function(x)
	{
		sprintf("min: %g, max: %g, mean: %g, median: %g, sd: %g",
			min(x, na.rm=TRUE), max(x, na.rm=TRUE),
			mean(x, na.rm=TRUE), median(x, na.rm=TRUE), sd(x, na.rm=TRUE))
	}

	rv <- list(mr.snp = hlaGenoMRate(geno), mr.samp = hlaGenoMRate_Samp(geno),
		maf = hlaGenoMFreq(geno),
		allele = table(geno$snp.allele))

	if (show)
	{
		cat("SNP Genotypes: \n")
		cat(sprintf("\t%d samples X %d SNPs\n",
			length(geno$sample.id), length(geno$snp.id)))
		cat(sprintf("\tSNPs range from %dbp to %dbp\n",
			min(geno$snp.position, na.rm=TRUE), max(geno$snp.position, na.rm=TRUE)))

		# missing rate for SNP
		cat(sprintf("Missing rate per SNP:\n\t%s\n", fn(rv$mr.snp)))
		# missing rate for sample
		cat(sprintf("Missing rate per sample:\n\t%s\n", fn(rv$mr.samp)))

		# minor allele frequency
		cat(sprintf("Minor allele frequency:\n\t%s\n", fn(rv$maf)))

		# allele information
		cat("Allele information:")
		print(rv$allele)
	}

	# return
	return(invisible(rv))
}


#######################################################################
# Summarize a "hlaSNPHaploClass" object
#

summary.hlaSNPHaploClass <- function(object, show=TRUE, ...)
{
	# check
	stopifnot(class(object) == "hlaSNPHaploClass")
	haplo <- object

	fn <- function(x)
	{
		sprintf("min: %g, max: %g, mean: %g, median: %g, sd: %g",
			min(x, na.rm=TRUE), max(x, na.rm=TRUE),
			mean(x, na.rm=TRUE), median(x, na.rm=TRUE), sd(x, na.rm=TRUE))
	}

	rv <- list(mr.snp = hlaGenoMRate(haplo), mr.samp = hlaGenoMRate_Samp(haplo),
		maf = hlaGenoMFreq(haplo),
		allele = table(haplo$snp.allele))

	if (show)
	{
		cat("SNP Haplotypes: \n")
		cat(sprintf("\t%d samples X %d SNPs\n",
			length(haplo$sample.id), length(haplo$snp.id)))
		cat(sprintf("\tSNPs range from %dbp to %dbp\n",
			min(haplo$snp.position, na.rm=TRUE), max(haplo$snp.position, na.rm=TRUE)))

		# missing rate for SNP
		cat(sprintf("Missing rate per SNP:\n\t%s\n", fn(rv$mr.snp)))
		# missing rate for sample
		cat(sprintf("Missing rate per sample:\n\t%s\n", fn(rv$mr.samp)))

		# minor allele frequency
		cat(sprintf("Minor allele frequency:\n\t%s\n", fn(rv$maf)))

		# allele information
		cat("Allele information:")
		print(rv$allele)
	}

	# return
	return(invisible(rv))
}




#######################################################################
#
# the function list for genotypes and haplotypes
#
#######################################################################

#######################################################################
# To the allele frequencies from genotypes or haplotypes
#

hlaGenoAFreq <- function(obj)
{
	# check
	stopifnot(class(obj) %in% c("hlaSNPGenoClass", "hlaSNPHaploClass"))
    if (class(obj) == "hlaSNPGenoClass")
    {
        rowMeans(obj$genotype, na.rm=TRUE) * 0.5
    } else {
        rowMeans(obj$haplotype, na.rm=TRUE)
    }
}


#######################################################################
# To the minor allele frequencies from genotypes or haplotypes
#

hlaGenoMFreq <- function(obj)
{
	# check
	stopifnot(class(obj) %in% c("hlaSNPGenoClass", "hlaSNPHaploClass"))
    if (class(obj) == "hlaSNPGenoClass")
    {
        F <- rowMeans(obj$genotype, na.rm=TRUE) * 0.5
    } else {
        F <- rowMeans(obj$haplotype, na.rm=TRUE)
    }
    pmin(F, 1-F)
}


#######################################################################
# To the missing rates from genotypes or haplotypes per SNP
#

hlaGenoMRate <- function(obj)
{
	# check
	stopifnot(class(obj) %in% c("hlaSNPGenoClass", "hlaSNPHaploClass"))
    if (class(obj) == "hlaSNPGenoClass")
    {
		rowMeans(is.na(obj$genotype))
	} else {
		rowMeans(is.na(obj$haplotype))
	}
}


#######################################################################
# To the missing rates from genotypes or haplotypes per sample
#

hlaGenoMRate_Samp <- function(obj)
{
	# check
	stopifnot(class(obj) %in% c("hlaSNPGenoClass", "hlaSNPHaploClass"))
    if (class(obj) == "hlaSNPGenoClass")
    {
        colMeans(is.na(obj$genotype))
    } else {
        F <- colMeans(is.na(obj$haplotype))
        F[seq(1, length(F), 2)] + F[seq(2, length(F), 2)]
    }
}






#######################################################################
#
# the function list for HLA types
#
#######################################################################


#######################################################################
# To get the starting and ending positions in basepair for HLA loci
#

hlaLociInfo <- function()
{
	ID <- c("A", "B", "C", "DRB1", "DRB5", "DQA1", "DQB1", "DPB1", "any")

	# starting position
	pos.HLA.start <- as.integer(c(30018310, 31429628, 31344508, 32654527, 32593129,
		32713161, 32735635, 33151738, NA))
	names(pos.HLA.start) <- ID

	# ending position
	pos.HLA.end <- as.integer(c(30021633, 31432914, 31347834, 32665559, 32605984,
		32719407, 32742444, 33162954, NA))
	names(pos.HLA.end) <- ID

	# return
	return(list(loci=ID, pos.HLA.start=pos.HLA.start, pos.HLA.end=pos.HLA.end))
}


#######################################################################
# Limit the resolution of HLA alleles
#

hlaAlleleDigit <- function(obj, max.resolution="4-digit")
{
	# check
	stopifnot((class(obj)=="hlaAlleleClass") | is.character(obj))
	if (is.character(obj)) stopifnot(is.vector(obj))
	stopifnot(max.resolution %in% c("2-digit", "4-digit", "6-digit", "8-digit",
		"allele", "protein", "2", "4", "6", "8", "full", ""))

	if (!(max.resolution %in% c("full", "")))
	{
		if (is.character(obj))
		{
			len <- c(1, 2, 3, 4, 1, 2, 1, 2, 3, 4)
			names(len) <- c("2-digit", "4-digit", "6-digit", "8-digit",
				"allele", "protein", "2", "4", "6", "8")
			maxlen <- len[[as.character(max.resolution)]]

			obj <- sapply(strsplit(obj, ":"), FUN =
					function(s, idx) {
						if (any(is.na(s)))
						{
							NA
						} else {
							if (length(idx) < length(s)) s <- s[idx]
							paste(s, collapse=":")
						}
					},
				idx = 1:maxlen)
		} else {
			rv <- list(locus = obj$locus,
				pos.start = obj$pos.start, pos.end = obj$pos.end,
				value = data.frame(sample.id = obj$value$sample.id,
					allele1 = hlaAlleleDigit(obj$value$allele1, max.resolution),
					allele2 = hlaAlleleDigit(obj$value$allele2, max.resolution),
					stringsAsFactors=FALSE)
			)
			if ("prob" %in% names(obj$value))
				rv$value$prob <- obj$value$prob
			class(rv) <- "hlaAlleleClass"
			obj <- rv
		}
	}

	return(obj)
}


#######################################################################
# To make a class of HLA alleles
#

hlaAllele <- function(sample.id, H1, H2, max.resolution="",
	locus="", locus.pos.start=NA, locus.pos.end=NA, prob=NULL, na.rm=TRUE)
{
	# check
	stopifnot(is.vector(sample.id))
	stopifnot(is.vector(H1) & is.character(H1))
	stopifnot(is.vector(H2) & is.character(H2))
	stopifnot(length(sample.id) == length(H1))
	stopifnot(length(sample.id) == length(H2))
	stopifnot(max.resolution %in% c("2-digit", "4-digit", "6-digit", "8-digit",
		"allele", "protein", "2", "4", "6", "8", "full", ""))

	HLAinfo <- hlaLociInfo()
	stopifnot(locus %in% HLAinfo$loci |
		(is.finite(locus.pos.start) & is.finite(locus.pos.end)))
	if (!is.null(prob))
		stopifnot(length(sample.id) == length(prob))

	# build
	H1[H1 == ""] <- NA; H1 <- hlaAlleleDigit(H1, max.resolution)
	H2[H2 == ""] <- NA; H2 <- hlaAlleleDigit(H2, max.resolution)
	if (locus %in% names(HLAinfo$pos.HLA.start))
	{
		if (!is.finite(locus.pos.start))
			locus.pos.start <- HLAinfo$pos.HLA.start[[locus]]
		if (!is.finite(locus.pos.end))
			locus.pos.end <- HLAinfo$pos.HLA.end[[locus]]
	}

	# remove missing values
	if (na.rm)
		flag <- (!is.na(H1)) & (!is.na(H2))
	else
		flag <- rep(TRUE, length(sample.id))

	# result
	rv <- list(locus = locus,
		pos.start = locus.pos.start, pos.end = locus.pos.end,
		value = data.frame(sample.id = sample.id[flag],
			allele1 = H1[flag], allele2 = H2[flag], stringsAsFactors=FALSE)
	)
	if (!is.null(prob))
		rv$value$prob <- prob[flag]
	class(rv) <- "hlaAlleleClass"
	return(rv)
}


#######################################################################
# To make a class of HLA alleles
#
# INPUT:
#   sample.id -- a vector of sample id
#   samp.sel -- a logical vector specifying selected samples
#

hlaAlleleSubset <- function(hla, samp.sel=NULL)
{
	# check
	stopifnot(class(hla) == "hlaAlleleClass")
	stopifnot(is.null(samp.sel) | is.logical(samp.sel) | is.integer(samp.sel))
	if (is.logical(samp.sel))
		stopifnot(length(samp.sel) == dim(hla$value)[1])
	if (is.integer(samp.sel))
		stopifnot(length(unique(samp.sel)) == length(samp.sel))

	# result
	if (is.null(samp.sel))
		samp.sel <- rep(TRUE, dim(hla$value)[1])
	rv <- list(locus = hla$locus,
		pos.start = hla$pos.start, pos.end = hla$pos.end,
		value = hla$value[samp.sel, ]
	)
	class(rv) <- "hlaAlleleClass"
	return(rv)
}


#######################################################################
# To combine two classes of HLA alleles
#
# INPUT:
#   H1 -- the first "hlaHLAAlleleClass" object
#   H2 -- the second "hlaHLAAlleleClass" object
#

hlaCombineAllele <- function(H1, H2)
{
	# check
	stopifnot(class(H1) == "hlaAlleleClass")
	stopifnot(class(H2) == "hlaAlleleClass")
	stopifnot(length(intersect(H1$sample.id, H2$sample.id)) == 0)
	stopifnot(H1$locus == H2$locus)
	stopifnot(H1$pos.start == H2$pos.start)
	stopifnot(H1$pos.end == H2$pos.end)

	id <- c("sample.id", "allele1", "allele2")

	# result
	rv <- list(locus = H1$locus,
		pos.start = H1$pos.start, pos.end = H1$pos.end,
		value = rbind(H1$value[, id], H2$value[, id])
	)
	if (!is.null(H1$value$prob) & !is.null(H2$value$prob))
	{
		rv$value$prob <- c(H1$value$prob, H2$value$prob)
	}
	class(rv) <- "hlaAlleleClass"
	return(rv)
}


#######################################################################
# To compare HLA alleles
#

hlaCompareAllele <- function(TrueHLA, PredHLA, allele.limit=NULL,
	call.threshold=NaN, max.resolution="", output.individual=FALSE, verbose=TRUE)
{
	# check
	stopifnot(class(TrueHLA) == "hlaAlleleClass")
	stopifnot(class(PredHLA) == "hlaAlleleClass")
	stopifnot(is.null(allele.limit) | is.vector(allele.limit) |
		class(allele.limit)=="hlaAttrBagClass" | class(allele.limit)=="hlaAttrBagObj")
	stopifnot(max.resolution %in% c("2-digit", "4-digit", "6-digit", "8-digit",
		"allele", "protein", "2", "4", "6", "8", "full", ""))
	stopifnot(is.logical(output.individual))
	stopifnot(is.logical(verbose))

	# get the common samples
	samp <- intersect(TrueHLA$value$sample.id, PredHLA$value$sample.id)
	if ((length(samp) != length(TrueHLA$value$sample.id)) |
		(length(samp) != length(PredHLA$value$sample.id)))
	{
		if (verbose)
		{
			cat("Call hlaCompareAllele: there are", length(samp),
				"samples.\n")
		}
	}
	# True HLA
	flag <- match(samp, TrueHLA$value$sample.id)
	if (length(samp) != length(TrueHLA$value$sample.id))
	{
		TrueHLA <- hlaAlleleSubset(TrueHLA, flag)
	} else {
		if (!all(flag == 1:length(TrueHLA$value$sample.id)))
			TrueHLA <- hlaAlleleSubset(TrueHLA, flag)
	}
	# Predicted HLA
	flag <- match(samp, PredHLA$value$sample.id)
	if (length(samp) != length(PredHLA$value$sample.id))
	{
		PredHLA <- hlaAlleleSubset(PredHLA, flag)
	} else {
		if (!all(flag == 1:length(PredHLA$value$sample.id)))
			PredHLA <- hlaAlleleSubset(PredHLA, flag)
	}

	# init
	flag <- !is.na(TrueHLA$value$allele1) & !is.na(TrueHLA$value$allele2) &
		!is.na(PredHLA$value$allele1) & !is.na(PredHLA$value$allele2)
	ts1 <- TrueHLA$value$allele1[flag]; ts2 <- TrueHLA$value$allele2[flag]
	ps1 <- PredHLA$value$allele1[flag]; ps2 <- PredHLA$value$allele2[flag]
	samp.id <- TrueHLA$value$sample.id[flag]

	# call threshold
	if (is.finite(call.threshold))
	{
		prob <- PredHLA$value$prob
		if (!is.null(prob)) prob <- prob[flag]
	} else {
		prob <- NULL
	}

	# allele limitation
	if (!is.null(allele.limit))
	{
		if (class(allele.limit) %in% c("hlaAttrBagClass", "hlaAttrBagObj"))
		{
			allele <- levels(factor(allele.limit$hla.allele))
			TrainFreq <- allele.limit$hla.freq
		} else {
			allele <- levels(factor(allele.limit))
			TrainFreq <- NULL
		}
	} else {
		allele <- levels(factor(c(ts1, ts2)))
		TrainFreq <- NULL
	}

	# max resolution
	if (!(max.resolution %in% c("full", "")))
	{
		ts1 <- hlaAlleleDigit(ts1, max.resolution)
		ts2 <- hlaAlleleDigit(ts2, max.resolution)
		ps1 <- hlaAlleleDigit(ps1, max.resolution)
		ps2 <- hlaAlleleDigit(ps2, max.resolution)
		tmp <- hlaAlleleDigit(allele, max.resolution)
		allele <- levels(factor(tmp))
		if ((length(tmp) != length(allele)) & !is.null(TrainFreq))
		{
			x <- rep(0, length(allele))
			for (i in 1:length(allele))
				x[i] <- sum(TrainFreq[tmp == allele[i]])
			TrainFreq <- x
		}
	}

	# allele filter
	flag <- (ts1 %in% allele) & (ts2 %in% allele)
	ts1 <- ts1[flag]; ts2 <- ts2[flag]
	ps1 <- ps1[flag]; ps2 <- ps2[flag]
	samp.id <- samp.id[flag]
	if (!is.null(prob)) prob <- prob[flag]

	# init ...
	cnt.ind <- 0; cnt.haplo <- 0; cnt.call <- 0
	n <- length(ts1)
	m <- length(allele)

	TrueNum <- rep(0, m); names(TrueNum) <- allele
	TrueNumAll <- rep(0, m); names(TrueNumAll) <- allele
	PredNum <- rep(0, m+1); names(PredNum) <- c(allele, "...")
	confusion <- matrix(0.0, nrow = m+1, ncol = m,
		dimnames = list(Predict=names(PredNum), True=names(TrueNum)))
	WrongTab <- NULL

	# for PredNum
	fn <- function(x, LT)
		{ if (x %in% LT) { return(x) } else { return("...") } }

	acc.array <- rep(NaN, n)
	if (n > 0)
	{
		for (i in 1:n)
		{
			# increase
			TrueNumAll[[ ts1[i] ]] <- TrueNumAll[[ ts1[i] ]] + 1
			TrueNumAll[[ ts2[i] ]] <- TrueNumAll[[ ts2[i] ]] + 1

			# probability cut-off
			if (is.null(prob))
				flag <- TRUE
			else
				flag <- prob[i]>=call.threshold
			if (flag)
			{
				# update TrueNum and PredNum
				TrueNum[[ ts1[i] ]] <- TrueNum[[ ts1[i] ]] + 1
				TrueNum[[ ts2[i] ]] <- TrueNum[[ ts2[i] ]] + 1
				PredNum[[ fn(ps1[i], allele) ]] <- PredNum[[ fn(ps1[i], allele) ]] + 1
				PredNum[[ fn(ps2[i], allele) ]] <- PredNum[[ fn(ps2[i], allele) ]] + 1

				# correct count of individuals
				if ( ((ts1[i]==ps1[i]) & (ts2[i]==ps2[i])) | ((ts2[i]==ps1[i]) & (ts1[i]==ps2[i])) )
				{
					cnt.ind <- cnt.ind + 1
				}

				# correct count of haplotypes
				s <- c(ts1[i], ts2[i]); p <- c(ps1[i], ps2[i])
				hnum <- 0
				if ((s[1]==p[1]) | (s[1]==p[2]))
				{
					if (s[1]==p[1]) { p[1] <- "" } else { p[2] <- "" }
					confusion[s[1], s[1]] <- confusion[s[1], s[1]] + 1
					cnt.haplo <- cnt.haplo + 1
					hnum <- hnum + 1
				}
				if ((s[2]==p[1]) | (s[2]==p[2]))
				{
					confusion[s[2], s[2]] <- confusion[s[2], s[2]] + 1
					cnt.haplo <- cnt.haplo + 1
					hnum <- hnum + 1
				}
				acc.array[i] <- 0.5*hnum

				# for confusion matrix
				s <- c(ts1[i], ts2[i]); p <- c(ps1[i], ps2[i])
				if (hnum == 1)
				{
					if ((s[1]==p[1]) | (s[1]==p[2]))
					{
						if (s[1]==p[1])
						{
							confusion[fn(p[2], allele), s[2]] <- confusion[fn(p[2], allele), s[2]] + 1
						} else {
							confusion[fn(p[1], allele), s[2]] <- confusion[fn(p[1], allele), s[2]] + 1
						}
					} else {
						if (s[2]==p[1])
						{
							confusion[fn(p[2], allele), s[1]] <- confusion[fn(p[2], allele), s[1]] + 1
						} else {
							confusion[fn(p[1], allele), s[1]] <- confusion[fn(p[1], allele), s[1]] + 1
						}
					}
				} else if (hnum == 0)
				{
					WrongTab <- cbind(WrongTab, c(s, fn(p[1], allele), fn(p[2], allele)))
				}

				# the number of calling
				cnt.call <- cnt.call + 1
			}
		}
	}

	# overall
	overall <- data.frame(total.num.ind = n,
		crt.num.ind = cnt.ind, crt.num.haplo = cnt.haplo,
		acc.ind = cnt.ind/cnt.call, acc.haplo = 0.5*cnt.haplo/cnt.call,
		call.threshold = call.threshold)
	if (is.finite(call.threshold))
	{
		overall$n.call <- cnt.call
		overall$call.rate <- cnt.call / n
	} else {
		overall$n.call <- n
		overall$call.rate <- 1.0
		overall$call.threshold <- 0
	}

	# confusion matrix
	if (is.null(WrongTab))
	{
		nw <- as.integer(0)
	} else {
		nw <- ncol(WrongTab)
	}
	rv <- .C("hlaConfusion", as.integer(m), confusion,
		nw, match(WrongTab, names(PredNum)) - as.integer(1),
		out = matrix(0.0, nrow=m+1, ncol=m, dimnames=list(Predict=names(PredNum), True=names(TrueNum))),
		tmp = double((m+1)*m),
		err = integer(1), NAOK = TRUE, PACKAGE = "HIBAG")
	if (rv$err != 0) stop(hlaErrMsg())
	confusion <- round(rv$out, 2)

	# detail -- sensitivity and specificity
	detail <- data.frame(allele = allele)
	if (!is.null(TrainFreq))
	{
		detail$train.num <- 2 * TrainFreq * allele.limit$n.samp
		detail$train.freq <- TrainFreq
	}
	detail$valid.num <- TrueNumAll
	detail$valid.freq <- TrueNumAll / sum(TrueNumAll)
	detail$call.rate <- TrueNum / TrueNumAll

	sens <- diag(confusion) / TrueNum
	spec <- 1 - (PredNum[1:m] - diag(confusion)) / (2*cnt.call - TrueNum)
	detail$accuracy <- (sens*TrueNum + spec*(2*cnt.call - TrueNum)) / (2*cnt.call)
	detail$sensitivity <- sens
	detail$specificity <- spec
	detail$ppv <- diag(confusion) / rowSums(confusion)[1:m]
	detail$npv <- 1 - (TrueNum - diag(confusion)) / (2*n - rowSums(confusion)[1:m])

	detail$call.rate[!is.finite(detail$call.rate)] <- 0
	detail[detail$call.rate<=0, c("sensitivity", "specificity", "ppv", "npv", "accuracy")] <- NaN

	# get miscall
	rv <- confusion; diag(rv) <- 0
	m.max <- apply(rv, 2, max); m.idx <- apply(rv, 2, which.max)
	s <- names(PredNum)[m.idx]; s[m.max<=0] <- NA
	p <- m.max / apply(rv, 2, sum)
	detail <- cbind(detail, miscall=s, miscall.prop=p, stringsAsFactors=FALSE)

	# output
	rv <- list(overall=overall, confusion=confusion, detail=detail)
	if (output.individual)
		rv$individual <- data.frame(sample.id=samp.id, accuracy=acc.array, stringsAsFactors=FALSE)
	return(rv)
}


#######################################################################
# Return a sample list satisfying the filter conditions
#

hlaSampleAllele <- function(TrueHLA, allele.limit = NULL, max.resolution="")
{
	# check
	stopifnot(class(TrueHLA) == "hlaAlleleClass")
	stopifnot(is.null(allele.limit) | is.vector(allele.limit) |
		class(allele.limit)=="hlaAttrBagClass" | class(allele.limit)=="hlaAttrBagObj")
	stopifnot(max.resolution %in% c("2-digit", "4-digit", "6-digit", "8-digit",
		"allele", "protein", "2", "4", "6", "8", "full", ""))

	# init
	flag <- !is.na(TrueHLA$value$allele1) & !is.na(TrueHLA$value$allele2)
	ts1 <- TrueHLA$value$allele1[flag]
	ts2 <- TrueHLA$value$allele2[flag]

	# max resolution
	if (!(max.resolution %in% c("full", "")))
	{
		ts1 <- hlaAlleleDigit(ts1, max.resolution)
		ts2 <- hlaAlleleDigit(ts2, max.resolution)
	}

	# allele limitation
	if (!is.null(allele.limit))
	{
		if (class(allele.limit) %in% c("hlaAttrBagClass", "hlaAttrBagObj"))
		{
			allele <- levels(factor(allele.limit$hla.allele))
		} else {
			allele <- levels(factor(allele.limit))
		}
		if (!(max.resolution %in% c("full", "")))
		{
			allele <- hlaAlleleDigit(allele, max.resolution)
		}
		flag[flag] <- (ts1 %in% allele) & (ts2 %in% allele)
	}

	# return
	return(TrueHLA$value$sample.id[flag])
}


#######################################################################
# Divide the list of HLA types to the training and validation sets
#
# INPUT:
#   HLA -- the HLA types, a "hlaHLAAlleleClass" object
#   train.prop -- the proportion of training samples
#

hlaSplitAllele <- function(HLA, train.prop=0.5)
{
	# check
	stopifnot(class(HLA) == "hlaAlleleClass")

	train.set <- NULL
	H <- HLA
	while (dim(H$value)[1] > 0)
	{
		v <- summary(H, show=FALSE)
		if (dim(v)[1] > 1)
		{
			v <- v[order(v[, "count"]), ]
		}

		allele <- rownames(v)[1]
		samp.id <- H$value$sample.id[H$value$allele1==allele | H$value$allele2==allele]
		samp.id <- as.character(samp.id)

		n.train <- ceiling(length(samp.id) * train.prop)
		train.sampid <- sample(samp.id, n.train)
		train.set <- c(train.set, train.sampid)

		H <- hlaAlleleSubset(H, samp.sel=
			match(setdiff(H$value$sample.id, samp.id), H$value$sample.id))
	}

	train.set <- train.set[order(train.set)]
	list(
		training = hlaAlleleSubset(HLA, samp.sel=match(train.set, HLA$value$sample.id)),
		validation = hlaAlleleSubset(HLA, samp.sel=match(setdiff(HLA$value$sample.id,
			train.set), HLA$value$sample.id))
	)
}


#######################################################################
# To select SNPs in the flanking region of a specified HLA locus
#
# INPUT:
#   snp.id -- a vector of snp id
#   position -- a vector of positions
#   hla.id -- the name of HLA locus
#   flank.bp -- the distance in basepair
#

hlaFlankingSNP <- function(snp.id, position, hla.id, flank.bp=500*1000)
{
	# init
	HLAInfo <- hlaLociInfo()

	# check
	stopifnot(length(snp.id) == length(position))
	stopifnot(is.character(hla.id))
	stopifnot(length(hla.id) == 1)
	stopifnot(hla.id %in% HLAInfo$loci)

	pos.start <- HLAInfo$pos.HLA.start[hla.id] - flank.bp
	pos.end <- HLAInfo$pos.HLA.end[hla.id] + flank.bp
	flag <- (pos.start <= position) & (position <= pos.end)
	return(snp.id[flag])
}


#######################################################################
# Summary a "hlaAlleleClass" object
#
# INPUT:
#   hla -- a "hlaAlleleClass" object
#

summary.hlaAlleleClass <- function(object, show=TRUE, ...)
{
	# check
	stopifnot(class(object) == "hlaAlleleClass")
	hla <- object

	HLA <- c(hla$value$allele1, hla$value$allele2)
	count <- table(HLA)
	freq <- prop.table(count)
	rv <- cbind(count=count, freq=freq)

	if (show)
	{
		cat(sprintf("HLA-%s, range [%dbp, %dbp]\n", hla$locus, hla$pos.start, hla$pos.end))
		cat(sprintf("There are %d samples.\n", dim(hla$value)[1]))
		cat(sprintf("There are %d different HLA alleles.\n", length(count)))
	}

	# return
	return(rv)
}




##########################################################################
##########################################################################
#
# Attribute Bagging (AB) method
#

##########################################################################
# To fit an attribute bagging model for predicting
#

hlaAttrBagging <- function(hla, genotype, nclassifier=100, mtry=c("sqrt", "all", "one"),
	prune=TRUE, rm.na=TRUE, verbose=TRUE, verbose.detail=FALSE)
{
	# check
	stopifnot(class(hla) == "hlaAlleleClass")
	stopifnot(class(genotype) == "hlaSNPGenoClass")
	stopifnot(is.character(mtry) | is.numeric(mtry))
	stopifnot(is.logical(verbose))
	stopifnot(is.logical(verbose.detail))
	if (verbose.detail) verbose <- TRUE

	# get the common samples
	samp.id <- intersect(hla$value$sample.id, genotype$sample.id)

	# hla types
	samp.flag <- match(samp.id, hla$value$sample.id)
	hla.allele1 <- hla$value$allele1[samp.flag]
	hla.allele2 <- hla$value$allele2[samp.flag]
	if (rm.na)
	{
		if (any(is.na(c(hla.allele1, hla.allele2))))
		{
			warning("There are missing HLA alleles, and the corresponding samples have been removed.")
			flag <- is.na(hla.allele1) | is.na(hla.allele2)
			samp.id <- setdiff(samp.id, hla$value$sample.id[samp.flag[flag]])
			samp.flag <- match(samp.id, hla$value$sample.id)
			hla.allele1 <- hla$value$allele1[samp.flag]
			hla.allele2 <- hla$value$allele2[samp.flag]
		}
	} else {
		if (any(is.na(c(hla.allele1, hla.allele2))))
		{
			stop("There are missing HLA alleles!")
		}
	}

	# SNP genotypes
	samp.flag <- match(samp.id, genotype$sample.id)
	snp.geno <- genotype$genotype[, samp.flag]
	storage.mode(snp.geno) <- "integer"

	tmp.snp.id <- genotype$snp.id
	tmp.snp.position <- genotype$snp.position
	tmp.snp.allele <- genotype$snp.allele

	# remove mono-SNPs
	snpsel <- rowMeans(snp.geno, na.rm=TRUE)
	snpsel[!is.finite(snpsel)] <- 0
	snpsel <- (0 < snpsel) & (snpsel < 2)
	if (sum(!snpsel) > 0)
	{
		snp.geno <- snp.geno[snpsel, ]
		if (verbose)
			cat(sprintf("%s monomorphic SNPs have been removed.\n", sum(!snpsel)))
		tmp.snp.id <- tmp.snp.id[snpsel]
		tmp.snp.position <- tmp.snp.position[snpsel]
		tmp.snp.allele <- tmp.snp.allele[snpsel]
	}

	# sampling probabilites of variable selection, it will be implemented in future
	prior.cutoff <- 0
	# possibly it is better to use the flat prior probability
	var.prior.prob <- "flat"

	if (var.prior.prob == "ld")
	{
		prior <- hlaGenoLD(hlaAlleleSubset(hla, match(samp.id, hla$value$sample.id)), snp.geno)
	} else if (var.prior.prob == "flat") {
		prior <- NULL
	} else {
		stop("Error `var.prior.prob'!")
	}
	if (!is.null(prior))
	{
		prior[!is.finite(prior)] <- 0
		prior <- prior / sum(prior)
		snpsel <- (prior > prior.cutoff)
		if (sum(!snpsel) > 0)
		{
			snp.geno <- snp.geno[snpsel, ]
			if (verbose)
				cat(sprintf("%s SNPs with low prior prob have been removed.\n", sum(!snpsel)))
			tmp.snp.id <- tmp.snp.id[snpsel]
			tmp.snp.position <- tmp.snp.position[snpsel]
			tmp.snp.allele <- tmp.snp.allele[snpsel]
			prior <- prior[snpsel]
			prior <- prior / sum(prior)
		}
		if (verbose)
		{
			cat("Sampling variables with prior probability:\n")
			print(summary(prior))
		}
	} else {
		if (verbose) cat("Sampling variables with flat prior probability.\n")
	}


	if (length(samp.id) <= 0)
		stop("There is no common sample between `hla' and `genotype'.")
	if (length(dim(snp.geno)[1]) <= 0)
		stop("There is no valid SNP markers.")


	###################################################################
	# initialize ...

	n.snp <- dim(snp.geno)[1]      # Num. of SNPs
	n.samp <- dim(snp.geno)[2]     # Num. of samples
	H <- factor(c(hla.allele1, hla.allele2))
	n.hla <- nlevels(H)
	H1 <- as.integer(H[1:n.samp]) - as.integer(1)
	H2 <- as.integer(H[(n.samp+1):(2*n.samp)]) - as.integer(1)

	# create an attribute bagging object
	rv <- .C("hlaAB_Model_Training", n.snp, n.samp, snp.geno, n.hla,
		H1, H2, AB=integer(1), err=integer(1),
		NAOK = TRUE, PACKAGE = "HIBAG")
	if (rv$err != 0) stop(hlaErrMsg())
	ABmodel <- rv$AB

	# number of variables randomly sampled as candidates at each split
	mtry <- mtry[1]
	if (is.character(mtry))
	{
		if (mtry == "sqrt")
		{
			mtry <- ceiling(sqrt(n.snp))
		} else if (mtry == "all")
		{
			mtry <- n.snp
		} else if (mtry == "one")
		{
			mtry <- as.integer(1)
		} else {
			stop("Invalid mtry!")
		}
	} else if (is.numeric(mtry))
	{
		if (is.finite(mtry))
		{
			if ((0 < mtry) & (mtry < 1)) mtry <- n.snp*mtry
			mtry <- ceiling(mtry)
			if (mtry > n.snp) mtry <- n.snp
		} else {
			mtry <- ceiling(sqrt(n.snp))
		}
	} else {
		stop("Invalid mtry value!")
	}
	if (mtry <= 0) mtry <- as.integer(1)

	if (verbose)
	{
		cat("Build an attribute bagging model with", nclassifier,
			"individual classifiers, and each classifier is created from", mtry,
			"randomly selected SNPs.\n")
		cat("\t# of SNPs:", n.snp, ", # of samples:", n.samp, "\n")
		cat("\t# of HLA alleles:", n.hla, "\n")
	}


	###################################################################
	# training ...
	# add new individual classifers
	rv <- .C("hlaAB_NewClassifiers", ABmodel, as.integer(nclassifier),
		as.integer(mtry), as.double(prior), !is.null(prior),
		as.logical(prune), verbose, verbose.detail, debug=FALSE,
		err=integer(1), NAOK = TRUE, PACKAGE = "HIBAG")
	if (rv$err != 0) stop(hlaErrMsg())

	# output
	rv <- list(n.samp = n.samp, n.snp = n.snp, sample.id = samp.id,
		snp.id = tmp.snp.id, snp.position = tmp.snp.position,
		snp.allele = tmp.snp.allele,
		snp.allele.freq = 0.5*rowMeans(snp.geno, na.rm=TRUE),
		hla.locus = hla$locus, hla.allele = levels(H), hla.freq = prop.table(table(H)),
		model = ABmodel)
	if (!is.null(prior)) rv$var.prior.prob <- prior
	class(rv) <- "hlaAttrBagClass"
	return(rv)
}


##########################################################################
# To fit an attribute bagging model for predicting
#

hlaClusterAttrBagging <- function(cl, hla, genotype, nclassifier=100,
	mtry=c("sqrt", "all", "one"), prune=TRUE, rm.na=TRUE,
	verbose=TRUE, verbose.detail=FALSE)
{
    if (!require(parallel))
    {
        if (!require(snow)) 
            stop("the `parallel' or `snow' package is needed.")
    }

	# check
    stopifnot(inherits(cl, "cluster"))
	stopifnot(class(hla) == "hlaAlleleClass")
	stopifnot(class(genotype) == "hlaSNPGenoClass")
	stopifnot(is.character(mtry) | is.numeric(mtry))
	stopifnot(is.logical(verbose))
	stopifnot(is.logical(verbose.detail))
	if (verbose.detail) verbose <- TRUE

	# split jobs
	clseq <- splitIndices(nclassifier, length(cl))
	clseqlen <- sapply(clseq, "length")

	if (verbose)
	{
		if (length(cl) <= 1)
		{
			cat("There is one job.\n")
		} else {
			cat("There are", length(cl), "jobs, and the numbers of classifiers created by jobs are:\n")
			cat(clseqlen)
			cat("\n")
		}
	}

	# set random number
	clusterSetRNGStream(cl)

	# run jobs
	ans <- clusterApply(cl, clseqlen, fun =
		function(nclass, hla, genotype, mtry, prune, rm.na, verbose, verbose.detail)
		{
			if (nclass > 0)
			{
				library(HIBAG)
				model <- hlaAttrBagging(hla=hla, genotype=genotype, nclassifier=nclass,
					mtry=mtry, prune=prune, rm.na=rm.na,
					verbose=verbose, verbose.detail=verbose.detail)
				rv <- hlaModelToObj(model)
				hlaClose(model)
			} else {
				rv <- NULL
			}
			return(rv)
		}, hla=hla, genotype=genotype, mtry=mtry, prune=prune, rm.na=rm.na,
			verbose=verbose, verbose.detail=verbose.detail)

	# combine
	ML <- NULL
	for (i in 1:length(ans))
	{
		v <- ans[[i]]
		if (!is.null(v))
		{
			if (is.null(ML))
			{
				ML <- v
			} else {
				ML <- hlaCombineModelObj(ML, v)
			}
		}
	}

	return(hlaModelfromObj(ML))
}


##########################################################################
# To fit an attribute bagging model for predicting
#

hlaClose <- function(model)
{
	# check
	stopifnot(class(model) == "hlaAttrBagClass")

	# class handler
	rv <- .C("hlaAB_Close", model$model, err=integer(1), NAOK = TRUE, PACKAGE = "HIBAG")
	if (rv$err != 0) stop(hlaErrMsg())

	# output
	return(invisible(NULL))
}


#######################################################################
# Check missing SNP predictors
#

hlaCheckSNPs <- function(model, object, verbose=TRUE)
{
	# check
	stopifnot(class(model) %in% c("hlaAttrBagClass", "hlaAttrBagObj"))
	stopifnot((is.vector(object) & is.character(object)) |
		(class(object) == "hlaSNPGenoClass"))

	# initialize
	if (class(model) == "hlaAttrBagClass")
		model <- hlaModelToObj(model)

	# show information
	cat("The HIBAG model:\n")
	cat(sprintf("\tThere are %d SNP predictors in total.\n", length(model$snp.id)))
	cat(sprintf("\tThere are %d individual classifiers.\n", length(model$classifiers)))

	if (is.vector(object))
	{
		target.snp <- as.character(object)
		src.snp <- model$snp.id
	} else {
		target.snp <- hlaSNPID(object)
		src.snp <- hlaSNPID(model)
	}

	NumOfSNP <- integer(length(model$classifiers))
	NumOfValidSNP <- integer(length(model$classifiers))

	# enumerate each classifier
	for (i in 1:length(model$classifiers))
	{
		v <- model$classifiers[[i]]
		flag <- src.snp[v$snpidx] %in% target.snp
		NumOfSNP[i] <- length(v$snpidx)
		NumOfValidSNP[i] <- sum(flag)
	}

	rv <- data.frame(NumOfValidSNP = NumOfValidSNP, NumOfSNP = NumOfSNP,
		fraction = NumOfValidSNP/NumOfSNP)

	if (verbose)
	{
		cat("Summarize the missing fractions of SNP predictors per classifier:\n")
		print(summary(1 - rv$fraction))
	}

	# output
	invisible(rv)
}


#######################################################################
# Predict HLA types from unphased SNP data
#

predict.hlaAttrBagClass <- function(object, genotypes, type=c("response", "prob"),
	allele.check=TRUE, verbose=TRUE, ...)
{
	# check
	stopifnot(class(object) == "hlaAttrBagClass")
	stopifnot(type[1] %in% c("response", "prob"))
	stopifnot(is.logical(allele.check))

	if (verbose)
		cat(sprintf("There are %d SNP predictors in the model.\n", length(object$snp.id)))

	if (class(genotypes) != "hlaSNPGenoClass")
	{
		stopifnot(is.numeric(genotypes))
		stopifnot(is.vector(genotypes) | is.matrix(genotypes))
		if (is.vector(genotypes))
		{
			stopifnot(length(genotypes) == object$n.snp)
			if (!is.null(names(genotypes)))
				stopifnot(all(names(genotypes) == object$snp.id))
			genotypes <- matrix(genotypes, ncol=1)
		} else {
			stopifnot(nrow(genotypes) == object$n.snp)
			if (!is.null(rownames(genotypes)))
				stopifnot(all(rownames(genotypes) == object$snp.id))
		}
		geno.sampid <- 1:ncol(genotypes)
	} else {
		geno.sampid <- genotypes$sample.id
		obj.id <- hlaSNPID(object)
		geno.id <- hlaSNPID(genotypes)

		# flag boolean
		flag <- FALSE
		if (length(obj.id) == length(geno.id))
		{
			if (all(obj.id == geno.id))
				flag <- TRUE
		}

		# check and switch A/B alleles
		if (flag)
		{
			if (allele.check)
			{
				genotypes <- hlaGenoSwitchStrand(genotypes, object, verbose)$genotype
			} else {
				genotypes <- genotypes$genotype
			}
		} else {
			# snp selection
			snp.sel <- match(obj.id, geno.id)
			# tmp variable
			tmp <- list(genotype = genotypes$genotype[snp.sel,],
				sample.id = genotypes$sample.id,
				snp.id = object$snp.id, snp.position = object$snp.position,
				snp.allele = genotypes$snp.allele[snp.sel])
			flag <- is.na(tmp$snp.allele)
			tmp$snp.allele[flag] <- object$snp.allele[match(tmp$snp.id[flag], object$snp.id)]
			class(tmp) <- "hlaSNPGenoClass"
			# verbose
			if (verbose)
			{
				cnt <- sum(flag)
				if (cnt > 0)
				{
					cat(sprintf("The test genotypes do not have %d SNP predictors.\n", cnt))
				}
			}
			# switch
			if (allele.check)
			{
				genotypes <- hlaGenoSwitchStrand(tmp, object, verbose)$genotype
			} else {
				genotypes <- genotypes$genotype
			}
		}
	}

	# initialize ...
	n.samp <- dim(genotypes)[2]
	n.hla <- length(object$hla.allele)

	# to predict HLA types
	if (type[1] == "response")
	{
		rv <- .C("hlaAB_Predict", object$model, as.integer(genotypes), n.samp,
			as.logical(verbose), H1=integer(n.samp), H2=integer(n.samp),
			prob=double(n.samp), err=integer(1), NAOK=TRUE, PACKAGE="HIBAG")
		if (rv$err != 0) stop(hlaErrMsg())

		res <- hlaAllele(geno.sampid,
			object$hla.allele[rv$H1 + 1], object$hla.allele[rv$H2 + 1],
			locus=object$hla.locus, prob=rv$prob)
	} else {
		rv <- .C("hlaAB_Predict_Prob", object$model, as.integer(genotypes), n.samp,
			as.logical(verbose), prob=matrix(NaN, nrow=n.hla*(n.hla+1)/2, ncol=n.samp),
			err=integer(1), NAOK=TRUE, PACKAGE="HIBAG")
		if (rv$err != 0) stop(hlaErrMsg())

		res <- rv$prob
		colnames(res) <- geno.sampid
		m <- outer(object$hla.allele, object$hla.allele, function(x, y) paste(x, y, sep="."))
		rownames(res) <- m[lower.tri(m, diag=TRUE)]
	}

	# return
	return(res)
}


#######################################################################
# summarize the "hlaAttrBagClass" object
#

summary.hlaAttrBagClass <- function(object, show=TRUE, ...)
{
	obj <- hlaModelToObj(object)
	summary(obj, show=show)
}


#######################################################################
# Save the parameters in a model of attribute bagging
#

hlaModelToObj <- function(model)
{
	# check
	stopifnot(class(model)=="hlaAttrBagClass")

	# call, get the number of classifiers
	rv <- .C("hlaAB_GetNumClassifiers", model$model, TreeNum = integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="HIBAG")
	if (rv$err != 0) stop(hlaErrMsg())

	# for each tree
	res <- vector("list", rv$TreeNum)
	for (i in 1:length(res))
	{
		# call, get the number of haplotypes
		rv <- .C("hlaAB_Idv_GetNumHaplo", model$model, as.integer(i),
			NumHaplo = integer(1), NumSNP = integer(1),
			err=integer(1), NAOK=TRUE, PACKAGE="HIBAG")
		if (rv$err != 0) stop(hlaErrMsg())

		# number of trios or samples
		if ("n.trio" %in% names(model))
			n.trio <- model$n.trio
		else
			n.trio <- model$n.samp

		# call, get freq. and haplotypes
		rv <- .C("hlaAB_Tree_GetHaplos", model$model, as.integer(i),
			freq=double(rv$NumHaplo), hla=integer(rv$NumHaplo), haplo=character(rv$NumHaplo),
			snpidx = integer(rv$NumSNP), samp.num = integer(n.trio), acc = double(1),
			err=integer(1), NAOK=TRUE, PACKAGE="HIBAG")
		if (rv$err != 0) stop(hlaErrMsg())

		res[[i]] <- list(
			samp.num = rv$samp.num,
			haplos = data.frame(freq=rv$freq, hla=model$hla.allele[rv$hla], haplo=rv$haplo,
				stringsAsFactors=FALSE),
			snpidx = rv$snpidx,
			outofbag.acc = rv$acc)
	}

	rv <- list(n.samp = model$n.samp, n.snp = model$n.snp)
	if ("n.trio" %in% names(model))
		rv$n.trio <- model$n.trio
	rv <- c(rv, list(
		sample.id = model$sample.id, snp.id = model$snp.id,
		snp.position = model$snp.position, snp.allele = model$snp.allele,
		snp.allele.freq = model$snp.allele.freq,
		hla.locus = model$hla.locus, hla.allele = model$hla.allele, hla.freq = model$hla.freq,
		classifiers = res))
	class(rv) <- "hlaAttrBagObj"
	return(rv)
}


#######################################################################
# To combine two model objects of attribute bagging
#

hlaCombineModelObj <- function(obj1, obj2)
{
	# check
	stopifnot(class(obj1) == "hlaAttrBagObj")
	stopifnot(class(obj2) == "hlaAttrBagObj")
	stopifnot(all(obj1$snp.id == obj2$snp.id))
	stopifnot(all(obj1$hla.allele == obj2$hla.allele))
	stopifnot(obj1$hla.locus == obj2$hla.locus)

	samp.id <- unique(c(obj1$sample.id, obj2$sample.id))

	rv <- list(n.samp = length(samp.id), n.snp = obj1$n.snp,
		sample.id = samp.id, snp.id = obj1$snp.id,
		snp.position = obj1$snp.position, snp.allele = obj1$snp.allele,
		snp.allele.freq = (obj1$snp.allele.freq + obj2$snp.allele.freq)*0.5,
		hla.locus = obj1$hla.locus,
		hla.allele = obj1$hla.allele, hla.freq = (obj1$hla.freq + obj2$hla.freq)*0.5,
		classifiers = c(obj1$classifiers, obj2$classifiers))
	if (!is.null(obj1$var.prior.prob) && !is.null(obj2$var.prior.prob))
		rv$var.prior.prob <- (obj1$var.prior.prob + obj2$var.prior.prob)*0.5
	class(rv) <- "hlaAttrBagObj"
	return(rv)
}


#######################################################################
# To get the top n individual classifiers
#

hlaSubModelObj <- function(obj, n)
{
	# check
	stopifnot(class(obj) == "hlaAttrBagObj")
	obj$classifiers <- obj$classifiers[1:n]
	return(obj)
}


#######################################################################
# To get a "hlaAttrBagClass" class
#

hlaModelfromObj <- function(obj)
{
	# check
	stopifnot(class(obj) == "hlaAttrBagObj")

	# create an attribute bagging object
	rv <- .C("hlaAB_Model_New",
		as.integer(obj$n.samp), as.integer(obj$n.snp), length(obj$hla.allele),
		model = integer(1), err=integer(1), NAOK=TRUE, PACKAGE="HIBAG")
	if (rv$err != 0) stop(hlaErrMsg())
	ABmodel <- rv$model

	# add individual classifiers
	for (tree in obj$classifiers)
	{
		hla <- match(tree$haplos$hla, obj$hla.allele) - 1
		if (any(is.na(hla)))
			stop("Invalid HLA alleles in the individual classifier.")
		if (is.null(tree$samp.num))
			snum <- rep(as.integer(1), obj$n.samp)
		else
			snum <- tree$samp.num
		rv <- .C("hlaAB_NewClassifierHaplo", ABmodel, length(tree$snpidx),
			as.integer(tree$snpidx-1), as.integer(snum), dim(tree$haplos)[1],
			as.double(tree$haplos$freq), as.integer(hla), as.character(tree$haplos$haplo),
			as.double(tree$outofbag.acc), err=integer(1), NAOK=TRUE, PACKAGE="HIBAG")
		if (rv$err != 0) stop(hlaErrMsg())
	}

	# output
	rv <- list(n.samp = obj$n.samp, n.snp = obj$n.snp)
	if ("n.trio" %in% names(obj))
		rv$n.trio <- obj$n.trio
	rv <- c(rv, list(
		sample.id = obj$sample.id, snp.id = obj$snp.id,
		snp.position = obj$snp.position, snp.allele = obj$snp.allele,
		snp.allele.freq = obj$snp.allele.freq,
		hla.locus = obj$hla.locus, hla.allele = obj$hla.allele, hla.freq = obj$hla.freq,
		model = ABmodel))
	class(rv) <- "hlaAttrBagClass"
	return(rv)
}


#######################################################################
# summarize the "hlaAttrBagObj" object
#

summary.hlaAttrBagObj <- function(object, show=TRUE, ...)
{
	# check
	stopifnot(class(object) == "hlaAttrBagObj")
	obj <- object

	if (show)
	{
		cat("HLA locus: ", obj$hla.locus, "\n", sep="")
		cat("Training dataset:", length(obj$sample.id), "samples X",
			length(obj$snp.id), "SNPs\n")
		cat("\t# of HLA alleles: ", length(obj$hla.allele), "\n", sep="")
	}

	# summarize ...
	snpset <- NULL
	outofbag.acc <- rep(NaN, length(obj$classifiers))
	numsnp <- rep(NA, length(obj$classifiers))
	numhaplo <- rep(NA, length(obj$classifiers))
	snp.hist <- rep(0, length(obj$snp.id))
	for (i in 1:length(obj$classifiers))
	{
		outofbag.acc[i] <- obj$classifiers[[i]]$outofbag.acc
		numsnp[i] <- length(obj$classifiers[[i]]$snpidx)
		numhaplo[i] <- length(obj$classifiers[[i]]$haplos$hla)
		snp.hist[obj$classifiers[[i]]$snpidx] <- snp.hist[obj$classifiers[[i]]$snpidx] + 1
		snpset <- unique(c(snpset, obj$classifiers[[i]]$snpidx))
	}
	snpset <- snpset[order(snpset)]
	outofbag.acc <- outofbag.acc * 100

	info <- data.frame(
		Mean = c(mean(numsnp), mean(numhaplo), mean(outofbag.acc)),
		SD = c(sd(numsnp), sd(numhaplo), sd(outofbag.acc)),
		Min = c(min(numsnp), min(numhaplo), min(outofbag.acc)),
		Max = c(max(numsnp), max(numhaplo), max(outofbag.acc))
	)
	rownames(info) <- c("num.snp", "num.haplo", "accuracy")

	if (show)
	{
		cat("\t# of individual classifiers: ", length(obj$classifiers), "\n", sep="")
		cat("\tTotal # of SNPs used: ", length(snpset), "\n", sep="")
		cat(sprintf("\tAverage # of SNPs in an individual classifier: %0.2f, sd: %0.2f, min: %d, max: %d\n",
			mean(numsnp), sd(numsnp), min(numsnp), max(numsnp)))
		cat(sprintf("\tAverage # of haplotypes in an individual classifier: %0.2f, sd: %0.2f, min: %d, max: %d\n",
			mean(numhaplo), sd(numhaplo), min(numhaplo), max(numhaplo)))
		cat(sprintf("\tAverage accuracy in an individual classifier: %0.2f%%, sd: %0.2f%%, min: %0.2f%%, max: %0.2f%%\n",
			mean(outofbag.acc), sd(outofbag.acc), min(outofbag.acc), max(outofbag.acc)))
	}

	rv <- list(num.classifier = length(obj$classifiers), num.snp = length(snpset),
		snp.id = obj$snp.id, snp.position = obj$snp.position,
		snp.hist = snp.hist, info = info)
	return(invisible(rv))
}


##########################################################################
# To calculate linkage disequilibrium between HLA locus and a SNP marker
#

hlaLD_HLA <- function(HLA, snp)
{
	# check
	stopifnot(class(HLA) == "hlaAlleleClass")
	stopifnot(is.vector(snp))
	stopifnot(is.numeric(snp))
	stopifnot(dim(HLA$value)[1] == length(snp))

	# HLA alleles
	flag <- is.finite(snp)
	alleles <- unique(c(HLA$value$allele1[flag], HLA$value$allele2[flag]))
	alleles <- alleles[order(alleles)]
	rv <- rep(NaN, length(alleles))
	for (i in 1:length(alleles))
	{
		x <- (HLA$value$allele1==alleles[i]) + (HLA$value$allele2==alleles[i])
		rv[i] <- cor(x, snp, use="complete.obs")^2
	}

	mean(rv, ra.rm=TRUE)
}


##########################################################################
# To calculate linkage disequilibrium between HLA locus and SNP markers
#

hlaGenoLD <- function(hla, geno)
{
	# check
	stopifnot(class(hla) == "hlaAlleleClass")
	if (class(geno) == "hlaSNPGenoClass")
	{
		stopifnot(dim(hla$value)[1] == length(geno$sample.id))
		geno <- geno$genotype
	} else if (is.matrix(geno))
	{
		stopifnot(dim(hla$value)[1] == dim(geno)[2])
	} else {
		stop("geno should be `hlaSNPGenoClass' or a matrix.")
	}

	# HLA alleles indicators
	alleles <- unique(c(hla$value$allele1, hla$value$allele2))
	alleles <- alleles[order(alleles)]
	allele.mat <- matrix(as.integer(0), nrow=length(hla$value$allele1), ncol=length(alleles))
	for (i in 1:length(alleles))
	{
		allele.mat[, i] <- (hla$value$allele1==alleles[i]) + (hla$value$allele2==alleles[i])
	}

	apply(geno, 1,
		function(x, allele.mat) {
			suppressWarnings(mean(cor(x, allele.mat, use="pairwise.complete.obs")^2, na.rm=TRUE))
		},
		allele.mat=allele.mat)
}


##########################################################################
# to get a model object of attribute bagging from a list of files
#

hlaAnonymize <- function(modelobj)
{
	# check
	stopifnot(class(modelobj) == "hlaAttrBagObj")
	
	modelobj$sample.id <- NULL
	for (i in 1:length(modelobj$classifiers))
	{
		modelobj$classifiers[[i]]$samp.num <- NULL
	}
	
	# output
	return(modelobj)
}


##########################################################################
# to get a model object of attribute bagging from a list of files
#

hlaModelFiles <- function(fn.list, action.missingfile=c("ignore", "stop"))
{
	# check
	stopifnot(is.character(fn.list))
	stopifnot(action.missingfile[1] %in% c("ignore", "stop"))

	# for-loop
	rv <- NULL
	for (fn in fn.list)
	{
		if (file.exists(fn))
		{
			tmp <- get(load(fn))
			if (is.null(rv))
			{
				rv <- tmp
			} else {
				rv <- hlaCombineModelObj(rv, tmp)
			}
		} else {
			if (action.missingfile[1] == "stop")
				stop(sprintf("Does not find '%s'.", fn))
		}
	}
	rv
}





##########################################################################
##########################################################################
#
# Visualization
#

##########################################################################
# To visualize an attribute bagging model
#

plot.hlaAttrBagClass <- function(x, ...)
{
	obj <- hlaModelToObj(x)
	plot(obj, ...)
}

print.hlaAttrBagClass <- function(x, ...)
{
	obj <- hlaModelToObj(x)
	print(obj)
	# output
	return(invisible(NULL))
}


##########################################################################
# To visualize an attribute bagging model
#

plot.hlaAttrBagObj <- function(x, xlab=NULL, ylab=NULL,
	locus.color="red", locus.lty=2, locus.cex=1.25, ...)
{
	# the starting and ending positions of HLA locus
	info <- hlaLociInfo()
	pos.start <- info$pos.HLA.start[[x$hla.locus]]/1000
	pos.end <- info$pos.HLA.end[[x$hla.locus]]/1000

	# summary of the attribute bagging model
	desp <- summary(x, show=FALSE)

	# x - label, y - label
	if (is.null(xlab)) xlab <- "snp position (kb)"
	if (is.null(ylab)) ylab <- "frequency of use"

	# draw
	plot(x$snp.position/1000, desp$snp.hist, xlab=xlab, ylab=ylab, ...)
	abline(v=pos.start, col=locus.color, lty=locus.lty)
	abline(v=pos.end, col=locus.color, lty=locus.lty)
	text((pos.start + pos.end)/2, max(desp$snp.hist), paste("HLA", x$hla.locus, sep="-"),
		col=locus.color, cex=locus.cex)
}

print.hlaAttrBagObj <- function(x, ...)
{
	summary(x)
	# output
	return(invisible(NULL))
}







#######################################################################
# To get the error message
#

hlaErrMsg <- function()
{
	rv <- .C("hlaErrMsg", msg=character(1), NAOK=TRUE, PACKAGE="HIBAG")
	rv$msg
}



#######################################################################
# Internal R library functions
#######################################################################

.onAttach <- function(lib, pkg)
{
	# load the dynamic-link library
	library.dynam("HIBAG", pkg, lib)
	# initialize HIBAG
	.C("hlaInit", PACKAGE="HIBAG")
	# set random number generator
	RNGkind("L'Ecuyer-CMRG")
	TRUE
}

.Last.lib <- function(libpath)
{
	# finalize HIBAG
	rv <- .C("hlaDone", PACKAGE="HIBAG")
	TRUE
}
