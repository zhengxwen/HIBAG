#######################################################################
#
# Package Name: HIBAG
# Description:
#   HIBAG -- HLA Genotype Imputation with Attribute Bagging
#
# HIBAG R package, HLA Genotype Imputation with Attribute Bagging
# Copyright (C) 2011-2015   Xiuwen Zheng (zhengx@u.washington.edu)
# All rights reserved.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

#######################################################################
#
# The internal functions
#
#######################################################################

#######################################################################
# get the name of HLA gene
#

.hla_gene_name_string <- function(geno.name)
{
    stopifnot(is.character(geno.name))
    i <- grep(paste("\\b",c(
            # HLA Classic I genes
            "A", "B", "C", "E", "F", "G",
            # HLA Classic II genes
            "DMA", "DMB", "DOA", "DOB",
            "DRA", "DRB\\d", "DQA\\d", "DQB\\d", "DPA\\d", "DPB\\d",
            # HLA Classic III genes
            "LTA", "TNF", "LTB", "HSPA1L", "HSPA1A", "HSPA1B",
            "C2", "BF", "C4A", "C4B"),
        "\\b", sep="", collapse="|"),
        geno.name)
    geno.name[i] <- paste("HLA -", geno.name[i])
    geno.name
}


#######################################################################
# get the name of HLA gene
#

.hla_assembly <- function(assembly =
    c("auto", "auto-silent", "hg18", "hg19", "hg20", "unknown"))
{
    assembly <- match.arg(assembly)
    if (assembly %in% c("auto", "auto-silent"))
    {
        if (assembly == "auto")
            message("using the default genome assembly (assembly=\"hg19\")")
        assembly <- "hg19"
    }
    assembly
}


#######################################################################
# call function in parallel
#

.DynamicClusterCall <- function(cl, fun, combine.fun, msg.fn, n,
    stop.cluster, ...)
{
    postNode <- function(con, type, value = NULL, tag = NULL)
    {
        eval(.SendData)
    }
    sendCall <- function(con, fun, args, return = TRUE, tag = NULL)
    {
        postNode(con, "EXEC",
            list(fun = fun, args = args, return = return, tag = tag))
        NULL
    }
    recvOneResult <- function(cl)
    {
        v <- eval(.RecvOneData)
        list(value = v$value$value, node = v$node, tag = v$value$tag)
    }


    #################################################################
    # check
    stopifnot(is.null(cl) | inherits(cl, "cluster"))
    stopifnot(is.function(fun))
    stopifnot(is.function(combine.fun))
    stopifnot(is.function(msg.fn))
    stopifnot(is.numeric(n))
    stopifnot(is.logical(stop.cluster))

    .SendData <- parse(text=
        "parallel:::sendData(con, list(type=type,data=value,tag=tag))")
    .RecvOneData <- parse(text="parallel:::recvOneData(cl)")

    val <- NULL
    if (!is.null(cl))
    {
        p <- length(cl)
        if (n > 0L && p)
        {
            ## **** this closure is sending to all nodes
            argfun <- function(i) c(i, list(...))

            submit <- function(node, job)
                sendCall(cl[[node]], fun, argfun(job), tag = job)

            for (i in 1:min(n, p)) submit(i, i)
            for (i in 1:n)
            {
                d <- recvOneResult(cl)
                j <- i + min(n, p)

                stopflag <- FALSE
                if (j <= n)
                {
                    submit(d$node, j)
                } else {
                    if (stop.cluster)
                    {
                        parallel::stopCluster(cl[d$node])
                        cl <- cl[-d$node]
                        stopflag <- TRUE
                    }
                }

                dv <- d$value
                if (inherits(dv, "try-error"))
                {
                    if (stop.cluster)
                        parallel::stopCluster(cl)
                    stop("One node produced an error: ", as.character(dv))
                }

                msg.fn(d$node, dv)
                if (stopflag)
                    message(sprintf("Stop \"job %d\".", d$node))

                val <- combine.fun(val, dv)
            }
        }
    } else {
        for (i in 1:n)
        {
            dv <- fun(i, ...)
            msg.fn(i, dv)
            val <- combine.fun(val, dv)
        }
    }

    val
}



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
#     assembly -- the human genome reference, such like "hg19"
#
#


#######################################################################
# To create a "hlaSNPGenoClass" object (SNP genotype object)
#

hlaMakeSNPGeno <- function(genotype, sample.id, snp.id, snp.position,
    A.allele, B.allele, assembly="auto")
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

    assembly <- .hla_assembly(assembly)

    rv <- list(genotype = genotype, sample.id = sample.id, snp.id = snp.id,
        snp.position = snp.position,
        snp.allele = paste(A.allele, B.allele, sep="/"),
        assembly = assembly)
    class(rv) <- "hlaSNPGenoClass"

    # valid snp.id
    flag <- is.na(rv$snp.id)
    if (any(flag))
    {
        warning("There is/are ", sum(flag),
            " SNP(s) with missing SNP id, and they have been removed.")
        rv <- hlaGenoSubset(rv, snp.sel=!flag)
    }
    # valid snp.position
    flag <- is.na(rv$snp.position)
    if (any(flag))
    {
        warning("There is/are ", sum(flag),
            " SNP(s) with missing SNP positions, and they have been removed.")
        rv <- hlaGenoSubset(rv, snp.sel=!flag)
    }

    rv
}


#######################################################################
# To select a subset of SNP genotypes
#

hlaGenoSubset <- function(genoobj, samp.sel=NULL, snp.sel=NULL)
{
    # check
    stopifnot(inherits(genoobj, "hlaSNPGenoClass"))
    stopifnot(is.null(samp.sel) | is.logical(samp.sel) | is.integer(samp.sel))
    if (is.logical(samp.sel))
        stopifnot(length(samp.sel) == length(genoobj$sample.id))
    stopifnot(is.null(snp.sel) | is.logical(snp.sel) | is.integer(snp.sel))
    if (is.logical(snp.sel))
        stopifnot(length(snp.sel) == length(genoobj$snp.id))
    if (is.integer(samp.sel))
    {
        stopifnot(!any(is.na(samp.sel)))
        stopifnot(length(unique(samp.sel)) == length(samp.sel))
    }
    if (is.integer(snp.sel))
    {
        stopifnot(!any(is.na(snp.sel)))
        stopifnot(length(unique(snp.sel)) == length(snp.sel))
    }

    # subset
    if (is.null(samp.sel))
        samp.sel <- rep(TRUE, length(genoobj$sample.id))
    if (is.null(snp.sel))
        snp.sel <- rep(TRUE, length(genoobj$snp.id))
    rv <- list(genotype = genoobj$genotype[snp.sel, samp.sel],
        sample.id = genoobj$sample.id[samp.sel],
        snp.id = genoobj$snp.id[snp.sel],
        snp.position = genoobj$snp.position[snp.sel],
        snp.allele = genoobj$snp.allele[snp.sel],
        assembly = genoobj$assembly
    )
    class(rv) <- "hlaSNPGenoClass"
    rv
}


#######################################################################
# To get the overlapping SNPs between target and template with
#   corrected strand.
#

hlaGenoSwitchStrand <- function(target, template,
    match.type=c("RefSNP+Position", "RefSNP", "Position"),
    same.strand=FALSE, verbose=TRUE)
{
    # check
    stopifnot(inherits(target, "hlaSNPGenoClass"))
    stopifnot(inherits(template, "hlaSNPGenoClass") |
        inherits(template, "hlaAttrBagClass") |
        inherits(template, "hlaAttrBagObj"))
    stopifnot(is.logical(same.strand))
    stopifnot(is.logical(verbose))
    match.type <- match.arg(match.type)

    # initialize
    s1 <- hlaSNPID(template, match.type)
    s2 <- hlaSNPID(target, match.type)

    flag <- TRUE
    if (length(s1) == length(s2))
    {
        if (all(s1 == s2, na.rm=TRUE))
        {
            s <- s1
            I1 <- I2 <- seq_len(length(s1))
            flag <- FALSE
        }
    }
    if (flag)
    {
        s <- intersect(s1, s2)
        if (length(s) <= 0) stop("There is no common SNP.")
        I1 <- match(s, s1); I2 <- match(s, s2)
    }

    # compute allele frequencies
    if (inherits(template, "hlaSNPGenoClass"))
    {
        template.afreq <- rowMeans(template$genotype, na.rm=TRUE) * 0.5
    } else {
        template.afreq <- template$snp.allele.freq
    }
    if (inherits(target, "hlaSNPGenoClass"))
    {
        target.afreq <- rowMeans(target$genotype, na.rm=TRUE) * 0.5
    } else {
        target.afreq <- rowMeans(target$haplotype, na.rm=TRUE)
    }

    # call
    gz <- .C(HIBAG_AlleleStrand,
        template$snp.allele, template.afreq, I1,
        target$snp.allele, target.afreq, I2,
        same.strand, length(s), out=logical(length(s)),
        out.n.ambiguity=integer(1), out.n.mismatching=integer(1),
        err=integer(1), NAOK=TRUE)
    if (gz$err != 0) stop(hlaErrMsg())

    if (verbose)
    {
        # switched allele pairs
        x <- sum(gz$out)
        if (x > 0)
        {
            if (x > 1)
            {
                a <- "are"; s <- "s"
            } else {
                a <- "is"; s <- ""
            }
            cat(sprintf(
    "There %s %d variant%s in total with switched allelic strand order%s.\n",
                a, x, s, s))
        } else {
            cat("No allelic strand orders are switched.\n")
        }

        # the number of ambiguity
        if (gz$out.n.ambiguity > 0)
        {
            if (gz$out.n.ambiguity > 1)
            {
                a <- "are"; s <- "s"
            } else {
                a <- "is"; s <- ""
            }
            cat("Due to stand ambiguity (such like C/G),",
                sprintf("the allelic strand order%s of %d variant%s %s",
                    s, gz$out.n.ambiguity, s, a),
                "determined by comparing allele frequencies.\n")
        }

        # the number of mismatching
        if (gz$out.n.mismatching > 0)
        {
            if (gz$out.n.mismatching > 1)
            {
                a <- "are"; s <- "s"
            } else {
                a <- "is"; s <- ""
            }
            cat("Due to mismatching alleles,",
                sprintf("the allelic strand order%s of %d variant%s %s",
                    s, gz$out.n.mismatching, s, a),
                "determined by comparing allele frequencies.\n")
        }
    }

    # result
    geno <- target$genotype[I2, ]
    if (is.vector(geno))
        geno <- matrix(geno, ncol=1)
    for (i in which(gz$out)) geno[i, ] <- 2 - geno[i, ]
    rv <- list(genotype = geno)
    rv$sample.id <- target$sample.id
    rv$snp.id <- target$snp.id[I2]
    rv$snp.position <- target$snp.position[I2]
    rv$snp.allele <- template$snp.allele[I1]
    rv$assembly <- template$assembly
    class(rv) <- "hlaSNPGenoClass"

    rv
}


#######################################################################
# To get the information of SNP ID and position
#

hlaSNPID <- function(obj, type=c("RefSNP+Position", "RefSNP", "Position"))
{
    stopifnot( inherits(obj, "hlaSNPGenoClass") |
        inherits(obj, "hlaAttrBagClass") | inherits(obj, "hlaAttrBagObj") )
    type <- match.arg(type)
    if (type == "RefSNP+Position")
        paste(obj$snp.id, obj$snp.position, sep="-")
    else if (type == "RefSNP")
        obj$snp.id
    else if (type == "Position")
        obj$snp.position
    else
        invisible()
}


#######################################################################
# To combine two SNP genotype dataset
#

hlaGenoCombine <- function(geno1, geno2,
    match.type=c("RefSNP+Position", "RefSNP", "Position"),
    allele.check=TRUE, same.strand=FALSE, verbose=TRUE)
{
    # check
    stopifnot(inherits(geno1, "hlaSNPGenoClass"))
    stopifnot(inherits(geno2, "hlaSNPGenoClass"))
    stopifnot(is.logical(allele.check))
    stopifnot(is.logical(same.strand))
    stopifnot(is.logical(verbose))
    match.type <- match.arg(match.type)

    if (allele.check)
    {
        tmp2 <- hlaGenoSwitchStrand(geno2, geno1, match.type,
            same.strand, verbose)
        tmp1 <- hlaGenoSubset(geno1, snp.sel=
            match(hlaSNPID(tmp2, match.type), hlaSNPID(geno1, match.type)))
    } else {
        s1 <- hlaSNPID(geno1, match.type)
        s2 <- hlaSNPID(geno2, match.type)
        set <- unique(intersect(s1, s2))
        tmp1 <- hlaGenoSubset(geno1, snp.sel=match(set, s1))
        tmp2 <- hlaGenoSubset(geno2, snp.sel=match(set, s2))
    }

    rv <- list(genotype = cbind(tmp1$genotype, tmp2$genotype),
        sample.id = c(tmp1$sample.id, tmp2$sample.id),
        snp.id = tmp1$snp.id, snp.position = tmp1$snp.position,
        snp.allele = tmp1$snp.allele,
        assembly = tmp1$assembly)
    colnames(rv$genotype) <- NULL
    rownames(rv$genotype) <- NULL

    class(rv) <- "hlaSNPGenoClass"
    rv
}


#######################################################################
# Convert to PLINK PED format
#

hlaGeno2PED <- function(geno, out.fn)
{
    # check
    stopifnot(inherits(geno, "hlaSNPGenoClass"))
    stopifnot(is.character(out.fn))

    # MAP file
    rv <- data.frame(chr=rep(6, length(geno$snp.id)), rs=geno$snp.id,
        morgan=rep(0, length(geno$snp.id)), bp=geno$snp.position,
        stringsAsFactors=FALSE)
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
        Paternal=rep(0, n), Maternal=rep(0, n),
        Sex=rep(0, n), Pheno=rep(-9, n), m)
    write.table(rv, file=paste(out.fn, ".ped", sep=""),
        row.names=FALSE, col.names=FALSE, quote=FALSE)

    # return
    invisible()
}


#######################################################################
# Convert from PLINK BED format
#

hlaBED2Geno <- function(bed.fn, fam.fn, bim.fn, rm.invalid.allele=FALSE,
    import.chr="xMHC", assembly="auto", verbose=TRUE)
{
    # check
    stopifnot(is.character(bed.fn) & (length(bed.fn)==1))
    stopifnot(is.character(fam.fn) & (length(fam.fn)==1))
    stopifnot(is.character(bim.fn) & (length(bim.fn)==1))
    stopifnot(is.character(import.chr))
    stopifnot(is.logical(rm.invalid.allele) & (length(rm.invalid.allele)==1))
    stopifnot(is.logical(verbose) & (length(verbose)==1))

    assembly <- .hla_assembly(assembly)

    # detect bed.fn
    bed <- .C(HIBAG_BEDFlag, bed.fn, snporder=integer(1), err=integer(1),
        NAOK=TRUE)
    if (bed$err != 0) stop(hlaErrMsg())
    if (verbose)
    {
        cat("Open \"", bed.fn, sep="")
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
        cat("Open \"", fam.fn, "\".\n", sep="")

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
    if (length(snp.id) != length(unique(snp.id)))
        stop("The SNP IDs in the PLINK binary file should be unique!")

    # snp allele
    snp.allele <- paste(bimD$allele1, bimD$allele2, sep="/")
    if (verbose)
        cat("Open \"", bim.fn, "\".\n", sep="")

    # SNP selection
    if (length(import.chr) == 1)
    {
        if (import.chr == "xMHC")
        {
            if (assembly %in% c("hg18", "NCBI36"))
            {
                snp.flag <- (chr==6) &
                    (25759242<=snp.pos) & (snp.pos<=33534827)
            } else if (assembly %in% c("hg19", "NCBI37"))
            {
                snp.flag <- (chr==6) &
                    (25651242<=snp.pos) & (snp.pos<=33544122)
            } else {
                stop("Invalid genome assembly.")
            }
            n.snp <- as.integer(sum(snp.flag))
            if (verbose)
            {
                cat(sprintf(
                    "Import %d SNPs within the xMHC region on chromosome 6.\n",
                    n.snp))
            }
            import.chr <- NULL
        } else if (import.chr == "")
        {
            n.snp <- length(snp.id)
            snp.flag <- rep(TRUE, n.snp)
            if (verbose)
                cat(sprintf("Import %d SNPs.\n", n.snp))
            import.chr <- NULL
        }
    }
    if (!is.null(import.chr))
    {
        snp.flag <- (chr %in% import.chr) & (snp.pos>0)
        n.snp <- as.integer(sum(snp.flag))
        if (verbose)
        {
            cat(sprintf("Import %d SNPs from chromosome %s.\n", n.snp,
                paste(import.chr, collapse=",")))
        }
    }
    if (n.snp <= 0) stop("There is no SNP imported.")

    # call the C function
    rv <- .C(HIBAG_ConvBED, bed.fn, length(sample.id), length(snp.id), n.snp,
        (bed$snporder==0), snp.flag, verbose,
        geno = matrix(as.integer(0), nrow=n.snp, ncol=length(sample.id)),
        err=integer(1), NAOK=TRUE)
    if (rv$err != 0) stop(hlaErrMsg())

    # result
    v <- list(genotype = rv$geno, sample.id = sample.id,
        snp.id = snp.id[snp.flag], snp.position = snp.pos[snp.flag],
        snp.allele = snp.allele[snp.flag], assembly = assembly)
    class(v) <- "hlaSNPGenoClass"

    # remove invalid snps
    if (rm.invalid.allele)
    {
        # check duplicated SNP ID
        flag <- duplicated(v$snp.id)
        if (any(flag))
        {
            if (verbose)
            {
                cat(sprintf("%d SNPs with duplicated ID have been removed.\n",
                    sum(flag)))
            }
            v <- hlaGenoSubset(v, snp.sel=!flag)
        }

        # check invalid alleles
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
        {
            cat(sprintf(
                "%d SNPs with invalid alleles have been removed.\n",
                sum(!flag)))
        }

        # get a subset
        v <- hlaGenoSubset(v, snp.sel=flag)
    }

    v
}


#######################################################################
# Convert from SNP GDS format (SNPRelate)
#

hlaGDS2Geno <- function(gds.fn, rm.invalid.allele=FALSE,
    import.chr="xMHC", assembly="auto", verbose=TRUE)
{
    # library
    if (!require(SNPRelate))
        stop("The SNPRelate package should be installed.")

    # check
    stopifnot(is.character(gds.fn) & is.vector(gds.fn))
    stopifnot(length(gds.fn) == 1)

    stopifnot(is.logical(rm.invalid.allele) & is.vector(rm.invalid.allele))
    stopifnot(length(rm.invalid.allele) == 1)

    stopifnot(is.character(import.chr))

    stopifnot(is.logical(verbose) & is.vector(verbose))
    stopifnot(length(verbose) == 1)

    assembly <- .hla_assembly(assembly)


    ####  open the GDS SNP file  ####

    chr <- NULL
    snp.pos <- NULL
    snp.id <- NULL

    gfile <- SNPRelate::snpgdsOpen(gds.fn)
    on.exit({ SNPRelate::snpgdsClose(gfile) })

    # snp.id
    snp.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gfile, "snp.id"))
    v <- gdsfmt::index.gdsn(gfile, "snp.rs.id", silent=TRUE)
    if (!is.null(v))
    {
        snp.rsid <- gdsfmt::read.gdsn(v)
    } else
        snp.rsid <- snp.id

    # chromosome
    chr <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gfile, "snp.chromosome"))

    # position
    snp.pos <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gfile, "snp.position"))
    snp.pos[!is.finite(snp.pos)] <- 0

    # SNP selection
    if (length(import.chr) == 1)
    {
        if (import.chr == "xMHC")
        {
            if (assembly %in% c("hg18", "NCBI36"))
            {
                snp.flag <- (chr==6) &
                    (25759242<=snp.pos) & (snp.pos<=33534827)
            } else if (assembly %in% c("hg19", "NCBI37"))
            {
                snp.flag <- (chr==6) &
                    (25651242<=snp.pos) & (snp.pos<=33544122)
            } else {
                stop("Invalid genome assembly.")
            }
            n.snp <- as.integer(sum(snp.flag))
            if (verbose)
            {
                cat(sprintf(
                    "Import %d SNPs within the xMHC region on chromosome 6.\n",
                    n.snp))
            }
            import.chr <- NULL
        } else if (import.chr == "")
        {
            n.snp <- length(snp.id)
            snp.flag <- rep(TRUE, n.snp)
            if (verbose)
                cat(sprintf("Import %d SNPs.\n", n.snp))
            import.chr <- NULL
        }
    }
    if (!is.null(import.chr))
    {
        snp.flag <- (chr %in% import.chr) & (snp.pos>0)
        n.snp <- as.integer(sum(snp.flag))
        if (verbose)
        {
            cat(sprintf("Import %d SNPs from chromosome %s.\n", n.snp,
                paste(import.chr, collapse=",")))
        }
    }
    if (n.snp <= 0) stop("There is no SNP imported.")

    # result
    v <- list(genotype = SNPRelate::snpgdsGetGeno(gfile,
            snp.id=snp.id[snp.flag], snpfirstdim=TRUE, verbose=FALSE),
        sample.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gfile, "sample.id")),
        snp.id = snp.rsid[snp.flag],
        snp.position = snp.pos[snp.flag],
        snp.allele = gdsfmt::read.gdsn(gdsfmt::index.gdsn(
            gfile, "snp.allele"))[snp.flag],
        assembly = assembly)
    class(v) <- "hlaSNPGenoClass"

    # remove invalid snps
    if (rm.invalid.allele)
    {
        # check duplicated SNP ID
        flag <- duplicated(v$snp.id)
        if (any(flag))
        {
            if (verbose)
            {
                cat(sprintf("%d SNPs with duplicated ID have been removed.\n",
                    sum(flag)))
            }
            v <- hlaGenoSubset(v, snp.sel=!flag)
        }

        # check invalid alleles
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
        {
            cat(sprintf("%d SNPs with invalid alleles have been removed.\n",
                sum(!flag)))
        }
        v <- hlaGenoSubset(v, snp.sel=flag)
    }

    v
}



#######################################################################
# Summarize a "hlaSNPGenoClass" object
#

summary.hlaSNPGenoClass <- function(object, show=TRUE, ...)
{
    # check
    stopifnot(inherits(object, "hlaSNPGenoClass"))
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
        cat("SNP genotypes: \n")
        cat(sprintf("\t%d samples X %d SNPs\n",
            length(geno$sample.id), length(geno$snp.id)))
        cat(sprintf("\tSNPs range from %dbp to %dbp",
            min(geno$snp.position, na.rm=TRUE),
            max(geno$snp.position, na.rm=TRUE)))
        if (!is.null(geno$assembly))
            cat(" on ", geno$assembly, "\n", sep="")
        else
            cat("\n")

        # missing rate for SNP
        cat(sprintf("Missing rate per SNP:\n\t%s\n", fn(rv$mr.snp)))
        # missing rate for sample
        cat(sprintf("Missing rate per sample:\n\t%s\n", fn(rv$mr.samp)))

        # minor allele frequency
        cat(sprintf("Minor allele frequency:\n\t%s\n", fn(rv$maf)))

        # allele information
        cat("Allelic information:")
        print(rv$allele)
    }

    # return
    invisible(rv)
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
    stopifnot(inherits(obj, "hlaSNPGenoClass"))
    rowMeans(obj$genotype, na.rm=TRUE) * 0.5
}


#######################################################################
# To the minor allele frequencies from genotypes or haplotypes
#

hlaGenoMFreq <- function(obj)
{
    # check
    stopifnot(inherits(obj, "hlaSNPGenoClass"))
    aF <- rowMeans(obj$genotype, na.rm=TRUE) * 0.5
    pmin(aF, 1 - aF)
}


#######################################################################
# To the missing rates from genotypes or haplotypes per SNP
#

hlaGenoMRate <- function(obj)
{
    # check
    stopifnot(inherits(obj, "hlaSNPGenoClass"))
    rowMeans(is.na(obj$genotype))
}


#######################################################################
# To the missing rates from genotypes or haplotypes per sample
#

hlaGenoMRate_Samp <- function(obj)
{
    # check
    stopifnot(inherits(obj, "hlaSNPGenoClass"))
    colMeans(is.na(obj$genotype))
}




#######################################################################
#
# the function list for HLA types
#
#######################################################################


#######################################################################
# To get the starting and ending positions in basepair for HLA loci
#

hlaLociInfo <- function(assembly =
    c("auto", "auto-silent", "hg18", "hg19", "hg20", "unknown"))
{
    # check
    assembly <- .hla_assembly(assembly)

    # file
    fn <- system.file("doc", sprintf("GeneInfo_%s.txt", assembly),
        package="HIBAG")
    if (file.exists(fn))
    {
        z <- read.table(fn, header=TRUE, stringsAsFactors=FALSE)
        rownames(z) <- z$name

        # output
        GenomicRanges::makeGRangesFromDataFrame(z)
    } else {
        if (assembly != "unknown")
            stop("Unknown human genome reference in 'assembly'!")
        invisible()
    }
}



#######################################################################
# Limit the resolution of HLA alleles
#

hlaAlleleDigit <- function(obj, max.resolution="4-digit", rm.suffix=FALSE)
{
    # check
    stopifnot(inherits(obj, "hlaAlleleClass") | is.character(obj))
    stopifnot(is.logical(rm.suffix))
    if (is.character(obj))
        stopifnot(is.vector(obj))
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
                            if (length(idx) == length(s))
                            {
                                if (rm.suffix & (nchar(s[length(s)])>0))
                                {
                                    z <- unlist(strsplit(s[length(s)], ""))
                                    for (i in 1:length(z))
                                    {
                                        if (!(z[i] %in% as.character(0:9)))
                                        {
                                            if (i > 1)
                                                z <- z[1:(i-1)]
                                            break
                                        }
                                    }
                                    s[length(s)] <- paste(z, collapse="")
                                }
                            }
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
                    stringsAsFactors=FALSE),
                assembly = obj$assembly
            )
            if ("prob" %in% names(obj$value))
                rv$value$prob <- obj$value$prob
            class(rv) <- "hlaAlleleClass"
            obj <- rv
        }
    }

    obj
}


#######################################################################
# Get unique HLA alleles
#

hlaUniqueAllele <- function(hla)
{
    # check
    stopifnot(is.character(hla) | inherits(hla, "hlaAlleleClass"))

    if (is.character(hla))
    {
        hla <- hla[!is.na(hla)]
        hla <- unique(hla)
        rv <- .C(HIBAG_SortAlleleStr, length(hla), hla,
            out = character(length(hla)),
            err = integer(1), NAOK=TRUE)
        if (rv$err != 0) stop(hlaErrMsg())
        rv$out
    } else {
        hlaUniqueAllele(as.character(c(hla$value$allele1, hla$value$allele2)))
    }
}


#######################################################################
# To make a class of HLA alleles
#

hlaAllele <- function(sample.id, H1, H2, max.resolution="", locus="any",
    assembly="auto", locus.pos.start=NA, locus.pos.end=NA, prob=NULL,
    na.rm=TRUE)
{
    # check
    stopifnot(is.vector(sample.id))
    stopifnot(is.vector(H1) & is.character(H1))
    stopifnot(is.vector(H2) & is.character(H2))
    stopifnot(length(sample.id) == length(H1))
    stopifnot(length(sample.id) == length(H2))
    stopifnot(max.resolution %in% c("2-digit", "4-digit", "6-digit",
        "8-digit", "allele", "protein", "2", "4", "6", "8", "full", ""))

    assembly <- .hla_assembly(assembly)
    HLAinfo <- hlaLociInfo(assembly)
    if (!is.null(prob))
        stopifnot(length(sample.id) == length(prob))

    # build
    H1[H1 == ""] <- NA
    H1 <- hlaAlleleDigit(H1, max.resolution)
    H2[H2 == ""] <- NA
    H2 <- hlaAlleleDigit(H2, max.resolution)

    if (locus %in% names(HLAinfo))
    {
        if (!is.finite(locus.pos.start))
            locus.pos.start <- BiocGenerics::start(HLAinfo[locus,])
        if (!is.finite(locus.pos.end))
            locus.pos.end <- BiocGenerics::end(HLAinfo[locus,])
    } else {
        locus.pos.start <- as.integer(NA)
        locus.pos.end <- as.integer(NA)
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
            allele1 = H1[flag], allele2 = H2[flag],
            stringsAsFactors=FALSE),
        assembly = assembly
    )
    if (!is.null(prob))
        rv$value$prob <- prob[flag]
    class(rv) <- "hlaAlleleClass"

    rv
}


#######################################################################
# To make a class of HLA alleles
#

hlaAlleleSubset <- function(hla, samp.sel=NULL)
{
    # check
    stopifnot(inherits(hla, "hlaAlleleClass"))
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
        value = hla$value[samp.sel, ],
        assembly = hla$assembly
    )

    if (!is.null(hla$postprob))
    {
        rv$postprob <- hla$postprob[, samp.sel]
        if (is.vector(rv$postprob))
        {
            rv$postprob <- matrix(rv$postprob, ncol=1)
        }
    }

    class(rv) <- "hlaAlleleClass"
    rv
}


#######################################################################
# To combine two classes of HLA alleles
#

hlaCombineAllele <- function(H1, H2)
{
    # check
    stopifnot(inherits(H1, "hlaAlleleClass"))
    stopifnot(inherits(H2, "hlaAlleleClass"))
    stopifnot(length(intersect(H1$sample.id, H2$sample.id)) == 0)
    stopifnot(H1$locus == H2$locus)
    stopifnot(H1$pos.start == H2$pos.start)
    stopifnot(H1$pos.end == H2$pos.end)

    id <- c("sample.id", "allele1", "allele2")

    # result
    rv <- list(locus = H1$locus,
        pos.start = H1$pos.start, pos.end = H1$pos.end,
        value = rbind(H1$value[, id], H2$value[, id]),
        assembly = H1$assembly
    )
    rownames(rv$value) <- NULL
    if (!is.null(H1$value$prob) & !is.null(H2$value$prob))
    {
        rv$value$prob <- c(H1$value$prob, H2$value$prob)
    }

    if (!is.null(H1$postprob) & !is.null(H2$postprob))
    {
        rv$postprob <- cbind(H1$postprob, H2$postprob)
    }

    class(rv) <- "hlaAlleleClass"
    rv
}


#######################################################################
# To compare HLA alleles
#

hlaCompareAllele <- function(TrueHLA, PredHLA, allele.limit=NULL,
    call.threshold=NaN, max.resolution="", output.individual=FALSE,
    verbose=TRUE)
{
    # check
    stopifnot(inherits(TrueHLA, "hlaAlleleClass"))
    stopifnot(inherits(PredHLA, "hlaAlleleClass"))
    stopifnot(is.null(allele.limit) |
        (is.vector(allele.limit) & is.character(allele.limit)) |
        inherits(allele.limit, "hlaAttrBagClass") |
        inherits(allele.limit, "hlaAttrBagObj"))
    stopifnot(max.resolution %in% c("2-digit", "4-digit", "6-digit",
        "8-digit", "allele", "protein", "2", "4", "6", "8", "full", ""))
    stopifnot(is.logical(output.individual))
    stopifnot(is.logical(verbose))

    # get the common samples
    samp <- intersect(TrueHLA$value$sample.id, PredHLA$value$sample.id)
    if ((length(samp) != length(TrueHLA$value$sample.id)) |
        (length(samp) != length(PredHLA$value$sample.id)))
    {
        if (verbose)
        {
            message("Calling 'hlaCompareAllele': there are ", length(samp),
                " individuals in common.\n")
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
        if (inherits(allele.limit, "hlaAttrBagClass") |
            inherits(allele.limit, "hlaAttrBagObj"))
        {
            allele <- hlaUniqueAllele(allele.limit$hla.allele)
            TrainFreq <- allele.limit$hla.freq
            TrainNum <- allele.limit$n.samp
        } else {
            allele <- hlaUniqueAllele(as.character(allele.limit))
            TrainFreq <- NULL
            TrainNum <- NaN
        }
    } else {
        allele <- hlaUniqueAllele(c(ts1, ts2))
        TrainFreq <- NULL
        TrainNum <- NaN
    }

    # max resolution
    if (!(max.resolution %in% c("full", "")))
    {
        ts1 <- hlaAlleleDigit(ts1, max.resolution)
        ts2 <- hlaAlleleDigit(ts2, max.resolution)
        ps1 <- hlaAlleleDigit(ps1, max.resolution)
        ps2 <- hlaAlleleDigit(ps2, max.resolution)
        tmp <- hlaAlleleDigit(allele, max.resolution)
        allele <- hlaUniqueAllele(tmp)
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
    cnt.ind <- 0L; cnt.haplo <- 0L; cnt.call <- 0L
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
        { if (x %in% LT) x else "..." }

    acc.array <- rep(NaN, n)
    ind.truehla <- character(n)
    ind.predhla <- character(n)

    if (n > 0L)
    {
        for (i in 1:n)
        {
            # increase
            TrueNumAll[[ ts1[i] ]] <- TrueNumAll[[ ts1[i] ]] + 1L
            TrueNumAll[[ ts2[i] ]] <- TrueNumAll[[ ts2[i] ]] + 1L

            # probability cut-off
            if (is.null(prob))
                flag <- TRUE
            else
                flag <- (prob[i] >= call.threshold)
            if (flag)
            {
                # update TrueNum and PredNum
                TrueNum[[ ts1[i] ]] <- TrueNum[[ ts1[i] ]] + 1L
                TrueNum[[ ts2[i] ]] <- TrueNum[[ ts2[i] ]] + 1L
                PredNum[[ fn(ps1[i], allele) ]] <-
                    PredNum[[ fn(ps1[i], allele) ]] + 1L
                PredNum[[ fn(ps2[i], allele) ]] <-
                    PredNum[[ fn(ps2[i], allele) ]] + 1L

                # correct count of individuals
                if ( ((ts1[i]==ps1[i]) & (ts2[i]==ps2[i])) |
                    ((ts2[i]==ps1[i]) & (ts1[i]==ps2[i])) )
                {
                    cnt.ind <- cnt.ind + 1
                }

                # correct count of haplotypes
                s <- c(ts1[i], ts2[i]); p <- c(ps1[i], ps2[i])
                ind.truehla[i] <- paste(s[order(s)], collapse="/")
                ind.predhla[i] <- paste(p[order(p)], collapse="/")
                
                hnum <- 0L
                if ((s[1]==p[1]) | (s[1]==p[2]))
                {
                    if (s[1]==p[1]) { p[1] <- "" } else { p[2] <- "" }
                    confusion[s[1], s[1]] <- confusion[s[1], s[1]] + 1L
                    cnt.haplo <- cnt.haplo + 1L
                    hnum <- hnum + 1L
                }
                if ((s[2]==p[1]) | (s[2]==p[2]))
                {
                    confusion[s[2], s[2]] <- confusion[s[2], s[2]] + 1L
                    cnt.haplo <- cnt.haplo + 1L
                    hnum <- hnum + 1L
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
                            confusion[fn(p[2], allele), s[2]] <-
                                confusion[fn(p[2], allele), s[2]] + 1L
                        } else {
                            confusion[fn(p[1], allele), s[2]] <-
                                confusion[fn(p[1], allele), s[2]] + 1L
                        }
                    } else {
                        if (s[2]==p[1])
                        {
                            confusion[fn(p[2], allele), s[1]] <-
                                confusion[fn(p[2], allele), s[1]] + 1L
                        } else {
                            confusion[fn(p[1], allele), s[1]] <-
                                confusion[fn(p[1], allele), s[1]] + 1L
                        }
                    }
                } else if (hnum == 0)
                {
                    WrongTab <- cbind(WrongTab,
                        c(s, fn(p[1], allele), fn(p[2], allele)))
                }

                # the number of calling
                cnt.call <- cnt.call + 1L
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
        overall$call.threshold <- 0L
    }

    # confusion matrix
    if (is.null(WrongTab))
    {
        nw <- as.integer(0)
    } else {
        nw <- ncol(WrongTab)
    }
    rv <- .C(HIBAG_Confusion, as.integer(m), confusion,
        nw, match(WrongTab, names(PredNum)) - as.integer(1),
        out = matrix(0.0, nrow=m+1, ncol=m, dimnames=
            list(Predict=names(PredNum), True=names(TrueNum))),
        tmp = double((m+1)*m),
        err = integer(1), NAOK = TRUE)
    if (rv$err != 0) stop(hlaErrMsg())
    confusion <- round(rv$out, 2)

    # detail -- sensitivity and specificity
    detail <- data.frame(allele = allele, stringsAsFactors=FALSE)
    if (!is.null(TrainFreq))
    {
        detail$train.num <- 2 * TrainFreq * TrainNum
        detail$train.freq <- TrainFreq
    }
    detail$valid.num <- TrueNumAll
    detail$valid.freq <- TrueNumAll / sum(TrueNumAll)
    detail$call.rate <- TrueNum / TrueNumAll

    sens <- diag(confusion) / TrueNum
    spec <- 1 - (PredNum[1:m] - diag(confusion)) / (2*cnt.call - TrueNum)
    detail$accuracy <- (sens*TrueNum + spec*(2*cnt.call - TrueNum)) /
        (2*cnt.call)
    detail$sensitivity <- sens
    detail$specificity <- spec
    detail$ppv <- diag(confusion) / rowSums(confusion)[1:m]
    detail$npv <- 1 - (TrueNum - diag(confusion)) /
        (2*n - rowSums(confusion)[1:m])

    detail$call.rate[!is.finite(detail$call.rate)] <- 0
    detail[detail$call.rate<=0,
        c("sensitivity", "specificity", "ppv", "npv", "accuracy")] <- NaN

    # get miscall
    rv <- confusion; diag(rv) <- 0
    m.max <- apply(rv, 2, max); m.idx <- apply(rv, 2, which.max)
    s <- names(PredNum)[m.idx]; s[m.max<=0] <- NA
    p <- m.max / apply(rv, 2, sum)
    detail <- cbind(detail, miscall=s, miscall.prop=p, stringsAsFactors=FALSE)
    rownames(detail) <- NULL

    # output
    rv <- list(overall=overall, confusion=confusion, detail=detail)
    if (output.individual)
    {
        rv$individual <- data.frame(sample.id=samp.id,
            true.hla=ind.truehla, pred.hla=ind.predhla,
            accuracy=acc.array, stringsAsFactors=FALSE)
    }
    rv
}


#######################################################################
# Return a sample list satisfying the filter conditions
#

hlaSampleAllele <- function(TrueHLA, allele.limit=NULL, max.resolution="")
{
    # check
    stopifnot(inherits(TrueHLA, "hlaAlleleClass"))
    stopifnot(is.null(allele.limit) | is.vector(allele.limit) |
        inherits(allele.limit, "hlaAttrBagClass") |
        inherits(allele.limit, "hlaAttrBagObj"))
    stopifnot(max.resolution %in% c("2-digit", "4-digit", "6-digit",
        "8-digit", "allele", "protein", "2", "4", "6", "8", "full", ""))

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
        if (inherits(allele.limit, "hlaAttrBagClass") |
            inherits(allele.limit, "hlaAttrBagObj"))
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
    TrueHLA$value$sample.id[flag]
}


#######################################################################
# Divide the list of HLA types to the training and validation sets
#

hlaSplitAllele <- function(HLA, train.prop=0.5)
{
    # check
    stopifnot(inherits(HLA, "hlaAlleleClass"))

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
        samp.id <- H$value$sample.id[ H$value$allele1==allele |
            H$value$allele2==allele ]
        samp.id <- as.character(samp.id)

        n.train <- ceiling(length(samp.id) * train.prop)
        train.sampid <- sample(samp.id, n.train)
        train.set <- c(train.set, train.sampid)

        H <- hlaAlleleSubset(H, samp.sel=
            match(setdiff(H$value$sample.id, samp.id), H$value$sample.id))
    }

    train.set <- train.set[order(train.set)]
    list(
        training = hlaAlleleSubset(HLA,
            samp.sel = match(train.set, HLA$value$sample.id)),
        validation = hlaAlleleSubset(HLA,
            samp.sel = match(setdiff(HLA$value$sample.id, train.set),
                HLA$value$sample.id)
            )
    )
}


#######################################################################
# To select SNPs in the flanking region of a specified HLA locus
#

hlaFlankingSNP <- function(snp.id, position, hla.id, flank.bp=500*1000,
    assembly="auto")
{
    # check
    stopifnot(length(snp.id) == length(position))
    stopifnot(is.character(hla.id))

    # init
    assembly <- .hla_assembly(assembly)
    HLAInfo <- hlaLociInfo(assembly)
    ID <- names(HLAInfo)

    stopifnot(length(hla.id) == 1)
    if (!(hla.id %in% ID))
        stop(paste("`hla.id' should be one of", paste(ID, collapse=",")))

    pos.start <- BiocGenerics::start(HLAInfo[hla.id,]) - flank.bp
    pos.end <- BiocGenerics::end(HLAInfo[hla.id,]) + flank.bp

    if (is.finite(pos.start) & is.finite(pos.end))
    {
        flag <- (pos.start <= position) & (position <= pos.end)
        snp.id[flag]
    } else
        stop("The position information is not available!")
}


#######################################################################
# Summary a "hlaAlleleClass" object
#

summary.hlaAlleleClass <- function(object, show=TRUE, ...)
{
    # check
    stopifnot(inherits(object, "hlaAlleleClass"))
    hla <- object

    HUA <- hlaUniqueAllele(c(hla$value$allele1, hla$value$allele2))
    HLA <- factor(match(c(hla$value$allele1, hla$value$allele2), HUA))
    levels(HLA) <- HUA  
    count <- table(HLA)
    freq <- prop.table(count)
    rv <- cbind(count=count, freq=freq)

    # get the number of unique genotypes
    m <- data.frame(a1 = hla$value$allele1,
        a2 = hla$value$allele2, stringsAsFactors=FALSE)
    lst <- apply(m, 1, FUN=function(x) {
        if (!is.na(x[1]) & !is.na(x[2]))
        {
            if (x[1] <= x[2])
                paste(x[1], x[2], sep="/")
            else
                paste(x[2], x[1], sep="/")
        } else
            NA
    })
    unique.n.geno <- nlevels(factor(lst))

    if (show)
    {
        cat("Gene: ", .hla_gene_name_string(hla$locus), "\n", sep="")
        cat(sprintf("Range: [%dbp, %dbp]", hla$pos.start, hla$pos.end))
        if (!is.null(hla$assembly))
            cat(" on ", hla$assembly, "\n", sep="")
        else
            cat("\n")
        cat(sprintf("# of samples: %d\n", dim(hla$value)[1]))
        cat(sprintf("# of unique HLA alleles: %d\n", length(count)))
        cat(sprintf("# of unique HLA genotypes: %d\n", unique.n.geno))

        p <- hla$value$prob
        if (!is.null(p))
        {
            z <- table(cut(p, breaks=c(0, 0.25, 0.5, 0.75, 1),
                right=FALSE, include.lowest=TRUE))
            z[] <- sprintf("%d (%0.1f%%)", z, prop.table(z)*100)
            names(attr(z, "dimnames")) <- "Posterior probability:"
            print(z)
        }
    }

    # return
    invisible(rv)
}


#######################################################################
# Check missing SNP predictors
#

hlaCheckSNPs <- function(model, object,
    match.type=c("RefSNP+Position", "RefSNP", "Position"), verbose=TRUE)
{
    # check
    stopifnot(inherits(model, "hlaAttrBagClass") |
        inherits(model, "hlaAttrBagObj"))
    stopifnot((is.vector(object) & is.character(object)) |
        inherits(object, "hlaSNPGenoClass"))
    match.type <- match.arg(match.type)
    stopifnot(is.logical(verbose))

    # initialize
    if (inherits(model, "hlaAttrBagClass"))
        model <- hlaModelToObj(model)

    # show information
    if (verbose)
    {
        cat("The HIBAG model:\n")
        cat(sprintf("\tThere are %d SNP predictors in total.\n",
            length(model$snp.id)))
        cat(sprintf("\tThere are %d individual classifiers.\n",
            length(model$classifiers)))
    }

    if (is.vector(object))
    {
        target.snp <- as.character(object)
        src.snp <- hlaSNPID(model, match.type)
    } else {
        target.snp <- hlaSNPID(object, match.type)
        src.snp <- hlaSNPID(model, match.type)
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
        cat("Summarize",
            "the missing fractions of SNP predictors per classifier:\n")
        print(summary(1 - rv$fraction))
    }

    # output
    invisible(rv)
}


##########################################################################
# to finalize the HIBAG model
#

hlaPublish <- function(mobj, platform=NULL, information=NULL, warning=NULL,
    rm.unused.snp=TRUE, anonymize=TRUE, verbose=TRUE)
{
    # check
    stopifnot(inherits(mobj, "hlaAttrBagObj") |
        inherits(mobj, "hlaAttrBagClass"))
    stopifnot(is.null(platform) | is.character(platform))
    stopifnot(is.null(information) | is.character(information))
    stopifnot(is.null(warning) | is.character(warning))
    stopifnot(is.logical(rm.unused.snp) & (length(rm.unused.snp)==1))
    stopifnot(is.logical(verbose))
    if (inherits(mobj, "hlaAttrBagClass"))
        mobj <- hlaModelToObj(mobj)

    # additional information
    if (is.null(platform))
        platform <- mobj$appendix$platform
    if (is.null(information))
        information <- mobj$appendix$information
    if (is.null(warning))
        warning <- mobj$appendix$warning
    mobj$appendix <- list(
        platform=platform, information=information, warning=warning)

    # remove unused SNPs
    if (rm.unused.snp)
    {
        # get frequency of use for SNPs
        snp.hist <- rep(0, length(mobj$snp.id))
        for (i in 1:length(mobj$classifiers))
        {
            idx <- mobj$classifiers[[i]]$snpidx
            snp.hist[idx] <- snp.hist[idx] + 1
        }

        flag <- (snp.hist > 0)
        if (sum(flag) < mobj$n.snp)
        {
            if (verbose)
            {
                cnt <- mobj$n.snp - sum(flag)
                cat(sprintf("Remove %d unused SNP%s.\n",
                    cnt, if (cnt>1) "s" else ""))
            }
        }
        mobj$n.snp <- sum(flag)
        mobj$snp.id <- mobj$snp.id[flag]
        mobj$snp.position <- mobj$snp.position[flag]
        mobj$snp.allele <- mobj$snp.allele[flag]
        mobj$snp.allele.freq <- mobj$snp.allele.freq[flag]

        idx.list <- rep(0, length(flag))
        idx.list[flag] <- 1:mobj$n.snp

        for (i in 1:length(mobj$classifiers))
        {
            mobj$classifiers[[i]]$snpidx <-
                idx.list[ mobj$classifiers[[i]]$snpidx ]
        }
    }

    # anonymize
    if (anonymize)
    {
        mobj$sample.id <- NULL
        for (i in 1:length(mobj$classifiers))
        {
            mobj$classifiers[[i]]$samp.num <- NULL
        }
    }

    # output
    mobj
}


##########################################################################
# to get a model object of attribute bagging from a list of files
#

hlaModelFiles <- function(fn.list, action.missingfile=c("ignore", "stop"),
    verbose=TRUE)
{
    # check
    stopifnot(is.character(fn.list))
    stopifnot(is.logical(verbose))
    action.missingfile <- match.arg(action.missingfile)

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
            s <- sprintf("There is no '%s'.", fn)
            if (action.missingfile == "stop")
            {
                stop(s)
            } else {
                if (verbose) message(s)
            }
        }
    }
    rv
}


##########################################################################
# Out-of-bag estimation of overall accuracy, per-allele sensitivity, etc
#

hlaOutOfBag <- function(model, hla, snp, call.threshold=NaN, verbose=TRUE)
{
    # check
    stopifnot(inherits(model, "hlaAttrBagObj") |
        inherits(model, "hlaAttrBagClass"))
    stopifnot(inherits(hla, "hlaAlleleClass"))
    stopifnot(inherits(snp, "hlaSNPGenoClass"))

    stopifnot(is.numeric(call.threshold) & is.vector(call.threshold))
    stopifnot(length(call.threshold) == 1)

    stopifnot(is.logical(verbose) & is.vector(verbose))
    stopifnot(length(verbose) == 1)


    ######################################################
    # initialize ...
    if (inherits(model, "hlaAttrBagClass"))
    {
        model <- hlaModelToObj(model)
        if (verbose) print(model)
    }

    # map samples
    if (is.null(model$sample.id))
        stop("There is no sample ID in the model.")
    samp.idx <- match(model$sample.id, snp$sample.id)
    if (any(is.na(samp.idx)))
        stop("Some of sample.id in the model do not exist in SNP genotypes.")
    hla.samp.idx <- match(model$sample.id, hla$value$sample.id)
    if (any(is.na(hla.samp.idx)))
        stop("Some of sample.id in the model do not exist in HLA types.")

    # map SNPs
    snp.idx <- match(model$snp.id, snp$snp.id)
    if (any(is.na(snp.idx)))
        stop("Some of snp.id in the model do not exist in SNP genotypes.")

    # genotypes and the number of classifiers
    geno <- snp$genotype[snp.idx, samp.idx]
    nclass <- length(model$classifiers)

    # the returned value
    ans <- NULL

    # the column names of details
    nm1 <- c("allele", "train.num", "train.freq")
    nm2 <- c("call.rate", "accuracy", "sensitivity", "specificity",
        "ppv", "npv")

    # for-loop
    for (i in 1:nclass)
    {
        mx <- model
        mx$classifiers <- mx$classifiers[i]
        s <- mx$classifiers[[1]]$samp.num
        if (is.null(s))
            stop("There is no bootstrap sample index.")

        tmp.model <- hlaModelFromObj(mx)
        tmp.geno <- geno[, s == 0]
        v <- predict(tmp.model, tmp.geno, verbose=FALSE)
        hlaClose(tmp.model)
        v$value$sample.id <- mx$sample.id[s == 0]
        pam <- hlaCompareAllele(hla, v, allele.limit=mx,
            call.threshold=call.threshold, verbose=FALSE)

        if (!is.null(ans))
        {
            # overall
            ans$overall <- ans$overall + pam$overall

            # confusion matrix
            ans$confusion <- ans$confusion + pam$confusion

            # details
            pam$detail <- pam$detail[, nm2]
            ans$n.detail <- ans$n.detail + !is.na(pam$detail)
            pam$detail[is.na(pam$detail)] <- 0
            ans$detail <- ans$detail + pam$detail
        } else {
            ans <- pam
            ans$detailhead <- ans$detail[, nm1]
            colnames(ans$detailhead) <- c("allele", "valid.num", "valid.freq")
            ans$detail <- ans$detail[, nm2]
            ans$n.detail <- !is.na(ans$detail)
            ans$detail[is.na(ans$detail)] <- 0
        }

        if (verbose)
        {
            cat(date(), sprintf(", passing the %d/%d classifiers.\n",
                i, nclass), sep="")
        }
    }

    # average
    ans$overall <- ans$overall / nclass
    ans$confusion <- ans$confusion / nclass
    ans$detail <- ans$detail / ans$n.detail

    # get miscall
    rv <- ans$confusion; diag(rv) <- 0
    m.max <- apply(rv, 2, max); m.idx <- apply(rv, 2, which.max)
    s <- rownames(ans$confusion)[m.idx]; s[m.max<=0] <- NA
    p <- m.max / apply(rv, 2, sum)

    # output
    ans$detail <- cbind(ans$detailhead, ans$detail,
        miscall=s, miscall.prop=p, stringsAsFactors=FALSE)
    ans$detailhead <- NULL
    ans$n.detail <- NULL
    ans
}


##########################################################################
# to create a report for evaluating accuracies
#

hlaReport <- function(object, export.fn="", type=c("txt", "tex", "html"),
    header=TRUE)
{
    # check
    stopifnot(is.list(object))
    stopifnot(is.data.frame(object$detail))
    stopifnot(is.character(export.fn))
    type <- match.arg(type)
    stopifnot(is.logical(header))

    # create an output file
    if (export.fn != "")
    {
        f <- file(export.fn, "wt")
        on.exit(close(f))
    } else {
        f <- ""
    }

    ###############################################################
    # text format

    d <- data.frame(Allele=object$detail$allele, stringsAsFactors=FALSE)

    if (!is.null(object$detail$train.num))
        d$NumTrain <- object$detail$train.num
    if (!is.null(object$detail$train.freq))
    {
        d$FreqTrain <- sprintf("%0.4f", object$detail$train.freq)
        d$FreqTrain[d$NumTrain <= 0] <- "0"

        L1 <- c("Allele", "Num.", "Freq.", "Num.", "Freq.",
            "CR", "ACC", "SEN", "SPE", "PPV", "NPV", "Miscall")
        L2 <- c("", "Train", "Train", "Valid.", "Valid.",
            "(%)", "(%)", "(%)", "(%)", "(%)", "(%)", "(%)")
        L2a <- c("", "Train", "Train", "Valid.", "Valid.",
            "(\\%)", "(\\%)", "(\\%)", "(\\%)", "(\\%)", "(\\%)", "(\\%)")
    } else {
        L1 <- c("Allele", "Num.", "Freq.",
            "CR", "ACC", "SEN", "SPE", "PPV", "NPV", "Miscall")
        L2 <- c("", "Valid.", "Valid.",
            "(%)", "(%)", "(%)", "(%)", "(%)", "(%)", "(%)")
        L2a <- c("", "Valid.", "Valid.",
            "(\\%)", "(\\%)", "(\\%)", "(\\%)", "(\\%)", "(\\%)", "(\\%)")
    }

    d$NumValid <- object$detail$valid.num
    d$FreqValid <- sprintf("%0.4f", object$detail$valid.freq)
    d$FreqValid[d$NumValid <= 0] <- "0"

    d$CR <- sprintf("%0.1f", object$detail$call.rate*100)
    d$CR[object$detail$call.rate < 0.0005] <- "--"

    d$ACC <- sprintf("%0.1f", object$detail$accuracy*100)
    d$ACC[!is.finite(object$detail$accuracy)] <- "--"

    d$SEN <- sprintf("%0.1f", object$detail$sensitivity*100)
    d$SEN[!is.finite(object$detail$sensitivity)] <- "--"

    d$SPE <- sprintf("%0.1f", object$detail$specificity*100)
    d$SPE[!is.finite(object$detail$specificity)] <- "--"

    d$PPV <- sprintf("%0.1f", object$detail$ppv*100)
    d$PPV[!is.finite(object$detail$ppv)] <- "--"

    d$NPV <- sprintf("%0.1f", object$detail$npv*100)
    d$NPV[!is.finite(object$detail$npv)] <- "--"

    d$Miscall <- sprintf("%s (%0.0f)",
        object$detail$miscall, object$detail$miscall.prop*100)
    d$Miscall[is.na(object$detail$miscall) |
        !is.finite(object$detail$miscall.prop)] <- "--"


    if (type == "txt")
    {
        cat(L1, file=f, sep="\t", append=TRUE)
        cat("\n", file=f, append=TRUE)
        cat(L2, file=f, sep="\t", append=TRUE)
        cat("\n", file=f, append=TRUE)
        cat("----\n", file=f, append=TRUE)
        cat(sprintf("Overall accuracy: %0.1f%%, Call rate: %0.1f%%\n",
            object$overall$acc.haplo*100, object$overall$call.rate*100),
            file=f, append=TRUE)
        write.table(d, file=f, append=TRUE, quote=FALSE, sep=" ",
            row.names=FALSE, col.names=FALSE)

    } else if (type == "tex")
    {
        if (header)
        {
            cat("\\title{Imputation Evaluation}", "",
                "\\documentclass[12pt]{article}", "",
                "\\usepackage{fullpage}",
                "\\usepackage{longtable}", "",
                "\\begin{document}", "",
                "\\maketitle", "",
                "\\setlength{\\LTcapwidth}{6.5in}", "",
                file=f, append=TRUE, sep="\n")
        }

        cat("% -------- BEGIN TABLE --------\n", file=f, append=TRUE)
        if (!is.null(object$detail$train.freq))
        {
            cat("\\begin{longtable}{rrrrr | rrrrrrl}\n", file=f, append=TRUE)
        } else {
            cat("\\begin{longtable}{rrr | rrrrrrl}\n", file=f, append=TRUE)
        }
        cat(c("\\caption{The sensitivity (SEN), specificity (SPE),",
            "positive predictive value (PPV), negative predictive value",
            "(NPV) and call rate (CR).}\n"), file=f, append=TRUE, sep=" ")
        cat("\\label{tab:accuracy} \\\\\n", file=f, append=TRUE)

        cat(L1, file=f, sep=" & ", append=TRUE)
        cat(" \\\\\n", file=f, append=TRUE)
        cat(L2a, file=f, sep=" & ", append=TRUE)
        cat(" \\\\\n", file=f, append=TRUE)
        cat("\\hline\\hline\n\\endfirsthead\n", file=f, append=TRUE)

        if (!is.null(object$detail$train.freq))
        {
            cat("\\multicolumn{12}{c}{{\\normalsize \\tablename\\",
                "\\thetable{} -- Continued from previous page}} \\\\\n",
                file=f, append=TRUE)
        } else {
            cat("\\multicolumn{10}{c}{{\\normalsize \\tablename\\",
                "\\thetable{} -- Continued from previous page}} \\\\\n",
                file=f, append=TRUE)
        }

        cat(L1, file=f, sep=" & ", append=TRUE)
        cat(" \\\\\n", file=f, append=TRUE)
        cat(L2a, file=f, sep=" & ", append=TRUE)
        cat(" \\\\\n", file=f, append=TRUE)
        cat("\\hline\\hline\n\\endhead\n\\hline\n", file=f, append=TRUE)

        if (!is.null(object$detail$train.freq))
        {
            cat("\\multicolumn{12}{r}{Continued on next page ...} \\\\\n",
                file=f, append=TRUE)
        } else {
            cat("\\multicolumn{10}{r}{Continued on next page ...} \\\\\n",
                file=f, append=TRUE)
        }

        cat("\\hline\n\\endfoot\n\\hline\\hline\n\\endlastfoot\n",
            file=f, append=TRUE)

        if (!is.null(object$detail$train.freq))
        {
            cat("\\multicolumn{12}{l}{\\it Overall accuracy:",
                sprintf("%0.1f\\%%, Call rate: %0.1f\\%%} \\\\\n",
                object$overall$acc.haplo*100, object$overall$call.rate*100),
                file=f, append=TRUE)
        } else {
            cat("\\multicolumn{10}{l}{\\it Overall accuracy:",
                sprintf("%0.1f\\%%, Call rate: %0.1f\\%%} \\\\\n",
                object$overall$acc.haplo*100, object$overall$call.rate*100),
                file=f, append=TRUE)
        }

        write.table(d, file=f, append=TRUE, quote=FALSE, sep=" & ",
            row.names=FALSE, col.names=FALSE, eol=" \\\\\n")

        cat("\\end{longtable}\n% -------- END TABLE --------\n",
            file=f, append=TRUE)

        if (header)
        {
            cat("\n\\end{document}\n", file=f, append=TRUE)
        }
    } else if (type == "html")
    {
        if (header)
        {
            cat("<!DOCTYPE html>",
                "<html>",
                "<head>",
                "  <title>Imputation Evaluation</title>",
                "</head>",
                "<body>",
                file=f, append=TRUE, sep="\n")
        }

        cat("<h1>Imputation Evaluation</h1>",
            "<p></p>",
            "<h3><b>Table 1:</b> The sensitivity (SEN), specificity (SPE),",
            "positive predictive value (PPV), negative predictive value (NPV)",
            "and call rate (CR).</h3>",
            file=f, append=TRUE, sep="\n")

        cat(
    "<table id=\"TB-Acc\" class=\"tabular\" border=\"1\"  CELLSPACING=\"1\">",
            "<tr>", 
            paste(paste("<th>", L1, " ", L2, "</th>", sep=""), collapse=" "),
            "</tr>",
            "<tr>",
            paste("<td colspan=\"", length(L1), "\">", sep=""),
            sprintf("<i> Overall accuracy: %0.1f%%, Call rate: %0.1f%% </i>",
                object$overall$acc.haplo*100, object$overall$call.rate*100),
            "</td>",
            "</tr>",
            file=f, append=TRUE, sep="\n")

        for (i in 1:nrow(d))
        {
            cat("<tr>", 
                paste(paste("<td>", d[i, ], "</td>", sep=""), collapse=" "),
                "</tr>",
                file=f, append=TRUE, sep="\n")
        }

        cat("</table>\n", file=f, append=TRUE)

        if (header)
        {
            cat("\n</body>\n</html>\n", file=f, append=TRUE)
        }
    }

    invisible()
}
