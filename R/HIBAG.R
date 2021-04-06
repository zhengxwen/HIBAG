#######################################################################
#
# Package Name: HIBAG
# Description:
#   HIBAG -- HLA Genotype Imputation with Attribute Bagging
#
# HIBAG R package, HLA Genotype Imputation with Attribute Bagging
# Copyright (C) 2011-2020   Xiuwen Zheng (zhengx@u.washington.edu)
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


##########################################################################
#
# Attribute Bagging method -- HIBAG algorithm
#

.printMatching <- function(p)
{
    w <- summary(p)
    class(w) <- "table"
    w <- c(w[1L],
        "0.1% Qu." = quantile(p, 0.001, na.rm=TRUE, names=FALSE),
        "1% Qu."   = quantile(p, 0.01, na.rm=TRUE, names=FALSE),
        w[2L], w[3L], w[5L], w[6L],
        w[4L], "SD" = sd(p, na.rm=TRUE))
    print(w)
}


##########################################################################
# Fit a HIBAG model for imputing HLA genotypes
#

hlaAttrBagging <- function(hla, snp, nclassifier=100L,
    mtry=c("sqrt", "all", "one"), prune=TRUE, na.rm=TRUE, mono.rm=TRUE, maf=NaN,
    nthread=1L, verbose=TRUE, verbose.detail=FALSE)
{
    # check
    stopifnot(inherits(hla, "hlaAlleleClass"))
    stopifnot(inherits(snp, "hlaSNPGenoClass"))
    stopifnot(is.numeric(nclassifier), length(nclassifier)==1L)
    stopifnot(is.character(mtry) | is.numeric(mtry), length(mtry)>0L)
    stopifnot(is.logical(prune), length(prune)==1L)
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    stopifnot(is.logical(mono.rm), length(mono.rm)==1L)
    stopifnot(is.numeric(maf), length(maf)==1L)
    stopifnot(is.numeric(nthread) | is.logical(nthread), length(nthread)==1L,
        !is.na(nthread))
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.detail), length(verbose.detail)==1L)
    if (verbose.detail) verbose <- TRUE

    if (is.na(nclassifier)) nclassifier <- 0L
    with.in.call <- nclassifier==0L
    with.matching <- (nclassifier > 0L)
    if (!with.matching)
    {
        nclassifier <- -nclassifier
        if (nclassifier == 0L) nclassifier <- 1L
    }

    # get the common samples
    samp.id <- intersect(hla$value$sample.id, snp$sample.id)

    # hla types
    samp.flag <- match(samp.id, hla$value$sample.id)
    hla.allele1 <- hla$value$allele1[samp.flag]
    hla.allele2 <- hla$value$allele2[samp.flag]
    if (na.rm)
    {
        if (any(is.na(c(hla.allele1, hla.allele2))))
        {
            warning("There are missing HLA alleles, ",
                "and the corresponding samples have been removed.")
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
    samp.flag <- match(samp.id, snp$sample.id)
    snp.geno <- snp$genotype[, samp.flag]
    if (!is.integer(snp.geno))
        storage.mode(snp.geno) <- "integer"
    if (with.in.call)
    {
        # save the variable, until remove it later
        .packageEnv$snp.geno <- snp.geno
    }

    tmp.snp.id <- snp$snp.id
    tmp.snp.position <- snp$snp.position
    tmp.snp.allele <- snp$snp.allele

    # remove mono-SNPs and MAF
    msg <- NULL
    if (mono.rm || is.finite(maf))
    {
        msg <- paste0("    MAF threshold: ", maf, "\n")
        mf <- rowMeans(snp.geno, na.rm=TRUE) * 0.5
        mf <- pmin(mf, 1-mf)
        mf[!is.finite(mf)] <- 0
        snpsel <- rep(TRUE, length(mf))
        if (mono.rm)
        {
            n0 <- sum(snpsel)
            snpsel <- snpsel & (mf > 0)
            n1 <- sum(snpsel)
            a <- n0 - n1
            if (a > 0L)
            {
                msg <- paste0(msg, "    excluding ", a, " monomorphic SNP",
                    .plural(a), "\n")
            }
        }
        if (is.finite(maf))
        {
            n0 <- sum(snpsel)
            snpsel <- snpsel & (mf >= maf)
            n1 <- sum(snpsel)
            a <- n0 - n1
            if (a > 0L)
            {
                msg <- paste0(msg, "    excluding ", a, " SNP", .plural(a),
                    " for MAF threshold\n")
            }
        }
        if (!all(snpsel))
        {
            tmp.snp.id <- tmp.snp.id[snpsel]
            tmp.snp.position <- tmp.snp.position[snpsel]
            tmp.snp.allele <- tmp.snp.allele[snpsel]
            snp.geno <- snp.geno[snpsel, , drop=FALSE]
        }
    }

    # check
    if (length(samp.id) <= 0L)
        stop("There is no common sample between 'hla' and 'snp'.")
    if (length(dim(snp.geno)[1L]) <= 0L)
        stop("There is no valid SNP markers.")


    ###################################################################
    # initialize ...

    n.snp <- dim(snp.geno)[1L]     # Num. of SNPs
    n.samp <- dim(snp.geno)[2L]    # Num. of samples
    HUA <- hlaUniqueAllele(c(hla.allele1, hla.allele2))
    H <- factor(match(c(hla.allele1, hla.allele2), HUA))
    levels(H) <- HUA
    n.hla <- nlevels(H)
    H1 <- as.integer(H[1L:n.samp]) - 1L
    H2 <- as.integer(H[(n.samp+1L):(2L*n.samp)]) - 1L

    # create an attribute bagging object (return an integer)
    ABmodel <- .Call(HIBAG_Training, n.snp, n.samp, snp.geno, n.hla, H1, H2)

    # number of variables randomly sampled as candidates at each split
    mtry <- mtry[1L]
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
            mtry <- 1L
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
    if (mtry <= 0) mtry <- 1L

    if (verbose)
    {
        cat(sprintf("Build a HIBAG model with %d individual classifier%s:\n",
            nclassifier, .plural(nclassifier)))
        cat(msg)
        cat("    # of SNPs randomly sampled as candidates for each selection: ",
            mtry, "\n", sep="")
        cat("    # of SNPs: ", n.snp, "\n", sep="")
        cat("    # of samples: ", n.samp, "\n", sep="")
        s <- ifelse(!grepl("^KIR", hla$locus), "HLA", "KIR")
        cat("    # of unique ", s, " alleles: ", n.hla, "\n", sep="")
        cat("CPU flags: ", .Call(HIBAG_Kernel_Version)[[2L]][1L], "\n", sep="")
    }

    # set the number of threads (initialized in HIBAG_NewClassifiers)
    if (isTRUE(nthread))
        nthread <- as.integer(defaultNumThreads())
    if (is.na(nthread) || (nthread<1L)) nthread <- 1L


    ###################################################################
    # training ...
    # add new individual classifers
    .Call(HIBAG_NewClassifiers, ABmodel, nclassifier, mtry, prune,
        nthread, verbose, verbose.detail, NULL)

    # output
    mod <- list(n.samp = n.samp, n.snp = n.snp, sample.id = samp.id,
        snp.id = tmp.snp.id, snp.position = tmp.snp.position,
        snp.allele = tmp.snp.allele,
        snp.allele.freq = 0.5*rowMeans(snp.geno, na.rm=TRUE),
        hla.locus = hla$locus, hla.allele = levels(H),
        hla.freq = prop.table(table(H)),
        assembly = as.character(snp$assembly)[1L],
        model = ABmodel,
        appendix = list())
    if (is.na(mod$assembly)) mod$assembly <- "unknown"
    if (with.in.call) mod$mtry <- mtry

    class(mod) <- "hlaAttrBagClass"


    ###################################################################
    # calculate matching proportion
    if (with.matching)
    {
        if (verbose)
            cat("Calculating matching proportion:\n")
        pd <- hlaPredict(mod, snp, cl=nthread, match.type="Pos+Allele", verbose=FALSE)
        mod$matching <- pd$value$matching
        if (verbose)
        {
            .printMatching(mod$matching)
            acc <- hlaCompareAllele(hla, pd, verbose=FALSE)$overall$acc.haplo
            cat(sprintf("Accuracy with training data: %.2f%%\n", acc*100))
            # out-of-bag accuracy
            mobj <- hlaModelToObj(mod)
            acc <- sapply(mobj$classifiers, function(x) x$outofbag.acc)
            cat(sprintf("Out-of-bag accuracy: %.2f%%\n", mean(acc)*100))
        }
    }

    # output
    mod
}


##########################################################################
# Fit a HIBAG model for imputing HLA genotypes in parallel
#

.show_model_obj <- function(mobj, autosave)
{
    z <- summary(mobj, show=FALSE)
    cat(ifelse(autosave, "==Saved==", " --"),
        paste0("#", length(mobj$classifier), ","), "avg oob acc:",
        sprintf("%0.2f%%, sd: %0.2f%%, min: %0.2f%%, max: %0.2f%%\n",
        z$info["accuracy", "Mean"], z$info["accuracy", "SD"],
        z$info["accuracy", "Min"], z$info["accuracy", "Max"]))
    invisible()
}

hlaParallelAttrBagging <- function(cl, hla, snp, auto.save="",
    nclassifier=100L, mtry=c("sqrt", "all", "one"), prune=TRUE, na.rm=TRUE,
    mono.rm=TRUE, maf=NaN, stop.cluster=FALSE, verbose=TRUE, verbose.detail=FALSE)
{
    # check
    stopifnot(is.null(cl) | is.logical(cl) | is.numeric(cl) |
        inherits(cl, "cluster"))
    stopifnot(inherits(hla, "hlaAlleleClass"))
    stopifnot(inherits(snp, "hlaSNPGenoClass"))
    stopifnot(is.character(auto.save), length(auto.save)==1L, !is.na(auto.save))
    if (auto.save!="" && is.na(.fn_obj_check(auto.save)))
        stop("'auto.save' should be a .rda/.RData or .rds file name.")
    stopifnot(is.numeric(nclassifier), length(nclassifier)==1L, nclassifier>0L)
    stopifnot(is.character(mtry) | is.numeric(mtry), length(mtry)>0L)
    stopifnot(is.logical(prune), length(prune)==1L)
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    stopifnot(is.logical(mono.rm), length(mono.rm)==1L)
    stopifnot(is.numeric(maf), length(maf)==1L)
    stopifnot(is.logical(stop.cluster), length(stop.cluster)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.detail), length(verbose.detail)==1L)

    if (verbose)
    {
        cat("Building a HIBAG model:\n")
        cat(sprintf("    %d individual classifier%s\n", nclassifier,
            .plural(nclassifier)))
        if (!is.null(cl))
        {
            cat(sprintf("    run in parallel with %d compute node%s\n",
                length(cl), .plural(length(cl))))
        }
        if (auto.save != "")
            cat("    autosave to ", sQuote(auto.save), "\n", sep="")
    }

    if (inherits(cl, "cluster"))
    {
        if (!requireNamespace("parallel", quietly=TRUE))
            stop("The `parallel' package should be installed.")
        if (verbose)
            cat("[-] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep="")
        # set random number for the cluster
        RNGkind("L'Ecuyer-CMRG")
        rand <- .Random.seed
        parallel::clusterSetRNGStream(cl)

        total <- 0L
        ans <- .DynamicClusterCall(cl,
            fun = function(job, hla, snp, mtry, prune, na.rm, mono.rm)
            {
                model <- hlaAttrBagging(hla=hla, snp=snp, nclassifier=-1L,
                    mtry=mtry, prune=prune, na.rm=na.rm, mono.rm=mono.rm, maf=maf,
                    verbose=FALSE, verbose.detail=FALSE)
                mobj <- hlaModelToObj(model)
                hlaClose(model)
                mobj
            },
            combine.fun = function(obj1, obj2)
            {
                if (is.null(obj1))
                    mobj <- obj2
                else if (is.null(obj2))
                    mobj <- obj1
                else
                    mobj <- hlaCombineModelObj(obj1, obj2)
                if (auto.save != "")
                    .fn_obj_save(auto.save, mobj)
                if (verbose & !is.null(mobj))
                    .show_model_obj(mobj, auto.save!="")
                mobj
            },
            msg.fn = function(job, obj)
            {
                if (verbose)
                {
                    z <- summary(obj, show=FALSE)
                    total <<- total + 1L
                    cat(sprintf(
                        "[%d] %s, worker%3d, # of SNPs: %g, # of haplo: %g, oob acc: %0.1f%%\n",
                        total, format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                        as.integer(job), z$info["num.snp", "Mean"],
                        z$info["num.haplo", "Mean"],
                        z$info["accuracy", "Mean"]))
                }
            },
            n=nclassifier, stop.cluster=stop.cluster,
            hla=hla, snp=snp, mtry=mtry, prune=prune,
            na.rm=na.rm, mono.rm=mono.rm
        )

        # the next random seed
        parallel::nextRNGStream(rand)
        parallel::nextRNGSubStream(rand)
        if (stop.cluster) cl <- NULL
        if (auto.save != "")
            ans <- .fn_obj_load(auto.save)
        mod <- hlaModelFromObj(ans)

    } else {
        # find the number of threads
        nthread <- 1L
        if (isTRUE(nthread))
            nthread <- as.integer(defaultNumThreads())
        if (is.numeric(cl)) nthread <- as.integer(cl)
        if (is.na(nthread) || (nthread<1L)) nthread <- 1L

        if (auto.save == "")
        {
            mod <- hlaAttrBagging(hla=hla, snp=snp, nclassifier=-nclassifier,
                mtry=mtry, prune=prune, na.rm=na.rm, mono.rm=mono.rm, maf=maf,
                nthread=nthread, verbose=verbose, verbose.detail=verbose.detail)
        } else {
            mod <- hlaAttrBagging(hla=hla, snp=snp, nclassifier=0L,
                mtry=mtry, prune=prune, na.rm=na.rm, mono.rm=mono.rm, maf=maf,
                nthread=nthread, verbose=verbose, verbose.detail=verbose.detail)
            mtry <- mod$mtry
            on.exit({ .packageEnv$snp.geno <- NULL }, add=TRUE)
            mobj <- hlaModelToObj(mod)
            .fn_obj_save(auto.save, mobj)
            if (verbose) .show_model_obj(mobj, TRUE)
            # add the remaining individual classifiers
            for (i in seq_len(nclassifier-1L))
            {
                .Call(HIBAG_NewClassifiers, mod$model, 1L, mtry, prune,
                    -nthread, verbose, verbose.detail, NULL)
                mobj <- hlaModelToObj(mod)
                .fn_obj_save(auto.save, mobj)
                if (verbose) .show_model_obj(mobj, TRUE)
            }
        }
    }

    # matching proportion
    if (verbose)
        cat("Calculating matching proportion:\n")
    if (is.null(cl)) cl <- FALSE
    pd <- hlaPredict(mod, snp, cl=cl, verbose=FALSE)
    mod$matching <- pd$value$matching
    mobj <- NULL
    if (auto.save != "")
    {
        mobj <- hlaModelToObj(mod)
        .fn_obj_save(auto.save, mobj)
    }
    if (verbose)
    {
        .printMatching(mod$matching)
        acc <- hlaCompareAllele(hla, pd, verbose=FALSE)$overall$acc.haplo
        cat(sprintf("Accuracy with training data: %.2f%%\n", acc*100))
        # out-of-bag accuracy
        if (is.null(mobj)) mobj <- hlaModelToObj(mod)
        acc <- sapply(mobj$classifiers, function(x) x$outofbag.acc)
        cat(sprintf("Out-of-bag accuracy: %.2f%%\n", mean(acc)*100))
    }

    # output
    if (auto.save != "") invisible() else mod
}


##########################################################################
# Close and dispose a HIBAG model
#

hlaClose <- function(model)
{
    stopifnot(inherits(model, "hlaAttrBagClass"))
    .Call(HIBAG_Close, model$model)
    invisible()
}


#######################################################################
# Predict HLA types using unphased SNP data
#

predict.hlaAttrBagClass <- function(object, snp, cl=FALSE,
    type=c("response", "dosage", "prob", "response+prob"), vote=c("prob", "majority"),
    allele.check=TRUE, match.type=c("Position", "Pos+Allele", "RefSNP+Position", "RefSNP"),
    same.strand=FALSE, verbose=TRUE, verbose.match=TRUE, ...)
{
    stopifnot(inherits(object, "hlaAttrBagClass"))
    hlaPredict(object, snp, cl, type, vote, allele.check, match.type,
        same.strand, verbose, verbose.match)
}

hlaPredict <- function(object, snp, cl=FALSE,
    type=c("response", "dosage", "prob", "response+prob"), vote=c("prob", "majority"),
    allele.check=TRUE, match.type=c("Position", "Pos+Allele", "RefSNP+Position", "RefSNP"),
    same.strand=FALSE, verbose=TRUE, verbose.match=TRUE)
{
    # check
    stopifnot(inherits(object, "hlaAttrBagClass"))
    stopifnot(is.logical(cl) | is.numeric(cl) | inherits(cl, "cluster"))
    stopifnot(is.logical(allele.check), length(allele.check)==1L)
    stopifnot(is.logical(same.strand), length(same.strand)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(is.logical(verbose.match), length(verbose.match)==1L)
    type <- match.arg(type)
    vote <- match.arg(vote)
    match.type <- match.arg(match.type)
    vote_method <- match(vote, c("prob", "majority"))

    if (inherits(cl, "cluster"))
    {
        if (!requireNamespace("parallel", quietly=TRUE))
            stop("The `parallel' package should be installed.")
        if (length(cl) <= 1L) cl <- FALSE
    }

    # if warning
    if (!is.null(object$appendix$warning))
        warning(object$appendix$warning, immediate.=TRUE)

    if (verbose)
    {
        # get the number of classifiers
        CNum <- .Call(HIBAG_GetNumClassifiers, object$model)
        s <- object$hla.allele
        if (length(s) > 3L) s <- c(s[1:3], "...")
        cat("HIBAG model for ", .hla_gene_name_string(object$hla.locus), ":\n",
            "    ", CNum, " individual classifier", .plural(CNum), "\n",
            "    ", length(object$snp.id), " SNPs\n",
            "    ", length(object$hla.allele), " unique HLA alleles: ",
            paste(s, collapse=", "), "\n", sep="")
        cat("Prediction:\n")
        if (vote_method == 1L)
            cat("    based on the averaged posterior probabilities\n")
        else
            cat("    by voting from all individual classifiers\n")
        if (inherits(cl, "cluster"))
        {
            cat("    run in parallel with ", length(cl), " compute node",
                .plural(length(cl)), "\n", sep="")
        }
    }

    if (!inherits(snp, "hlaSNPGenoClass"))
    {
        # it should be a vector or a matrix
        stopifnot(is.numeric(snp), is.vector(snp) | is.matrix(snp))
        if (is.vector(snp))
        {
            stopifnot(length(snp) == object$n.snp)
            snp <- matrix(snp, ncol=1L)
        } else {
            stopifnot(nrow(snp) == object$n.snp)
        }
        geno.sampid <- 1L:ncol(snp)
        assembly <- "auto-silent"

    } else {

        ##################################################
        # a 'hlaSNPGenoClass' object, check assembly first

        model.assembly <- as.character(object$assembly)[1L]
        if (is.na(model.assembly))
            model.assembly <- "unknown"
        geno.assembly <- as.character(snp$assembly)[1L]
        if (is.na(geno.assembly))
            geno.assembly <- "unknown"
        refstr <- sprintf("Model assembly: %s, SNP assembly: %s",
            model.assembly, geno.assembly)
        if (verbose)
            cat(refstr, "\n", sep="")

        if (model.assembly != geno.assembly)
        {
            if (any(c(model.assembly, geno.assembly) %in% "unknown"))
            {
                if (verbose)
                    message("The human genome references might not match!")
                if (geno.assembly == "unknown")
                    assembly <- model.assembly
                else
                    assembly <- geno.assembly
            } else {
                warning("The human genome references do not match! ",
                    refstr, ".", immediate.=TRUE)
                assembly <- model.assembly
            }
        } else {
            if (model.assembly != "unknown")
                assembly <- model.assembly
            else
                assembly <- "auto"
        }

        ##################################################

        # verbose for different matching methods
        if (verbose && verbose.match)
        {
            cat("Matching the SNPs between the model and the test data:\n")
            tab <- NULL
            for (tp in c("Position", "Pos+Allele", "RefSNP+Position", "RefSNP"))
            {
                s1 <- hlaSNPID(object, tp)
                s2 <- hlaSNPID(snp, tp)
                s <- unique(intersect(s1, s2))
                mcnt <- length(s1) - length(s)
                d <- data.frame(c1=paste0("   ", tp),
                    c2=sprintf("%d (%.1f%%)", mcnt, mcnt/length(s1)*100),
                    c3=ifelse(tp==match.type, "*being used", ""),
                    stringsAsFactors=FALSE)
                tab <- rbind(tab, d)
            }
            names(tab) <- c("match.type=\"--\"", "  missing SNPs #", "")
            tab[1L,3L] <- paste(tab[1L,3L], "[1]")
            tab[2L,3L] <- paste(tab[2L,3L], "[2]")
            print(tab, row.names=FALSE)
            cat("      [1]: useful if ambiguous strands on array-based platforms\n")
            cat("      [2]: suggested if the model and test data have been matched to the same reference genome\n")
            s <- object$appendix$platform
            if (is.null(s)) s <- "not applicable" else s <- paste(s, collapse=",")
            cat("    Model platform: ", s, "\n", sep="")
        } else if (verbose)
        {
            cat("Using match.type='", match.type, "' for SNP matching\n", sep="")
        }

        geno.sampid <- snp$sample.id
        obj.id <- hlaSNPID(object, match.type)
        geno.id <- hlaSNPID(snp, match.type)

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
                snp <- hlaGenoSwitchStrand(snp, object,
                    match.type, same.strand, verbose)$genotype
            } else {
                snp <- snp$genotype
            }
        } else {

            # snp selection
            snp.sel <- match(obj.id, geno.id)
            snp.sel[duplicated(snp.sel)] <- NA_integer_
            snp.allele <- snp$snp.allele
            snp.allele[is.na(snp.allele)] <- ""

            # temporary variable
            tmp <- list(genotype = snp$genotype[snp.sel, , drop=FALSE],
                sample.id = snp$sample.id,
                snp.id = object$snp.id, snp.position = object$snp.position,
                snp.allele = snp.allele[snp.sel],
                assembly = snp$assembly)
            flag <- is.na(tmp$snp.allele)
            tmp$snp.allele[flag] <- object$snp.allele[
                match(tmp$snp.id[flag], object$snp.id)]
            class(tmp) <- "hlaSNPGenoClass"

            # the total number of missing snp
            missing.cnt <- sum(flag)

            if (missing.cnt == length(obj.id))
            {
                stop("There is no overlapping of SNPs!")
            } else if (missing.cnt > 0.5*length(obj.id))
            {
                warning("More than 50% of SNPs are missing!", immediate.=TRUE)
            }

            # switch
            if (allele.check)
            {
                snp <- hlaGenoSwitchStrand(tmp, object, match.type,
                    same.strand, verbose)$genotype
            } else {
                snp <- tmp$genotype
            }
        }
    }

    # check
    if (dim(snp)[1L] != object$n.snp)
    {
        stop("The number of SNPs is not valid, and it maybe due to duplicated 'snp.id' ",
            "or incorrect dimension of genotype matrix.")
    }

    # initialize ...
    n.samp <- dim(snp)[2L]
    n.hla <- length(object$hla.allele)
    if (verbose)
    {
        cat(sprintf("# of samples: %d\n", n.samp))
        cat("CPU flags: ", .Call(HIBAG_Kernel_Version)[[2L]][1L], "\n", sep="")
    }

    # parallel units
    if (!inherits(cl, "cluster"))
    {
        # the number of threads
        nthread <- 1L
        if (isTRUE(cl)) nthread <- as.integer(defaultNumThreads())
        if (is.numeric(cl)) nthread <- cl[1L]
        if (is.na(nthread) || nthread<1L) nthread <- 1L

        # pointer to functions for an extensible component
        pm <- attr(cl, "proc_ptr")

        # predict HLA genotypes
        if (type %in% c("response", "dosage", "response+prob"))
        {
            # the best-guess prediction
            if (type == "response")
            {
                rv <- .Call(HIBAG_Predict_Resp, object$model, as.integer(snp),
                    n.samp, vote_method, nthread, verbose, pm)
                names(rv) <- c("H1", "H2", "prob", "matching")
            } else if (type == "dosage")
            {
                rv <- .Call(HIBAG_Predict_Dosage, object$model, as.integer(snp),
                    n.samp, vote_method, nthread, verbose, pm)
                names(rv) <- c("H1", "H2", "prob", "matching", "dosage")
            } else {
                rv <- .Call(HIBAG_Predict_Resp_Prob, object$model,
                    as.integer(snp), n.samp, vote_method, nthread, verbose, pm)
                names(rv) <- c("H1", "H2", "prob", "matching", "postprob")
            }

            # output object
            res <- hlaAllele(geno.sampid,
                H1 = object$hla.allele[rv$H1 + 1L],
                H2 = object$hla.allele[rv$H2 + 1L],
                locus = object$hla.locus, prob = rv$prob,
                na.rm = FALSE, assembly = assembly)
            res$value$matching <- rv$matching
            if (!is.null(rv$dosage))
            {
                res$dosage <- rv$dosage
                rownames(res$dosage) <- object$hla.allele
                colnames(res$dosage) <- geno.sampid
            }
            if (!is.null(rv$postprob))
            {
                res$postprob <- rv$postprob
                colnames(res$postprob) <- geno.sampid
                m <- outer(object$hla.allele, object$hla.allele,
                    function(x, y) paste(y, x, sep="/"))
                rownames(res$postprob) <- m[lower.tri(m, diag=TRUE)]
            }
            NA.cnt <- sum(is.na(res$value$allele1) | is.na(res$value$allele2))

        } else {
            # all probabilities
            rv <- .Call(HIBAG_Predict_Resp_Prob, object$model,
                as.integer(snp), n.samp, vote_method, nthread, verbose, pm)
            names(rv) <- c("H1", "H2", "prob", "matching", "postprob")

            # output object
            res <- rv$postprob
            colnames(res) <- geno.sampid
            m <- outer(object$hla.allele, object$hla.allele,
                function(x, y) paste(y, x, sep="/"))
            rownames(res) <- m[lower.tri(m, diag=TRUE)]
            NA.cnt <- sum(colSums(res) <= 0L, na.rm=TRUE)
        }
    } else {

        # run in parallel
        rv <- parallel::clusterApply(cl=cl,
            parallel::splitIndices(n.samp, length(cl)),
            fun = function(i, mobj, snp, type, vote)
            {
                if (length(i) > 0L)
                {
                    library(HIBAG)
                    m <- hlaModelFromObj(mobj)
                    on.exit(hlaClose(m))
                    hlaPredict(m, snp[,i], type=type, vote=vote, verbose=FALSE)
                } else
                    NULL
            },
            mobj=hlaModelToObj(object), snp=snp, type=type, vote=vote
        )

        # merge the results
        if (type %in% c("response", "dosage", "response+prob"))
        {
            res <- rv[[1L]]
            for (i in 2L:length(rv))
            {
                if (!is.null(rv[[i]]))
                    res <- hlaCombineAllele(res, rv[[i]])
            }
            res$value$sample.id <- geno.sampid
            if (!is.null(res$dosage))
                colnames(res$dosage) <- geno.sampid
            if (!is.null(res$postprob))
                colnames(res$postprob) <- geno.sampid
            NA.cnt <- sum(is.na(res$value$allele1) | is.na(res$value$allele2))
        } else {
            res <- rv[[1L]]
            for (i in 2L:length(rv))
            {
                if (!is.null(rv[[i]]))
                    res <- cbind(res, rv[[i]])
            }
            colnames(res) <- geno.sampid
            NA.cnt <- sum(colSums(res) <= 0L, na.rm=TRUE)
        }
    } 

    if (NA.cnt > 0L)
    {
        warning("No prediction output for ", NA.cnt, " individual",
            .plural(NA.cnt), " (possibly due to missing SNPs).",
            immediate.=TRUE)
    }
    # return
    res
}


#######################################################################
# Merge predictions by voting
#

hlaPredMerge <- function(..., weight=NULL, equivalence=NULL)
{
    # check "..."
    pdlist <- list(...)
    if (length(pdlist) <= 0L)
        stop("No object is passed to 'hlaPredMerge'.")
    for (i in 1:length(pdlist))
    {
        if (!inherits(pdlist[[i]], "hlaAlleleClass"))
        {
            stop("The object(s) passed to 'hlaPredMerge' ",
                "should be 'hlaAlleleClass'.")
        }
        if (is.null(pdlist[[i]]$postprob))
        {
            stop("The object(s) passed to 'hlaPredMerge' should have ",
                "a field of 'postprob' returned from ",
                "'hlaPredict(..., type=\"response+prob\", vote=\"majority\")'.")
        }
    }

    # check equivalence
    stopifnot(is.null(equivalence) | is.data.frame(equivalence))
    if (is.data.frame(equivalence))
    {
        if (ncol(equivalence) != 2L)
        {
            stop("'equivalence' should have two columns: ",
                "the first for new equivalent alleles, and ",
                "the second for the alleles possibly existed ",
                "in the object(s) passed to 'hlaPredMerge'.")
        }
    }

    # check locus and sample.id
    samp.id <- pdlist[[1L]]$value$sample.id
    locus <- pdlist[[1L]]$locus
    for (i in 1L:length(pdlist))
    {
        if (!identical(samp.id, pdlist[[i]]$value$sample.id))
            stop("The sample IDs should be the same.")
        if (!identical(locus, pdlist[[i]]$locus))
            stop("The locus should be the same.")
    }

    # check weight
    if (!is.null(weight))
    {
        stopifnot(is.numeric(weight) & is.vector(weight))
        if (length(pdlist) != length(weight))
            stop("Invalid 'weight'.")
        weight <- abs(weight)
        weight <- weight / sum(weight)
    } else {
        weight <- rep(1/length(pdlist), length(pdlist))
    }

    #############################################################
    # replace function
    replace <- function(allele)
    {
        if (!is.null(equivalence))
        {
            i <- match(allele, equivalence[, 2L])
            flag <- !is.na(i)
            i <- i[flag]
            allele[flag] <- equivalence[i, 1L]
        }
        allele
    }


    # all different alleles
    hla.allele <- NULL
    for (i in 1L:length(pdlist))
    {
        h <- unique(unlist(strsplit(rownames(pdlist[[i]]$postprob), "/")))
        hla.allele <- unique(c(hla.allele, replace(h)))
    }
    hla.allele <- hlaUniqueAllele(hla.allele)
    n.hla <- length(hla.allele)
    n.samp <- length(samp.id)

    prob <- matrix(0.0, nrow=n.hla*(n.hla+1L)/2, ncol=n.samp)
    m <- outer(hla.allele, hla.allele, function(x, y) paste(x, y, sep="/"))
    m <- m[lower.tri(m, diag=TRUE)]

    # for-loop
    for (i in 1L:length(pdlist))
    {
        p <- pdlist[[i]]$postprob
        h <- replace(unlist(strsplit(rownames(p), "/")))

        h1 <- h[seq(1L, length(h), 2L)]; h2 <- h[seq(2L, length(h), 2L)]
        j1 <- match(paste(h1, h2, sep="/"), m)
        j2 <- match(paste(h2, h1, sep="/"), m)
        j1[is.na(j1)] <- j2[is.na(j1)]

        # check
        stopifnot(!any(is.na(j1)))

        p <- p * weight[i]
        for (j in seq_len(length(j1)))
            prob[j1[j], ] <- prob[j1[j], ] + p[j, ]
    }
    colnames(prob) <- samp.id
    rownames(prob) <- m

    pb <- apply(prob, 2L, max)
    pt <- unlist(strsplit(m[apply(prob, 2L, which.max)], "/"))
    assembly <- pdlist[[1L]]$assembly
    if (is.null(assembly)) assembly <- "auto"

    rv <- hlaAllele(samp.id,
        H1 = pt[seq(2L, length(pt), 2L)],
        H2 = pt[seq(1L, length(pt), 2L)],
        locus = locus,
        locus.pos.start = pdlist[[1L]]$pos.start,
        locus.pos.end = pdlist[[1L]]$pos.end,
        prob = pb, na.rm = FALSE,
        assembly = assembly)
    rv$postprob <- prob
    rv
}


#######################################################################
# Summarize the "hlaAttrBagClass" object
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
    stopifnot(inherits(model, "hlaAttrBagClass"))

    # get a list of classifiers
    clr <- .Call(HIBAG_GetClassifierList, model$model, model$hla.allele)

    # output
    rv <- list(n.samp = model$n.samp, n.snp = model$n.snp,
        sample.id = model$sample.id, snp.id = model$snp.id,
        snp.position = model$snp.position, snp.allele = model$snp.allele,
        snp.allele.freq = model$snp.allele.freq,
        hla.locus = model$hla.locus,
        hla.allele = model$hla.allele, hla.freq = model$hla.freq,
        assembly = model$assembly,
        classifiers = clr,
        matching = model$matching,
        appendix = model$appendix)
    class(rv) <- "hlaAttrBagObj"
    rv
}


#######################################################################
# Combine two HIBAG models
#

hlaCombineModelObj <- function(obj1, obj2)
{
    # check
    stopifnot(inherits(obj1, "hlaAttrBagObj"))
    stopifnot(inherits(obj2, "hlaAttrBagObj"))
    stopifnot(identical(obj1$hla.locus, obj2$hla.locus))
    stopifnot(identical(obj1$snp.id, obj2$snp.id))
    stopifnot(identical(obj1$hla.allele, obj2$hla.allele))
    stopifnot(identical(obj1$assembly, obj2$assembly))

    samp.id <- unique(c(obj1$sample.id, obj2$sample.id))
    if (!is.null(obj1$appendix) | !is.null(obj2$appendix))
    {
        appendix <- list(
            platform =
                unique(c(obj1$appendix$platform, obj2$appendix$platform)),
            information =
                unique(c(obj1$appendix$information,
                    obj2$appendix$information)),
            warning =
                unique(c(obj1$appendix$warning, obj2$appendix$warning))
        )
    } else {
        appendix <- NULL
    }

    rv <- list(n.samp = length(samp.id), n.snp = obj1$n.snp,
        sample.id = samp.id, snp.id = obj1$snp.id,
        snp.position = obj1$snp.position, snp.allele = obj1$snp.allele,
        snp.allele.freq = (obj1$snp.allele.freq + obj2$snp.allele.freq)*0.5,
        hla.locus = obj1$hla.locus, hla.allele = obj1$hla.allele,
        hla.freq = (obj1$hla.freq + obj2$hla.freq)*0.5,
        assembly = obj1$assembly,
        classifiers = c(obj1$classifiers, obj2$classifiers),
        matching = c(obj1$matching, obj2$matching),
        appendix = appendix)
    if (identical(obj1$sample.id, obj2$sample.id) && !is.null(obj1$sample.id))
    {
        n1 <- length(obj1$classifiers)
        n2 <- length(obj2$classifiers)
        n <- n1 + n2
        rv$matching <- n1/n*obj1$matching + n2/n*obj2$matching
    }
    class(rv) <- "hlaAttrBagObj"
    rv
}


#######################################################################
# Get a HIBAG model with the top n individual classifiers
#

hlaSubModelObj <- function(obj, n)
{
    # check
    stopifnot(inherits(obj, "hlaAttrBagObj"))
    stopifnot(is.numeric(n), length(n)==1L)
    obj$classifiers <- obj$classifiers[1L:n]
    obj
}


#######################################################################
# Create a "hlaAttrBagClass" class from an R object
#

hlaModelFromObj <- function(obj)
{
    # check
    stopifnot(inherits(obj, "hlaAttrBagObj"))

    # create an attribute bagging object
    ABmodel <- .Call(HIBAG_New, obj$n.samp, obj$n.snp, length(obj$hla.allele))

    # add individual classifiers
    for (tree in obj$classifiers)
    {
        hla <- match(tree$haplos$hla, obj$hla.allele) - 1L
        if (any(is.na(hla)))
            stop("Invalid HLA alleles in the individual classifier.")
        if (is.null(tree$samp.num))
            snum <- rep.int(1L, obj$n.samp)
        else
            snum <- as.integer(tree$samp.num)

        # create a new classifier
        .Call(HIBAG_NewClassifierHaplo, ABmodel, as.integer(tree$snpidx - 1L),
            snum, as.double(tree$haplos$freq), hla,
            as.character(tree$haplos$haplo),
            tree$outofbag.acc)
    }

    # output
    rv <- list(n.samp = obj$n.samp, n.snp = obj$n.snp,
        sample.id = obj$sample.id, snp.id = obj$snp.id,
        snp.position = obj$snp.position, snp.allele = obj$snp.allele,
        snp.allele.freq = obj$snp.allele.freq,
        hla.locus = obj$hla.locus, hla.allele = obj$hla.allele,
        hla.freq = obj$hla.freq,
        assembly = as.character(obj$assembly)[1L],
        model = ABmodel,
        matching = obj$matching,
        appendix = obj$appendix)
    if (is.na(rv$assembly)) rv$assembly <- "unknown"

    class(rv) <- "hlaAttrBagClass"
    rv
}


#######################################################################
# Summarize the "hlaAttrBagObj" object
#

summary.hlaAttrBagObj <- function(object, show=TRUE, ...)
{
    # check
    stopifnot(inherits(object, "hlaAttrBagObj"))
    obj <- object

    if (show)
    {
        cat("Gene: ", .hla_gene_name_string(obj$hla.locus), "\n", sep="")
        cat("Training dataset:", obj$n.samp, "samples X",
            length(obj$snp.id), "SNPs\n")
        cat("    # of HLA alleles: ", length(obj$hla.allele), "\n", sep="")
    }

    # summarize ...
    snpset <- NULL
    outofbag.acc <- rep(NaN, length(obj$classifiers))
    numsnp <- rep(NA, length(obj$classifiers))
    numhaplo <- rep(NA, length(obj$classifiers))
    snp.hist <- rep(0, length(obj$snp.id))
    for (i in 1L:length(obj$classifiers))
    {
        outofbag.acc[i] <- obj$classifiers[[i]]$outofbag.acc
        numsnp[i] <- length(obj$classifiers[[i]]$snpidx)
        numhaplo[i] <- length(obj$classifiers[[i]]$haplos$hla)
        snp.hist[obj$classifiers[[i]]$snpidx] <-
            snp.hist[obj$classifiers[[i]]$snpidx] + 1L
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
        cat("    # of individual classifiers: ", length(obj$classifiers),
            "\n", sep="")
        cat("    total # of SNPs used: ", length(snpset), "\n", sep="")
        cat("    avg. # of SNPs in an individual classifier:",
            sprintf("%0.2f\n        (sd: %0.2f, min: %d, max: %d, median: %0.2f)\n",
            mean(numsnp), sd(numsnp), min(numsnp), max(numsnp),
            median(numsnp)))
        cat("    avg. # of haplotypes in an individual classifier:",
            sprintf("%0.2f\n        (sd: %0.2f, min: %d, max: %d, median: %0.2f)\n",
            mean(numhaplo), sd(numhaplo), min(numhaplo), max(numhaplo),
            median(numhaplo)))
        cat("    avg. out-of-bag accuracy:",
            sprintf("%0.2f%%\n        (sd: %0.2f%%, min: %0.2f%%, max: %0.2f%%, median: %0.2f%%)\n",
            mean(outofbag.acc), sd(outofbag.acc), min(outofbag.acc),
            max(outofbag.acc), median(outofbag.acc)))

        p <- obj$matching
        if (!is.null(p))
        {
            cat("Matching proportion:\n")
            .printMatching(p)
        }

        if (is.null(obj$assembly))
            cat("Genome assembly: unknown\n")
        else
            cat("Genome assembly: ", obj$assembly, "\n", sep="")

        if (!is.null(obj$appendix$platform))
            cat("Platform:", obj$appendix$platform, "\n")
        if (!is.null(obj$appendix$information))
            cat("Information:", obj$appendix$information, "\n")
        if (!is.null(obj$appendix$warning))
            message(obj$appendix$warning)
    }

    rv <- list(num.classifier = length(obj$classifiers),
        num.snp = length(snpset),
        snp.id = obj$snp.id, snp.position = obj$snp.position,
        snp.hist = snp.hist, info = info)
    invisible(rv)
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

    stopifnot(is.numeric(call.threshold), length(call.threshold)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)


    ######################################################
    # initialize ...
    if (inherits(model, "hlaAttrBagClass"))
    {
        model <- hlaModelToObj(model)
        if (verbose)
            print(model)
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
    for (i in 1L:nclass)
    {
        mx <- model
        mx$classifiers <- mx$classifiers[i]
        s <- mx$classifiers[[1L]]$samp.num
        if (is.null(s))
            stop("There is no bootstrap sample index.")

        tmp.model <- hlaModelFromObj(mx)
        tmp.geno <- geno[, s == 0L]
        v <- hlaPredict(tmp.model, tmp.geno, verbose=FALSE)
        hlaClose(tmp.model)
        v$value$sample.id <- mx$sample.id[s == 0L]
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
    m.max <- apply(rv, 2L, max); m.idx <- apply(rv, 2L, which.max)
    s <- rownames(ans$confusion)[m.idx]; s[m.max<=0] <- NA
    p <- m.max / apply(rv, 2L, sum)

    # output
    ans$detail <- cbind(ans$detailhead, ans$detail,
        miscall=s, miscall.prop=p, stringsAsFactors=FALSE)
    ans$detailhead <- NULL
    ans$n.detail <- NULL
    ans
}


##########################################################################
##########################################################################
#
# Linkage Disequilibrium
#

##########################################################################
# Calculate linkage disequilibrium between HLA locus and SNP markers
#

hlaGenoLD <- function(hla, geno)
{
    # check
    stopifnot(inherits(hla, "hlaAlleleClass"))
    if (inherits(geno, "hlaSNPGenoClass"))
    {
        stopifnot(dim(hla$value)[1L] == length(geno$sample.id))
        if (any(hla$value$sample.id != geno$sample.id))
        {
            hla <- hlaAlleleSubset(hla, samp.sel =
                match(geno$sample.id, hla$value$sample.id))
        }
        geno <- geno$genotype
    } else if (is.matrix(geno))
    {
        stopifnot(is.numeric(geno))
        stopifnot(dim(hla$value)[1L] == dim(geno)[2L])
    } else if (is.vector(geno))
    {
        stopifnot(is.numeric(geno))
        stopifnot(dim(hla$value)[1L] == length(geno))
        geno <- matrix(geno, ncol=1L)
    } else {
        stop("geno should be `hlaSNPGenoClass', a vector or a matrix.")
    }

    # HLA alleles indicators
    alleles <- unique(c(hla$value$allele1, hla$value$allele2))
    alleles <- alleles[order(alleles)]
    allele.mat <- matrix(0L,
        nrow=length(hla$value$allele1), ncol=length(alleles))
    for (i in 1L:length(alleles))
    {
        allele.mat[, i] <- (hla$value$allele1==alleles[i]) +
            (hla$value$allele2==alleles[i])
    }

    apply(geno, 1L,
        function(x, allele.mat) {
            suppressWarnings(
                mean(
                    cor(x, allele.mat, use="pairwise.complete.obs")^2,
                    na.rm=TRUE
                )
            )
        },
        allele.mat=allele.mat)
}


##########################################################################
# SNP LD visualization
#

hlaLDMatrix <- function(geno, loci=NULL, maf=0.01, assembly="auto",
    draw=TRUE, verbose=TRUE)
{
    # check
    stopifnot(inherits(geno, "hlaSNPGenoClass"))
    stopifnot(is.numeric(maf), length(maf)==1L)
    stopifnot(is.null(loci) | is.vector(loci))
    stopifnot(is.logical(draw), length(draw)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # maf filter
    if (is.na(maf)) maf <- 0
    if (maf > 0)
    {
        af <- rowMeans(geno$genotype, na.rm=TRUE) * 0.5
        af <- pmin(af, 1-af)
        xx <- af >= maf
        xx[is.na(xx)] <- FALSE
        if (sum(xx) < length(af))
        {
            if (verbose)
            {
                cat("MAF filter (>=", maf, "), excluding ",
                    length(af)-sum(xx), " SNP(s)\n", sep="")
            }
            geno <- hlaGenoSubset(geno, snp.sel=xx)
        }
    }

    if (!is.null(loci))
    {
        if (assembly=="auto") assembly <- geno$assembly
        assembly <- .hla_assembly(assembly)
        info <- hlaLociInfo(assembly)
        if (!all(loci %in% rownames(info)))
        {
            stop("'loci' should be one of ", paste(rownames(info),
                collapse=", "))
        }
        info <- info[loci, ]
        st <- sapply(info$start, function(p) {
            mean(c(which(p <= geno$snp.position)[1L],
            rev(which(p >= geno$snp.position))[1L])) })
        ed <- sapply(info$end, function(p) {
            mean(c(which(p <= geno$snp.position)[1L],
            rev(which(p >= geno$snp.position))[1L])) })
    }

    # calculate the LD matrix
    ld <- suppressWarnings(cor(t(geno$genotype), use="na.or.complete")^2)

    if (isTRUE(draw))
    {
        Var1 <- Var2 <- value <- NULL
        dat <- reshape2::melt(ld)
        p <- ggplot2::ggplot(dat, ggplot2::aes(x=Var2, y=Var1)) +
            ggplot2::geom_raster(ggplot2::aes(fill=value)) +
            ggplot2::scale_fill_gradient2(low="grey95", mid="orange", high="red", midpoint=0.5)
        rg <- range(geno$snp.position)
        p <- p + ggplot2::labs(
            x=sprintf("%d SNPs in total, MAF >= %g", length(geno$snp.id), maf),
            y="SNP index", title=expression("Linkage Disequilibrium r"^2))
        ii <- c(0, 0.25, 0.5, 0.75, 1)
        p <- p + ggplot2::scale_x_continuous(
            breaks=unname(quantile(1:length(geno$snp.position), ii)),
            labels=sprintf("%.3fMb", c(quantile(geno$snp.position, ii)/10^6)))
        for (i in seq_along(loci))
        {
            md <- (st[i]+ed[i])*0.5
            if (is.finite(md))
            {
                p <- p +
                    ggplot2::geom_vline(xintercept=st[i], size=0.25, colour="blue", alpha=0.5) +
                    ggplot2::geom_vline(xintercept=ed[i], size=0.25, colour="blue", alpha=0.5) +
                    ggplot2::geom_vline(xintercept=md, linetype=2, colour="blue") +
                    ggplot2::annotate("label", x=md,
                        y=length(geno$snp.position), vjust="inward",
                        label=ifelse(grepl("^KIR", loci[i]), loci[i],
                            paste("HLA", loci[i], sep="-")))
            }
        }
        return(p)
    } else {
        return(ld)
    }
}


##########################################################################
# Calculate the distances among different HLA alleles
#

hlaDistance <- function(model)
{
    # check
    if (inherits(model, "hlaAttrBagClass"))
        model <- hlaModelToObj(model)

    # the total number of classifiers
    n <- length(model$classifiers)
    lst <- vector("list", n)
    hla <- model$hla.allele
    num <- matrix(0L, nrow=length(hla), ncol=length(hla))
    for (i in seq_len(n))
    {
        z <- model$classifiers[[i]]$haplos
        m <- .Call(HIBAG_Distance, length(hla), match(z$hla, hla),
            z$freq, z$haplo)
        num <- num + !is.na(m)
        m[is.na(m)] <- 0
        lst[[i]] <- m
    }

    # output
    rv <- Reduce("+", lst) / num
    colnames(rv) <- rownames(rv) <- hla
    rv
}


##########################################################################
##########################################################################
#
# Visualization
#

##########################################################################
# Plot an attribute bagging model
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
    invisible()
}



##########################################################################
# Plot an attribute bagging model
#

plot.hlaAttrBagObj <- function(x, snp.col="gray33", snp.pch=1, snp.sz=1,
    locus.col="blue", locus.lty=1L, locus.lty2=2L, addplot=NULL,
    assembly="auto", ...)
{
    # check
    stopifnot(inherits(x, "hlaAttrBagObj"))
    assembly <- match.arg(assembly)

    # the starting and ending positions of HLA locus
    if (assembly == "auto")
    {
        if (!is.null(x$assembly))
            assembly <- x$assembly
    }
    info <- hlaLociInfo(assembly)
    pos.start <- info[x$hla.locus, "start"] / 1000L
    pos.end <- info[x$hla.locus, "end"] / 1000L

    # summary of the attribute bagging model
    desp <- summary(x, show=FALSE)

    # draw
    dat <- data.frame(pos=x$snp.position/1000L, ht=desp$snp.hist)
    pos <- ht <- NULL

    if (is.null(addplot))
    {
        p <- ggplot2::ggplot(dat, ggplot2::aes(x=pos, y=ht)) +
            ggplot2::xlab("SNP Position (KB)") +
            ggplot2::ylab("Frequency of Use")
    } else {
        p <- addplot
    }
    p <- p +
        ggplot2::geom_vline(xintercept=pos.start, colour=locus.col,
            linetype=locus.lty, size=0.25, alpha=0.5) +
        ggplot2::geom_vline(xintercept=pos.end, colour=locus.col,
            linetype=locus.lty, size=0.25, alpha=0.5) +
        ggplot2::geom_vline(xintercept=(pos.start+pos.end)*0.5,
            linetype=locus.lty2, colour=locus.col) +
        ggplot2::annotate("label", x=(pos.start + pos.end)/2,
            y=max(dat$ht), vjust="inward",
            label=ifelse(grepl("^KIR", x$hla.locus), x$hla.locus,
                paste("HLA", x$hla.locus, sep="-")))
    if (is.null(addplot))
    {
        p <- p + ggplot2::geom_point(color=snp.col, shape=snp.pch, size=snp.sz)
    } else {
        p <- p + ggplot2::geom_point(ggplot2::aes(x=pos, y=ht), dat,
            color=snp.col, shape=snp.pch, size=snp.sz)
    }
    p
}

print.hlaAttrBagObj <- function(x, ...)
{
    summary(x)
    invisible()
}



#######################################################################
# Export stardard R library function(s)
#######################################################################

hlaSetKernelTarget <- function(cpu=c("max", "auto.avx2", "base",
    "sse2", "sse4", "avx", "avx2", "avx512f", "avx512bw"))
{
    cpu <- match.arg(cpu)
    .Call(HIBAG_Kernel_SetTarget, cpu)
    .Call(HIBAG_Kernel_Version)[[2L]]
}


.onAttach <- function(lib, pkg)
{
    # get version and CPU information
    info <- .Call(HIBAG_Kernel_Version)

    # information
    packageStartupMessage(
        "HIBAG (HLA Genotype Imputation with Attribute Bagging)")
    packageStartupMessage(
        sprintf("Kernel Version: v%d.%d (%s)", info[[1L]][1L], info[[1L]][2L],
        info[[2L]][1L]))
    if (is.na(info[[3L]]))
        packageStartupMessage("No Intel Threading Building Blocks (TBB)")
    TRUE
}

.hibag_data_list <- list(
	data.frame = "data.frame",
	clr_nm = c("samp.num", "haplos", "snpidx", "outofbag.acc"),
	clr_haplo_nm = c("freq", "hla", "haplo"))

.onLoad <- function(lib, pkg)
{
    .Call(HIBAG_Init, .hibag_data_list)
    TRUE
}
