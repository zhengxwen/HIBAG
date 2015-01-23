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


##########################################################################
#
# Attribute Bagging method -- HIBAG algorithm
#

##########################################################################
# To fit an attribute bagging model for predicting
#

hlaAttrBagging <- function(hla, snp, nclassifier=100,
    mtry=c("sqrt", "all", "one"), prune=TRUE, rm.na=TRUE,
    verbose=TRUE, verbose.detail=FALSE)
{
    # check
    stopifnot(inherits(hla, "hlaAlleleClass"))
    stopifnot(inherits(snp, "hlaSNPGenoClass"))
    stopifnot(is.character(mtry) | is.numeric(mtry))
    stopifnot(is.logical(verbose))
    stopifnot(is.logical(verbose.detail))
    if (verbose.detail) verbose <- TRUE

    # get the common samples
    samp.id <- intersect(hla$value$sample.id, snp$sample.id)

    # hla types
    samp.flag <- match(samp.id, hla$value$sample.id)
    hla.allele1 <- hla$value$allele1[samp.flag]
    hla.allele2 <- hla$value$allele2[samp.flag]
    if (rm.na)
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
    storage.mode(snp.geno) <- "integer"

    tmp.snp.id <- snp$snp.id
    tmp.snp.position <- snp$snp.position
    tmp.snp.allele <- snp$snp.allele

    # remove mono-SNPs
    snpsel <- rowMeans(snp.geno, na.rm=TRUE)
    snpsel[!is.finite(snpsel)] <- 0
    snpsel <- (0 < snpsel) & (snpsel < 2)
    if (sum(!snpsel) > 0)
    {
        snp.geno <- snp.geno[snpsel, ]
        if (verbose)
        {
            a <- sum(!snpsel)
            cat(sprintf("%s monomorphic SNP%s ha%s been removed.\n",
                a, if (a>1) "s" else "", if (a>1) "ve" else "s"))
        }
        tmp.snp.id <- tmp.snp.id[snpsel]
        tmp.snp.position <- tmp.snp.position[snpsel]
        tmp.snp.allele <- tmp.snp.allele[snpsel]
    }

    if (length(samp.id) <= 0)
        stop("There is no common sample between 'hla' and 'snp'.")
    if (length(dim(snp.geno)[1]) <= 0)
        stop("There is no valid SNP markers.")


    ###################################################################
    # initialize ...

    n.snp <- dim(snp.geno)[1]      # Num. of SNPs
    n.samp <- dim(snp.geno)[2]     # Num. of samples
    HUA <- hlaUniqueAllele(c(hla.allele1, hla.allele2))
    H <- factor(match(c(hla.allele1, hla.allele2), HUA))
    levels(H) <- HUA
    n.hla <- nlevels(H)
    H1 <- as.integer(H[1:n.samp]) - as.integer(1)
    H2 <- as.integer(H[(n.samp+1):(2*n.samp)]) - as.integer(1)

    # create an attribute bagging object
    rv <- .C(HIBAG_Training, n.snp, n.samp, snp.geno, n.hla,
        H1, H2, AB=integer(1), err=integer(1), NAOK=TRUE)
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
        cat("Build a HIBAG model with", nclassifier,
            "individual classifier(s):\n")
        cat("# of SNPs randomly sampled as candidates for each selection: ",
            mtry, "\n", sep="")
        cat("# of SNPs: ", n.snp, ", # of samples: ", n.samp, "\n", sep="")
        cat("# of unique HLA alleles: ", n.hla, "\n", sep="")
    }


    ###################################################################
    # training ...
    # add new individual classifers
    rv <- .C(HIBAG_NewClassifiers, ABmodel, as.integer(nclassifier),
        as.integer(mtry), as.logical(prune), verbose, verbose.detail,
        err=integer(1), NAOK=TRUE)
    if (rv$err != 0) stop(hlaErrMsg())

    # output
    rv <- list(n.samp = n.samp, n.snp = n.snp, sample.id = samp.id,
        snp.id = tmp.snp.id, snp.position = tmp.snp.position,
        snp.allele = tmp.snp.allele,
        snp.allele.freq = 0.5*rowMeans(snp.geno, na.rm=TRUE),
        hla.locus = hla$locus, hla.allele = levels(H),
        hla.freq = prop.table(table(H)),
        assembly = as.character(snp$assembly)[1],
        model = ABmodel,
        appendix = list())
    if (is.na(rv$assembly)) rv$assembly <- "unknown"

    class(rv) <- "hlaAttrBagClass"
    rv
}


##########################################################################
# To fit an attribute bagging model for predicting
#

hlaParallelAttrBagging <- function(cl, hla, snp, auto.save="",
    nclassifier=100, mtry=c("sqrt", "all", "one"), prune=TRUE, rm.na=TRUE,
    stop.cluster=FALSE, verbose=TRUE)
{
    # check
    stopifnot(is.null(cl) | inherits(cl, "cluster"))
    stopifnot(inherits(hla, "hlaAlleleClass"))
    stopifnot(inherits(snp, "hlaSNPGenoClass"))

    stopifnot(is.character(auto.save) & (length(auto.save)==1))
    stopifnot(is.numeric(nclassifier))
    stopifnot(is.character(mtry) | is.numeric(mtry))
    stopifnot(is.logical(prune))
    stopifnot(is.logical(rm.na))
    stopifnot(is.logical(stop.cluster))
    stopifnot(is.logical(verbose))

    if (!is.null(cl))
    {
        if (!require(parallel, warn.conflicts=FALSE))
            stop("The `parallel' package should be installed.")
    }

    if (verbose)
    {
        if (!is.null(cl))
        {
            cat("Build a HIBAG model of", sprintf(
                "%d individual classifier%s in parallel with %d node%s:\n",
                nclassifier, if (nclassifier>1) "s" else "",
                length(cl), if (length(cl)>1) "s" else ""
            ))
        } else {
            cat(sprintf(
                "Build a HIBAG model of %d individual classifier%s:\n",
                nclassifier, if (nclassifier>1) "s" else ""
            ))
        }
        if (auto.save != "")
            cat("The model is autosaved in '", auto.save, "'.\n", sep="")
    }

    # set random number
    if (!is.null(cl))
    {
        RNGkind("L'Ecuyer-CMRG")
        rand <- eval(parse(text=".Random.seed"))
        parallel::clusterSetRNGStream(cl)
    }

    ans <- local({
        total <- 0L

        .DynamicClusterCall(cl,
            fun = function(job, hla, snp, mtry, prune, rm.na)
            {
                eval(parse(text="library(HIBAG)"))
                model <- hlaAttrBagging(hla=hla, snp=snp, nclassifier=1,
                    mtry=mtry, prune=prune, rm.na=rm.na,
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
                    save(mobj, file=auto.save)
                if (verbose & !is.null(mobj))
                {
                    z <- summary(mobj, show=FALSE)
                    cat("  --  average out-of-bag accuracy:", sprintf(
                        "%0.2f%%, sd: %0.2f%%, min: %0.2f%%, max: %0.2f%%\n",
                        z$info["accuracy", "Mean"], z$info["accuracy", "SD"],
                        z$info["accuracy", "Min"], z$info["accuracy", "Max"]))
                }
                mobj
            },
            msg.fn = function(job, obj)
            {
                if (verbose)
                {
                    z <- summary(obj, show=FALSE)
                    eval(parse(text="total <<- total + 1L"))
                    cat(date(), sprintf(
    ", %4d, job %3d, # of SNPs: %g, # of haplotypes: %g, accuracy: %0.1f%%\n",
                        total, as.integer(job), z$info["num.snp", "Mean"],
                        z$info["num.haplo", "Mean"],
                        z$info["accuracy", "Mean"]), sep="")
                }
            },
            n = nclassifier, stop.cluster = stop.cluster,
            hla=hla, snp=snp, mtry=mtry, prune=prune, rm.na=rm.na
        )
    })

    if (!is.null(cl) & !stop.cluster)
    {
        parallel::nextRNGStream(rand)
        parallel::nextRNGSubStream(rand)
    }

    # return
    if (auto.save == "")
        hlaModelFromObj(ans)
    else
        invisible()
}


##########################################################################
# To fit an attribute bagging model for predicting
#

hlaClose <- function(model)
{
    # check
    stopifnot(inherits(model, "hlaAttrBagClass"))

    # class handler
    rv <- .C(HIBAG_Close, model$model, err=integer(1), NAOK=TRUE)
    if (rv$err != 0) stop(hlaErrMsg())

    # output
    invisible()
}


#######################################################################
# Predict HLA types from unphased SNP data
#

predict.hlaAttrBagClass <- function(object, snp, cl=NULL,
    type=c("response", "prob", "response+prob"), vote=c("prob", "majority"),
    allele.check=TRUE, match.type=c("RefSNP+Position", "RefSNP", "Position"),
    same.strand=FALSE, verbose=TRUE, ...)
{
    # check
    stopifnot(inherits(object, "hlaAttrBagClass"))
    stopifnot(is.null(cl) | inherits(cl, "cluster"))
    stopifnot(is.logical(allele.check))
    stopifnot(is.logical(same.strand))
    stopifnot(is.logical(verbose))

    type <- match.arg(type)
    vote <- match.arg(vote)
    match.type <- match.arg(match.type)
    vote_method <- match(vote, c("prob", "majority"))
    if (!is.null(cl))
    {
        if (!require(parallel, warn.conflicts=FALSE))
            stop("The `parallel' package should be installed.")
        if (length(cl) <= 1) cl <- NULL
    }

    # if warning
    if (!is.null(object$appendix$warning))
    {
        message(object$appendix$warning)
        warning(object$appendix$warning)
    }

    if (verbose)
    {
        # call, get the number of classifiers
        rv <- .C(HIBAG_GetNumClassifiers, object$model, CNum = integer(1),
            err=integer(1), NAOK=TRUE)
        if (rv$err != 0) stop(hlaErrMsg())

        if (rv$CNum > 1) { s <- "s" } else { s <- "" }
        cat(sprintf(
"HIBAG model: %d individual classifier%s, %d SNPs, %d unique HLA alleles.\n",
            rv$CNum, s, length(object$snp.id), length(object$hla.allele)))

        if (vote_method == 1)
        {
            cat("Predicting based on the averaged posterior probabilities",
                "from all individual classifiers.\n")
        } else
            cat("Predicting by voting from all individual classifiers.\n")

        if (!is.null(cl))
        {
            cat(sprintf(
                "Run in parallel with %d computing node%s.\n",
                length(cl), if (length(cl)>1) "s" else ""
            ))
        }
    }

    if (!inherits(snp, "hlaSNPGenoClass"))
    {
        # it should be a vector or a matrix
        stopifnot(is.numeric(snp))
        stopifnot(is.vector(snp) | is.matrix(snp))

        if (is.vector(snp))
        {
            stopifnot(length(snp) == object$n.snp)
            snp <- matrix(snp, ncol=1)
        } else {
            stopifnot(nrow(snp) == object$n.snp)
        }
        geno.sampid <- 1:ncol(snp)
        assembly <- "auto-silent"

    } else {

        # a 'hlaSNPGenoClass' object

        ##################################################
        # check assembly first

        model.assembly <- as.character(object$assembly)[1]
        if (is.na(model.assembly))
            model.assembly <- "unknown"
        geno.assembly <- as.character(snp$assembly)[1]
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
                if (verbose)
                    message("The human genome references do not match!")
                warning("The human genome references do not match! ",
                    refstr, ".")
                assembly <- model.assembly
            }
        } else {
            if (model.assembly != "unknown")
                assembly <- model.assembly
            else
                assembly <- "auto"
        }

        ##################################################

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
            snp.allele <- snp$snp.allele
            snp.allele[is.na(snp.allele)] <- ""

            # tmp variable
            tmp <- list(genotype = snp$genotype[snp.sel,],
                sample.id = snp$sample.id,
                snp.id = object$snp.id, snp.position = object$snp.position,
                snp.allele = snp.allele[snp.sel],
                assembly = snp$assembly)

            flag <- is.na(tmp$snp.allele)
            tmp$snp.allele[flag] <- object$snp.allele[
                match(tmp$snp.id[flag], object$snp.id)]

            if (is.vector(tmp$genotype))
                tmp$genotype <- matrix(tmp$genotype, ncol=1)

            class(tmp) <- "hlaSNPGenoClass"

            # the total number of missing snp
            missing.cnt <- sum(flag)

            # verbose
            if (verbose)
            {
                if (missing.cnt > 0)
                {
                    cat(sprintf("There %s %d missing SNP%s (%0.1f%%).\n",
                        if (missing.cnt > 1) "are" else "is",
                        missing.cnt,
                        if (missing.cnt > 1) "s" else "",
                        100*missing.cnt/length(obj.id)))
                }
            }

            # try alternative matching if possible
            if (match.type == "RefSNP+Position")
            {
                # match by RefSNP IDs
                s1 <- hlaSNPID(object, "RefSNP")
                s2 <- hlaSNPID(snp, "RefSNP")
                mcnt.id <- length(s1) - length(intersect(s1, s2))

                # match by positions
                s1 <- hlaSNPID(object, "Position")
                s2 <- hlaSNPID(snp, "Position")
                mcnt.pos <- length(s1) - length(intersect(s1, s2))

                if (((mcnt.id < missing.cnt) |
                    (mcnt.pos < missing.cnt)) & verbose)
                {
                    message(
                        "Hint:\n",
                        "The current SNP matching requires both RefSNP IDs ",
                        "and positions, and lower missing fraction(s) could ",
                        "be gained:\n",
                        sprintf("\t%0.1f%% by matching RefSNP IDs only, ",
                            100*mcnt.id/length(s1)),
                        "call 'predict(..., match.type=\"RefSNP\")'\n",
                        sprintf("\t%0.1f%% by matching positions only, ",
                            100*mcnt.pos/length(s1)),
                        "call 'predict(..., match.type=\"Position\")'\n",
                        "Any concern about SNP mismatching should be emailed ",
                        "to the genotyping platform provider.")

                    if (!is.null(object$appendix$platform))
                    {
                        message("The supported platform(s): ",
                            object$appendix$platform)
                    }
                }
            }

            if (missing.cnt == length(obj.id))
            {
                stop("There is no overlapping of SNPs!")
            } else if (missing.cnt > 0.5*length(obj.id))
            {
                warning("More than 50% of SNPs are missing!")
            }

            # switch
            if (allele.check)
            {
                snp <- hlaGenoSwitchStrand(tmp, object, match.type,
                    same.strand, verbose)$genotype
            } else {
                snp <- snp$genotype
            }
        }
    }

    # check
    if (dim(snp)[1] != object$n.snp)
        stop("Internal error! Please contact the author.")

    # initialize ...
    n.samp <- dim(snp)[2]
    n.hla <- length(object$hla.allele)
    if (verbose)
        cat(sprintf("The number of samples: %d.\n", n.samp))

    # parallel units
    if (is.null(cl))
    {
        # to predict HLA types
        if (type %in% c("response", "response+prob"))
        {
            # the best-guess prediction

            if (type == "response")
            {
                rv <- .C(HIBAG_Predict_Resp, object$model, as.integer(snp),
                    n.samp, as.integer(vote_method), as.logical(verbose),
                    H1=integer(n.samp), H2=integer(n.samp),
                    prob=double(n.samp), err=integer(1), NAOK=TRUE)
            } else {
                rv <- .C(HIBAG_Predict_Resp_Prob, object$model,
                    as.integer(snp), n.samp, as.integer(vote_method),
                    as.logical(verbose),
                    H1=integer(n.samp), H2=integer(n.samp),
                    prob=double(n.samp),
                    postprob=matrix(NaN, nrow=n.hla*(n.hla+1)/2, ncol=n.samp),
                    err=integer(1), NAOK=TRUE)
            }
            if (rv$err != 0) stop(hlaErrMsg())

            res <- hlaAllele(geno.sampid,
                H1=object$hla.allele[rv$H1 + 1],
                H2=object$hla.allele[rv$H2 + 1],
                locus = object$hla.locus, prob = rv$prob, na.rm = FALSE,
                assembly = assembly)
            if (!is.null(rv$postprob))
            {
                res$postprob <- rv$postprob
                colnames(res$postprob) <- geno.sampid
                m <- outer(object$hla.allele, object$hla.allele,
                    function(x, y) paste(x, y, sep="/"))
                rownames(res$postprob) <- m[lower.tri(m, diag=TRUE)]
            }

            NA.cnt <- sum(is.na(res$value$allele1) | is.na(res$value$allele2))

        } else {
            # all probabilites

            rv <- .C(HIBAG_Predict_Prob, object$model, as.integer(snp),
                n.samp, as.integer(vote_method), as.logical(verbose),
                prob=matrix(NaN, nrow=n.hla*(n.hla+1)/2, ncol=n.samp),
                err=integer(1), NAOK=TRUE)
            if (rv$err != 0) stop(hlaErrMsg())

            res <- rv$prob
            colnames(res) <- geno.sampid
            m <- outer(object$hla.allele, object$hla.allele,
                function(x, y) paste(x, y, sep="/"))
            rownames(res) <- m[lower.tri(m, diag=TRUE)]

            NA.cnt <- sum(colSums(res) <= 0)
        }
    } else {

        # in parallel
        rv <- parallel::clusterApply(cl=cl,
            parallel::splitIndices(n.samp, length(cl)),
            fun = function(idx, mobj, snp, type, vote)
            {
                if (length(idx) > 0)
                {
                    eval(parse(text="library(HIBAG)"))
                    m <- hlaModelFromObj(mobj)
                    pd <- predict(m, snp[,idx], type=type, vote=vote,
                        verbose=FALSE)
                    hlaClose(m)
                    pd
                } else
                    NULL
            },
            mobj=hlaModelToObj(object), snp=snp, type=type, vote=vote
        )

        if (type %in% c("response", "response+prob"))
        {
            res <- rv[[1]]
            for (i in 2:length(rv))
            {
                if (!is.null(rv[[i]]))
                    res <- hlaCombineAllele(res, rv[[i]])
            }
            res$value$sample.id <- geno.sampid
            if (!is.null(res$postprob))
                colnames(res$postprob) <- geno.sampid
            NA.cnt <- sum(is.na(res$value$allele1) | is.na(res$value$allele2))
        } else {
            res <- rv[[1]]
            for (i in 2:length(rv))
            {
                if (!is.null(rv[[i]]))
                    res <- cbind(res, rv[[i]])
            }
            colnames(res) <- geno.sampid
            NA.cnt <- sum(colSums(res) <= 0)
        }
    } 

    if (NA.cnt > 0)
    {
        if (NA.cnt > 1) s <- "s" else s <- ""
        warning(sprintf(
            "No prediction output%s for %d individual%s ",
            s, NA.cnt, s),
            "(possibly due to missing SNPs.)")
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
    if (length(pdlist) <= 0)
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
                "'predict(..., type=\"response+prob\", vote=\"majority\")'.")
        }
    }

    # check equivalence
    stopifnot(is.null(equivalence) | is.data.frame(equivalence))
    if (is.data.frame(equivalence))
    {
        if (ncol(equivalence) != 2)
        {
            stop("'equivalence' should have two columns: ",
                "the first for new equivalent alleles, and ",
                "the second for the alleles possibly existed ",
                "in the object(s) passed to 'hlaPredMerge'.")
        }
    }

    # check locus and sample.id
    samp.id <- pdlist[[1]]$value$sample.id
    locus <- pdlist[[1]]$locus
    for (i in 1:length(pdlist))
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
            i <- match(allele, equivalence[, 2])
            flag <- !is.na(i)
            i <- i[flag]
            allele[flag] <- equivalence[i, 1]
        }
        allele
    }


    # all different alleles
    hla.allele <- NULL
    for (i in 1:length(pdlist))
    {
        h <- unique(unlist(strsplit(rownames(pdlist[[i]]$postprob), "/")))
        hla.allele <- unique(c(hla.allele, replace(h)))
    }
    hla.allele <- hlaUniqueAllele(hla.allele)
    n.hla <- length(hla.allele)
    n.samp <- length(samp.id)

    prob <- matrix(0.0, nrow=n.hla*(n.hla+1)/2, ncol=n.samp)
    m <- outer(hla.allele, hla.allele, function(x, y) paste(x, y, sep="/"))
    m <- m[lower.tri(m, diag=TRUE)]

    # for-loop
    for (i in 1:length(pdlist))
    {
        p <- pdlist[[i]]$postprob
        h <- replace(unlist(strsplit(rownames(p), "/")))

        h1 <- h[seq(1, length(h), 2)]; h2 <- h[seq(2, length(h), 2)]
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

    pb <- apply(prob, 2, max)
    pt <- unlist(strsplit(m[apply(prob, 2, which.max)], "/"))
    assembly <- pdlist[[1]]$assembly
    if (is.null(assembly)) assembly <- "auto"

    rv <- hlaAllele(samp.id,
        H1 = pt[seq(2, length(pt), 2)],
        H2 = pt[seq(1, length(pt), 2)],
        locus = locus,
        locus.pos.start = pdlist[[1]]$pos.start,
        locus.pos.end = pdlist[[1]]$pos.end,
        prob = pb, na.rm = FALSE,
        assembly = assembly)
    rv$postprob <- prob
    rv
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
    stopifnot(inherits(model, "hlaAttrBagClass"))

    # call, get the number of classifiers
    rv <- .C(HIBAG_GetNumClassifiers, model$model, CNum = integer(1),
        err=integer(1), NAOK=TRUE)
    if (rv$err != 0) stop(hlaErrMsg())

    # for each tree
    res <- vector("list", rv$CNum)
    for (i in 1:length(res))
    {
        # call, get the number of haplotypes
        rv <- .C(HIBAG_Idv_GetNumHaplo, model$model, as.integer(i),
            NumHaplo = integer(1), NumSNP = integer(1),
            err=integer(1), NAOK=TRUE)
        if (rv$err != 0) stop(hlaErrMsg())

        # call, get freq. and haplotypes
        rv <- .C(HIBAG_Classifier_GetHaplos, model$model, as.integer(i),
            freq=double(rv$NumHaplo), hla=integer(rv$NumHaplo),
            haplo=character(rv$NumHaplo), snpidx = integer(rv$NumSNP),
            samp.num = integer(model$n.samp), acc = double(1),
            err=integer(1), NAOK=TRUE)
        if (rv$err != 0) stop(hlaErrMsg())

        res[[i]] <- list(
            samp.num = rv$samp.num,
            haplos = data.frame(freq = rv$freq, hla = model$hla.allele[rv$hla],
                haplo = rv$haplo, stringsAsFactors=FALSE),
            snpidx = rv$snpidx,
            outofbag.acc = rv$acc)
    }

    rv <- list(n.samp = model$n.samp, n.snp = model$n.snp,
        sample.id = model$sample.id, snp.id = model$snp.id,
        snp.position = model$snp.position, snp.allele = model$snp.allele,
        snp.allele.freq = model$snp.allele.freq,
        hla.locus = model$hla.locus,
        hla.allele = model$hla.allele, hla.freq = model$hla.freq,
        assembly = model$assembly,
        classifiers = res,
        appendix <- model$appendix)
    class(rv) <- "hlaAttrBagObj"
    rv
}


#######################################################################
# To combine two model objects of attribute bagging
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
        appendix = appendix)
    class(rv) <- "hlaAttrBagObj"
    rv
}


#######################################################################
# To get the top n individual classifiers
#

hlaSubModelObj <- function(obj, n)
{
    # check
    stopifnot(inherits(obj, "hlaAttrBagObj"))
    obj$classifiers <- obj$classifiers[1:n]
    obj
}


#######################################################################
# To get a "hlaAttrBagClass" class
#

hlaModelFromObj <- function(obj)
{
    # check
    stopifnot(inherits(obj, "hlaAttrBagObj"))

    # create an attribute bagging object
    rv <- .C(HIBAG_New,
        as.integer(obj$n.samp), as.integer(obj$n.snp), length(obj$hla.allele),
        model = integer(1), err=integer(1), NAOK=TRUE)
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
        rv <- .C(HIBAG_NewClassifierHaplo, ABmodel, length(tree$snpidx),
            as.integer(tree$snpidx-1), as.integer(snum), dim(tree$haplos)[1],
            as.double(tree$haplos$freq), as.integer(hla),
            as.character(tree$haplos$haplo), as.double(tree$outofbag.acc),
            err=integer(1), NAOK=TRUE)
        if (rv$err != 0) stop(hlaErrMsg())
    }

    # output
    rv <- list(n.samp = obj$n.samp, n.snp = obj$n.snp,
        sample.id = obj$sample.id, snp.id = obj$snp.id,
        snp.position = obj$snp.position, snp.allele = obj$snp.allele,
        snp.allele.freq = obj$snp.allele.freq,
        hla.locus = obj$hla.locus, hla.allele = obj$hla.allele,
        hla.freq = obj$hla.freq,
        assembly = as.character(obj$assembly)[1],
        model = ABmodel,
        appendix = obj$appendix)
    if (is.na(rv$assembly)) rv$assembly <- "unknown"

    class(rv) <- "hlaAttrBagClass"
    rv
}


#######################################################################
# summarize the "hlaAttrBagObj" object
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
        cat("\t# of individual classifiers: ", length(obj$classifiers),
            "\n", sep="")
        cat("\ttotal # of SNPs used: ", length(snpset), "\n", sep="")
        cat("\taverage # of SNPs in an individual classifier:",
            sprintf("%0.2f, sd: %0.2f, min: %d, max: %d\n",
            mean(numsnp), sd(numsnp), min(numsnp), max(numsnp)))
        cat("\taverage # of haplotypes in an individual classifier:",
            sprintf("%0.2f, sd: %0.2f, min: %d, max: %d\n",
            mean(numhaplo), sd(numhaplo), min(numhaplo), max(numhaplo)))
        cat("\taverage out-of-bag accuracy:",
            sprintf("%0.2f%%, sd: %0.2f%%, min: %0.2f%%, max: %0.2f%%\n",
            mean(outofbag.acc), sd(outofbag.acc), min(outofbag.acc),
            max(outofbag.acc)))

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
##########################################################################
#
# Linkage Disequilibrium
#

##########################################################################
# To calculate linkage disequilibrium between HLA locus and SNP markers
#

hlaGenoLD <- function(hla, geno)
{
    # check
    stopifnot(inherits(hla, "hlaAlleleClass"))
    if (inherits(geno, "hlaSNPGenoClass"))
    {
        stopifnot(dim(hla$value)[1] == length(geno$sample.id))
        if (any(hla$value$sample.id != geno$sample.id))
        {
            hla <- hlaAlleleSubset(hla, samp.sel =
                match(geno$sample.id, hla$value$sample.id))
        }
        geno <- geno$genotype
    } else if (is.matrix(geno))
    {
        stopifnot(is.numeric(geno))
        stopifnot(dim(hla$value)[1] == dim(geno)[2])
    } else if (is.vector(geno))
    {
        stopifnot(is.numeric(geno))
        stopifnot(dim(hla$value)[1] == length(geno))
        geno <- matrix(geno, ncol=1)
    } else {
        stop("geno should be `hlaSNPGenoClass', a vector or a matrix.")
    }

    # HLA alleles indicators
    alleles <- unique(c(hla$value$allele1, hla$value$allele2))
    alleles <- alleles[order(alleles)]
    allele.mat <- matrix(as.integer(0),
        nrow=length(hla$value$allele1), ncol=length(alleles))
    for (i in 1:length(alleles))
    {
        allele.mat[, i] <- (hla$value$allele1==alleles[i]) +
            (hla$value$allele2==alleles[i])
    }

    apply(geno, 1,
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
    invisible()
}



##########################################################################
# To visualize an attribute bagging model
#

plot.hlaAttrBagObj <- function(x, xlab=NULL, ylab=NULL,
    locus.color="red", locus.lty=2, locus.cex=1.25, assembly="auto", ...)
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
    pos.start <- info$pos.HLA.start[[x$hla.locus]]/1000
    pos.end <- info$pos.HLA.end[[x$hla.locus]]/1000

    # summary of the attribute bagging model
    desp <- summary(x, show=FALSE)

    # x - label, y - label
    if (is.null(xlab)) xlab <- "SNP Position (KB)"
    if (is.null(ylab)) ylab <- "Frequency of Use"

    # draw
    plot(x$snp.position/1000, desp$snp.hist, xlab=xlab, ylab=ylab, ...)
    abline(v=pos.start, col=locus.color, lty=locus.lty)
    abline(v=pos.end, col=locus.color, lty=locus.lty)
    text((pos.start + pos.end)/2, max(desp$snp.hist),
        paste("HLA", x$hla.locus, sep="-"),
        col=locus.color, cex=locus.cex)
}

print.hlaAttrBagObj <- function(x, ...)
{
    summary(x)
    invisible()
}



#######################################################################
# To get the error message
#

hlaErrMsg <- function()
{
    .Call(HIBAG_ErrMsg)
}



#######################################################################
# Internal R library functions
#######################################################################

.onAttach <- function(lib, pkg)
{
    # initialize HIBAG
    SSE.Flag <- .Call(HIBAG_Init)

    # information
    packageStartupMessage(
        "HIBAG (HLA Genotype Imputation with Attribute Bagging)")
    if (SSE.Flag == 1)
        s <- "Supported by Streaming SIMD Extensions (SSE2)"
    else if (SSE.Flag == 2)
        s <- "Supported by Streaming SIMD Extensions (SSE4.2 + hardware POPCNT)"
    else
        s <- ""
    if (s != "") packageStartupMessage(s)

    TRUE
}

.Last.lib <- function(libpath)
{
    # finalize HIBAG
    .Call(HIBAG_Done)
    TRUE
}
