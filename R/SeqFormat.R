#######################################################################
#
# Package Name: HIBAG
# Description:
#   HIBAG -- HLA Genotype Imputation with Attribute Bagging
#
# HIBAG R package, HLA Genotype Imputation with Attribute Bagging
# Copyright (C) 2011-2017   Xiuwen Zheng (zhengx@u.washington.edu)
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

# Package-wide variable
.packageEnv <- new.env()


# Get P Codes
.pcode <- function(release)
{
    varnm <- paste0(release, ".pcode")
    if (!exists(varnm, envir=.packageEnv))
    {
        fn <- system.file("extdata", release, "hla_nom_p.txt.xz",
            package="HIBAG")
        s <- readLines(fn)
        s <- s[substr(s, 1L, 1L) != "#"]  # remove comments

        z <- strsplit(s, ";", fixed=TRUE)
        a1 <- sapply(z, `[`, i=1L)
        a2 <- sapply(z, `[`, i=2L)
        a3 <- sapply(z, `[`, i=3L)
        a3[is.na(a3)] <- a2[is.na(a3)]
        a1 <- paste0(a1, a3)

        pcode <- data.frame(code=a1, allele=a2, stringsAsFactors=FALSE)
        assign(varnm, pcode, envir=.packageEnv)
    } else {
        pcode <- get(varnm, envir=.packageEnv)
    }
    pcode
}


# Get G Codes
.gcode <- function(release)
{
    varnm <- paste0(release, ".gcode")
    if (!exists(varnm, envir=.packageEnv))
    {
        fn <- system.file("extdata", release, "hla_nom_g.txt.xz",
            package="HIBAG")
        s <- readLines(fn)
        s <- s[substr(s, 1L, 1L) != "#"]  # remove comments

        z <- strsplit(s, ";", fixed=TRUE)
        a1 <- sapply(z, `[`, i=1L)
        a2 <- sapply(z, `[`, i=2L)
        a3 <- sapply(z, `[`, i=3L)
        a3[is.na(a3)] <- a2[is.na(a3)]
        a1 <- paste0(a1, a3)

        gcode <- data.frame(code=a1, allele=a2, stringsAsFactors=FALSE)
        assign(varnm, gcode, envir=.packageEnv)
    } else {
        gcode <- get(varnm, envir=.packageEnv)
    }
    gcode
}


# Get feature information
.feature <- function(release)
{
    varnm <- paste0(release, ".feature.info")
    if (!exists(varnm, envir=.packageEnv))
    {
        fn <- system.file("extdata", release, "FeatureInfo.txt",
            package="HIBAG")
        feature <- read.table(fn, header=TRUE, sep="\t", quote="",
            stringsAsFactors=FALSE)
        assign(varnm, feature, envir=.packageEnv)
    } else {
        feature <- get(varnm, envir=.packageEnv)
    }
    feature
}


# Get Protein Sequences
.protein <- function(hla.id, release)
{
    varnm <- paste0(release, ".", hla.id, "_prot")
    if (!exists(varnm, envir=.packageEnv))
    {
        fn <- system.file("extdata", release, "SeqAlign",
            sprintf("%s_prot.txt.xz", tolower(hla.id)), package="HIBAG")

        s <- readLines(fn)
        s1 <- trimws(s[7L], which="right")
        s2 <- trimws(s[8L], which="right")
        stopifnot(substr(s1, nchar(s1), nchar(s1)) == "1")

        # identify the start position
        x <- scan(text=s[9L], what=character(), n=1, quiet=TRUE)
        ss <- gsub(x, paste(rep(" ", nchar(x)), collapse=""), s[9L], fixed=TRUE)
        ss <- substr(ss, 1, nchar(s2))
        ss <- gsub(" ", "", ss)
        start <- nchar(ss)

        # filter
        shead <- sprintf(" %s*", hla.id)
        flag <- (substr(s, 1L, nchar(shead)) == shead)
        s <- s[flag]
        s <- gsub(shead, "", s, fixed=TRUE)

        # get a matrix, mat[,1] -- allele, mat[,2] -- sequence
        mat <- t(sapply(s, function(x) {
            v <- scan(text=x, what=character(), quiet=TRUE)
            c(v[1L], paste(v[-1L], collapse=""))
        }, USE.NAMES=FALSE))

        # merge
        id <- mat[, 1L]
        ss <- t(sapply(unique(id), function(x) {
            c(x, paste(mat[id == x, 2L], collapse=""))
        }, USE.NAMES=FALSE))
        reference <- ss[1L,2L]
        ss[1L,2L] <- paste(rep("-", nchar(ss[1L,2L])), collapse="")

        # remove dots
        if (hla.id != "DQB1")
        {
            # DQB1 reference has deletions
            .Call(HIBAG_SeqRmDot, reference, ss)
        }

        # get feature information
        fea <- .feature(release)
        fea <- fea[fea$id==hla.id, ]
        fea <- fea[fea$name %in% paste("Exon", 1:16), ]
        len <- fea$end - fea$start + 1L
        v <- cumsum(len)
        end <- (v %/% 3L) + (v %% 3L)
        v <- c(1L, v[-length(v)] + 1L)
        st <- (v + 2L) %/% 3L

        prot <- list(reference=reference, start=start, allele=ss[,1L],
            sequence=ss[,2L],
            feature = data.frame(id=fea$name, start=st, end=end,
                stringsAsFactors=FALSE))
        assign(varnm, prot, envir=.packageEnv)
    } else {
        prot <- get(varnm, envir=.packageEnv)
    }
    prot
}


# region
.subset <- function(locus, region, seq)
{
    if (region %in% c("P.code", "G.code"))
    {
        if (locus %in% c("A", "B", "C"))
        {
            # classical I
            start <- seq$feature$start[2L]
            end <- seq$feature$end[3L]
        } else {
            # classical II
            start <- seq$feature$start[2L]
            end <- seq$feature$end[2L]
        }
        c(start, end)
    } else
        NULL
}



##########################################################################
# Convert HLA Alleles to Amino Acid Sequences
#

hlaConvSequence <- function(hla=character(), locus=NULL,
    method=c("protein", "protein_reference"),
    code=c("exact", "P.code", "G.code", "P.code.merge", "G.code.merge"),
    region=c("auto", "all", "P.code", "G.code"),
    release=c("v3.22.0"), replace=NULL)
{
    stopifnot(is.character(hla) | inherits(hla, "hlaAlleleClass"))
    method <- match.arg(method)
    code <- match.arg(code)
    region <- match.arg(region)
    release <- match.arg(release)
    stopifnot(is.null(replace) | is.character(replace))
    if (is.character(replace))
    {
        stopifnot(is.vector(replace))
        if (is.null(names(replace)))
        {
            stop("'replace' should be a character vector with names, ",
                "like c(\"09:02\"=\"107:01\")")
        }
    }

    if (region == "auto")
    {
        if (code == "exact")
            region <- "all"
        else if (code %in% c("P.code", "P.code.merge"))
            region <- "P.code"
        else if (code %in% c("G.code", "G.code.merge"))
            region <- "G.code"
    }

    hlalist <- c("A", "B", "C", "DRB1", "DQA1", "DQB1", "DPB1", "DPA1")

    if (inherits(hla, "hlaAlleleClass"))
    {
        if (!is.null(locus))
            warning("'locus' should be NULL.")
        if (code %in% c("P.code", "G.code"))
            stop("'code' should be 'exact', 'P.code.merge' or 'G.code.merge'.")
        locus <- hla$locus

        if (!(locus %in% hlalist))
            stop("`hla$locus` should be one of ", paste(hlalist, collapse=", "))

        p <- .protein(locus, release)
        s <- hlaConvSequence(c(hla$value$allele1, hla$value$allele2),
            locus=hla$locus, method=method, code=code, region=region,
            release=release, replace=replace)

        v <- .subset(locus, region, p)
        if (is.null(v)) v <- c(1L, 1000000L)

        n <- length(s)
        ans <- list(locus = hla$locus,
            pos.start = hla$locus.pos.start, pos.end = hla$locus.pos.end,
            value = data.frame(sample.id = hla$value$sample.id,
                allele1 = s[seq.int(1L, n/2)], allele2 = s[seq.int(n/2+1L, n)],
                stringsAsFactors=FALSE),
            assembly = hla$assembly,
            start.position = p$start - v[1L] + 1L,
            reference = substring(p$reference, v[1L], v[2L])
        )
        if (!is.null(hla$value$prob))
            ans$value$prob <- hla$value$prob
        class(ans) <- "hlaAASeqClass"

    } else {

        stopifnot(length(locus) == 1L)
        stopifnot(is.vector(hla))
        if (!(locus %in% hlalist))
            stop("`locus` should be one of ", paste(hlalist, collapse=", "))

        if (method == "protein_reference")
        {
            if (length(hla) > 0L)
                warning("the argument 'hla' is ignored.")
            p <- .protein(locus, release)
            p$feature <- cbind(p$feature, sequence=apply(p$feature[, -1L], 1L,
                function(x) substr(p$reference, x[1L], x[2L])))
            ans <- list(reference = p$reference, start.position = p$start,
                feature=p$feature)

        } else if (method == "protein")
        {
            # replace
            if (is.character(replace))
            {
                i <- match(hla, names(replace))
                hla[!is.na(i)] <- replace[i[!is.na(i)]]
            }

            unihla <- unique(hla)
            unihla <- unihla[!is.na(unihla)]

            seq <- .protein(locus, release)
            ss <- seq$sequence[match(unihla, seq$allele)]
            ans <- sapply(seq_along(ss), FUN=function(i) {
                if (is.na(ss[i])) NULL else ss[i]
            }, simplify=FALSE, USE.NAMES=TRUE)
            names(ans) <- unihla

            if (code %in% c("P.code", "P.code.merge"))
            {
                pcode <- .pcode(release)
                if (anyNA(ss))
                {
                    h1 <- paste0(locus, "*", unihla[is.na(ss)])
                    h2 <- paste0(h1, "P")
                    i1 <- match(h1, pcode$code)
                    i2 <- match(h2, pcode$code)
                    i1[is.na(i1)] <- i2[is.na(i1)]
                    mat <- sapply(i1, FUN=function(i) {
                        if (!is.na(i))
                        {
                            a <- unlist(strsplit(pcode$allele[i], "/", fixed=TRUE))
                            m <- seq$sequence[match(a, seq$allele)]
                            names(m) <- a
                            m
                        } else
                            NULL
                    }, simplify=FALSE, USE.NAMES=FALSE)
                    ans[is.na(ss)] <- mat
                }
            } else if (code %in% c("G.code", "G.code.merge"))
            {
                gcode <- .gcode(release)
                if (anyNA(ss))
                {
                    h1 <- paste0(locus, "*", unihla[is.na(ss)])
                    h2 <- paste0(h1, "G")
                    i1 <- match(h1, gcode$code)
                    i2 <- match(h2, gcode$code)
                    i1[is.na(i1)] <- i2[is.na(i1)]
                    mat <- sapply(i1, FUN=function(i) {
                        if (!is.na(i))
                        {
                            a <- unlist(strsplit(gcode$allele[i], "/", fixed=TRUE))
                            m <- seq$sequence[match(a, seq$allele)]
                            names(m) <- a
                            m
                        } else
                            NULL
                    }, simplify=FALSE, USE.NAMES=FALSE)
                    ans[is.na(ss)] <- mat
                }
            }

            # any warning?
            len <- lengths(ans)
            if (sum(len > 1L) > 0L)
            {
                message("Allelic ambiguity: ",
                    paste(unihla[len > 1L], collapse=", "))
            }
            if (sum(len == 0L) > 0L)
            {
                warning("No matching: ",
                    paste(unihla[len == 0L], collapse=", "),
                    "\n  See: http://hla.alleles.org/alleles/text_index.html",
                    immediate.=TRUE)
            }

            # region
            v <- .subset(locus, region, seq)
            if (!is.null(v))
            {
                ans <- sapply(ans, function(x) {
                    if (!is.null(x)) substr(x, v[1L], v[2L]) else NULL
                }, simplify=FALSE, USE.NAMES=TRUE)
            }

            # merge ?
            if (code %in% c("exact", "P.code.merge", "G.code.merge"))
            {
                ans <- unname(sapply(ans, FUN=function(x)
                    .Call(HIBAG_SeqMerge, x),
                    simplify=TRUE, USE.NAMES=FALSE))
            }

            # output
            ans <- ans[match(hla, unihla)]

        } else
            ans <- NULL
    }

    ans
}



#######################################################################
# Summary a "hlaAASeqClass" object
#

.matrix_sequence <- function(aa)
{
    lt <- sapply(aa, function(s) as.integer(charToRaw(s)),
        simplify=FALSE, USE.NAMES=FALSE)
    n <- max(lengths(lt))
    ix <- seq_len(n)
    lt <- sapply(lt, function(a) a[ix])
    matrix(unlist(lt), nrow=n)
}

summary.hlaAASeqClass <- function(object, poly.only=TRUE, head=0L,
    verbose=TRUE, ...)
{
    # check
    stopifnot(inherits(object, "hlaAASeqClass"))
    stopifnot(is.logical(poly.only), length(poly.only)==1L)
    stopifnot(is.numeric(head), length(head)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    m <- .matrix_sequence(c(object$value$allele1, object$value$allele2))
    level <- unique(c(m))
    level <- level[!is.na(level)]
    level <- level[order(level)]
    levelstr <- sapply(level, function(x) rawToChar(as.raw(x)))

    mt <- apply(m, 1, function(x)
        c(sum(is.finite(x), na.rm=TRUE),
        sapply(level, function(y) sum(x==y, na.rm=TRUE))))
    mt <- t(matrix(unlist(mt), nrow=length(level)+1L))
    colnames(mt) <- c("Num", levelstr)
    mt <- cbind(Pos=seq_len(nrow(mt))-object$start.position+1L, mt)

    if (isTRUE(poly.only))
    {
        i1 <- match("Num", colnames(mt))
        i2 <- match("-", colnames(mt))
        if (!is.na(i1) & !is.na(i2))
            mt <- mt[mt[,i1] != mt[,i2], ]
    }

    if (verbose)
    {
        z <- mt
        storage.mode(z) <- "character"
        z[z == "0"] <- "."
        z <- rbind(c("Pos", "Num", levelstr), z)
        z <- format(z, justify="right")
        if (head < 1L) head <- .Machine$integer.max - 1L
        head <- head + 1L
        for (i in 1L:min(head, nrow(z)))
            cat(z[i, ], "\n")
        if (nrow(z) > head)
            cat("......\n")
    }

    # return
    invisible(mt)
}
