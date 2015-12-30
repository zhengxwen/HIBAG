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

        mat <- t(sapply(s, function(x) {
            v <- scan(text=x, what=character(), quiet=TRUE)
            c(v[1L], paste(v[-1L], collapse=""))
        }, USE.NAMES=FALSE))

        # combine
        i <- which(mat[1L,1L] == mat[,1L])
        len <- i[2L] - i[1L]
        ss <- mat[seq.int(1L, len), ]
        for (j in seq_along(i[-1L]))
        {
            k <- seq.int(len*j + 1L, length.out=len)
            stopifnot(identical(ss[,1L], mat[k, 1L]))
            ss[,2L] <- paste0(ss[,2L], mat[k, 2L])
        }
        reference <- ss[1L,2L]
        ss[1L,2L] <- paste(rep("-", nchar(ss[1L,2L])), collapse="")

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




##########################################################################
# Convert HLA Alleles to Amino Acid Sequences
#

hlaConvSequence <- function(hla=character(), locus=NULL,
    method=c("protein", "protein_reference"),
    code=c("exact", "P.code", "G.code", "P.code.merge", "G.code.merge"),
    region=c("all", "P.code", "G.code"),
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

        n <- length(s)
        ans <- list(locus = hla$locus,
            pos.start = hla$locus.pos.start, pos.end = hla$locus.pos.end,
            value = data.frame(sample.id = hla$value$sample.id,
                allele1 = s[seq.int(1L, n/2)], allele2 = s[seq.int(n/2+1L, n)],
                stringsAsFactors=FALSE),
            assembly = hla$assembly,
            start.position = ifelse(region=="all", p$start, 1L)
        )
        if (!is.null(hla$value$prob))
            ans$value$prob <- hla$value$prob
        class(ans) <- "hlaSeqClass"

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
            n1 <- sum(len == 0L)
            n2 <- sum(len > 1L)
            if ((n1 > 0L) & (n2 > 0L))
            {
                warning("\nAllelic ambiguity: ",
                    paste(unihla[len > 1L], collapse=", "),
                    "\nNo matching: ",
                    paste(unihla[len == 0L], collapse=", "),
                    call.=FALSE, immediate.=TRUE)
            } else {
                if (n2 > 0L)
                {
                    warning("\nAllelic ambiguity: ",
                        paste(unihla[len > 1L], collapse=", "),
                        call.=FALSE, immediate.=TRUE)
                }
                if (n1 > 0L)
                {
                    warning("\nNo matching: ",
                        paste(unihla[len == 0L], collapse=", "),
                        call.=FALSE, immediate.=TRUE)
                }
            }

            # region
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

                ans <- sapply(ans, function(x) {
                    if (!is.null(x)) substr(x, start, end) else NULL
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
