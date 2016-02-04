#######################################################################
#
# Package Name: HIBAG
# Description:
#   HIBAG -- HLA Genotype Imputation with Attribute Bagging
#
# HIBAG R package, HLA Genotype Imputation with Attribute Bagging
# Copyright (C) 2015   Xiuwen Zheng (zhengx@u.washington.edu)
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
# Association Tests
#

# Define a generic function
hlaAssocTest <- function(hla, ...)
{
    UseMethod("hlaAssocTest", hla)
}


# Print out the results
.assoc_show <- function(mat, pval.idx, show.all)
{
    v <- mat
    p <- as.matrix(v[, pval.idx])
    x <- sprintf("%.3f ", p); dim(x) <- dim(p)
    x[p < 0.001] <- "<0.001*"

    flag <- (p >= 0.001) & (p <= 0.05)
    flag[is.na(flag)] <- FALSE
    if (any(flag, na.rm=TRUE))
        x[flag] <- gsub(" ", "*", x[flag], fixed=TRUE)
    flag <- !is.finite(p)
    if (any(flag, na.rm=TRUE))
        x[flag] <- "."

    v[, pval.idx] <- x
    s <- format(v, digits=4L)
    s[is.na(v)] <- "."
    for (i in pval.idx)
        s[, i] <- as.character(s[, i])

    f <- apply(p, 1L, function(x) any(x <= 0.05, na.rm=TRUE))
    v <- NULL
    if (sum(f) > 0L)
    {
        v <- rbind(v, s[f, ])
        if ((sum(f) < nrow(s)) & (sum(!f) > 0L) & show.all)
            v <- rbind(v, "-----"=rep("", ncol(s)))
    }
    if ((sum(!f) > 0L) & show.all)
        v <- rbind(v, s[!f, ])
    print(v)

    invisible()
}



##########################################################################
# Fit statistical models in assocation tests for HLA alleles
#

# Association tests are applied to HLA classical alleles
hlaAssocTest.hlaAlleleClass <- function(hla, formula, data,
    model=c("dominant", "additive", "recessive", "genotype"),
    model.fit=c("glm"), prob.threshold=NaN, use.prob=FALSE, showOR=FALSE,
    verbose=TRUE, ...)
{
    stopifnot(inherits(hla, "hlaAlleleClass"))
    stopifnot(inherits(formula, "formula"))
    model <- match.arg(model)
    model.fit <- match.arg(model.fit)
    stopifnot(is.logical(use.prob), length(use.prob)==1L)
    stopifnot(is.logical(showOR), length(showOR)==1L)
    if (missing(data))
    {
        data <- environment(formula)
    } else if (is.data.frame(data))
    {
        if (nrow(data) != nrow(hla$value))
            stop("'hla' and 'data' must have the same length.")
    }

    stopifnot(is.numeric(prob.threshold), length(prob.threshold)==1L)
    if (is.finite(prob.threshold))
    {
        if (!is.data.frame(data))
            stop("'data' should be a data.frame, if 'prob.threshold' is used.")
        p <- hla$value$prob
        if (is.null(p))
            stop("No posterior probability in 'hla' ('hla$value$prob' should be available).")
        flag <- (p >= prob.threshold)
        flag[is.na(flag)] <- FALSE
        hla <- hlaAlleleSubset(hla, flag)
        data <- data[flag, ]
        if (verbose)
        {
            m <- sum(!flag)
            cat("Exclude ", m, " individual", .plural(m),
                " from the study due to the call threshold (",
                prob.threshold, ")\n", sep="")
        }
    }


    # formula
    fa <- format(formula)
    s <- unlist(strsplit(fa, "~", fixed=TRUE))
    if (length(s) != 2L)
        stop("Invalid `formula`: ", fa)
    if (s[1L] == "")
        stop("No dependent variable in `formula`: ", fa)

    yv <- trimws(s[1L])
    y <- data[[yv]]
    if (is.null(y))
        stop(sprintf("Dependent variable '%s' does not exist in `data`.", yv))
    else if (length(y) != nrow(hla$value))
        stop(sprintf("'hla' and '%s' must have the same length.", yv))

    if (isTRUE(use.prob))
    {
        if (is.null(hla$value$prob))
            stop("There is no posterior probability.")
    }


    # need proportion (%)
    flag <- is.factor(y)
    if (flag)
    {
        flag <- (nlevels(y) == 2L)
        if (flag) yy <- as.integer(y) - 1L
    } else {
        if (all(y %in% c(0,1), na.rm=TRUE))
        {
            y <- as.factor(y)
            flag <- (nlevels(y) == 2L)
            if (flag) yy <- as.integer(y) - 1L
        }
    }

    # ans
    allele <- with(hla$value, hlaUniqueAllele(c(allele1, allele2)))

    # genotype distribution: dominant, additive, recessive, genotype
    mat <- mat2 <- NULL
    for (i in seq_along(allele))
    {
        s <- allele[i]
        suppressWarnings(switch(model,
            dominant = {
                    v <- with(hla$value, (allele1==s) | (allele2==s))
                    x <- c(sum(v==FALSE, na.rm=TRUE), sum(v==TRUE, na.rm=TRUE))
                    if (flag)
                        b <- sapply(c(FALSE, TRUE), function(i) mean(yy[v==i], na.rm=TRUE))
                },
            additive = {
                    x <- with(hla$value, c(
                        sum(c(allele1, allele2) != s, na.rm=TRUE),
                        sum(c(allele1, allele2) == s, na.rm=TRUE)))
                    if (flag)
                    {
                        z <- with(hla$value, c(allele1, allele2) == s)
                        b <- c(mean(c(yy,yy)[!z], na.rm=TRUE),
                            mean(c(yy,yy)[z], na.rm=TRUE))
                    }
                },
            recessive = {
                    v <- with(hla$value, (allele1==s) & (allele2==s))
                    x <- c(sum(v==FALSE, na.rm=TRUE), sum(v==TRUE, na.rm=TRUE))
                    if (flag)
                        b <- sapply(c(FALSE, TRUE), function(i) mean(yy[v==i], na.rm=TRUE))
                },
            genotype = {
                    v <- with(hla$value, (allele1==s) + (allele2==s))
                    x <- c(sum(v==0L, na.rm=TRUE), sum(v==1L, na.rm=TRUE),
                        sum(v==2L, na.rm=TRUE))
                    if (flag)
                        b <- sapply(0:2, function(i) mean(yy[v==i], na.rm=TRUE))
                }
        ))
        if (flag)
        {
            b <- round(b * 100.0, 1)
            b[!is.finite(b)] <- NaN
            mat2 <- rbind(mat2, b)
        }
        mat <- rbind(mat, x)
    }

    switch(model,
        dominant  = { colnames(mat) <- c("[-/-]", "[-/h,h/h]") },
        additive  = { colnames(mat) <- c("[-]", "[h]") },
        recessive = { colnames(mat) <- c("[-/-,-/h]", "[h/h]") },
        genotype  = { colnames(mat) <- c("[-/-]", "[-/h]", "[h/h]") }
    )
    ans <- as.data.frame(mat)
    rownames(ans) <- allele
    pidx <- NULL

    if (!is.null(mat2))
    {
        colnames(mat2) <- paste0("%.", colnames(mat))
        ans <- cbind(ans, mat2)
    }

    # chi-sq tests
    if (is.factor(y))
    {
        w1 <- w2 <- w3 <- rep(NaN, length(allele))
        yy <- y
        if (model == "additive") yy <- rep(yy, 2L)
        for (i in seq_along(allele))
        {
            s <- allele[i]
            x <- switch(model,
                dominant = with(hla$value, (allele1==s) | (allele2==s)),
                additive = with(hla$value, c(allele1, allele2) == s),
                recessive = with(hla$value, (allele1==s) & (allele2==s)),
                genotype = {
                    v <- with(hla$value, (allele1==s) + (allele2==s)) + 1L
                    attr(v, "levels") <- c("-/-", "-/h", "h/h")
                    attr(v, "class") <- "factor"
                    v
                }  
            )
            x <- as.factor(x)
            a <- try(v <- suppressWarnings(chisq.test(x, yy)), silent=TRUE)
            if (!inherits(a, "try-error"))
            {
                w1[i] <- v$statistic
                w2[i] <- v$p.value
            }
            a <- try(v <- fisher.test(x, yy), silent=TRUE)
            if (!inherits(a, "try-error"))
                w3[i] <- a$p.value
        }
        ans$chisq.st <- w1
        ans$chisq.p <- w2
        ans$fisher.p <- w3
        pidx <- c(pidx, ncol(ans)-1L, ncol(ans))

    } else {
        if (model == "dominant")
        {
            mat <- matrix(NaN, nrow=length(allele), ncol=3L)
            colnames(mat) <- c("avg.[-/-]", "avg.[-/h,h/h]", "ttest.p")
            for (i in seq_along(allele))
            {
                s <- allele[i]
                x <- with(hla$value, (allele1==s) | (allele2==s))
                mat[i, 1L] <- suppressWarnings(mean(y[x==FALSE], na.rm=TRUE))
                mat[i, 2L] <- suppressWarnings(mean(y[x==TRUE], na.rm=TRUE))
                a <- try(v <- t.test(y[x==FALSE], y[x==TRUE]), silent=TRUE)
                if (!inherits(a, "try-error"))
                    mat[i, 3L] <- v$p.value
            }
        } else if (model == "recessive")
        {
            mat <- matrix(NaN, nrow=length(allele), ncol=3L)
            colnames(mat) <- c("avg.[-/-,-/h]", "avg.[h/h]", "ttest.p")
            for (i in seq_along(allele))
            {
                s <- allele[i]
                x <- with(hla$value, (allele1==s) & (allele2==s))
                mat[i, 1L] <- suppressWarnings(mean(y[x==FALSE], na.rm=TRUE))
                mat[i, 2L] <- suppressWarnings(mean(y[x==TRUE], na.rm=TRUE))
                a <- try(v <- t.test(y[x==FALSE], y[x==TRUE]), silent=TRUE)
                if (!inherits(a, "try-error"))
                    mat[i, 3L] <- v$p.value
            }
        } else {
            mat <- matrix(NaN, nrow=length(allele), ncol=4L)
            colnames(mat) <- c("avg.[-/-]", "avg.[-/h]", "avg.[h/h]", "anova.p")
            for (i in seq_along(allele))
            {
                s <- allele[i]
                x <- with(hla$value, (allele1==s) + (allele2==s))
                mat[i, 1L] <- suppressWarnings(mean(y[x==0L], na.rm=TRUE))
                mat[i, 2L] <- suppressWarnings(mean(y[x==1L], na.rm=TRUE))
                mat[i, 3L] <- suppressWarnings(mean(y[x==2L], na.rm=TRUE))
                x <- as.factor(x)
                a <- try(v <- aov(y ~ x), silent=TRUE)
                if (!inherits(a, "try-error"))
                {
                    v <- summary(v)
                    mat[i, 4L] <- v[[1L]]$`Pr(>F)`[1L]
                }
            }
        }

        ans <- cbind(ans, mat)
        pidx <- c(pidx, ncol(ans))
    }

    # regression
    vars <- attr(terms(formula), "term.labels")
    if (length(vars) > 0L)
    {
        if (!is.element("h", vars))
        {
            stop("Independent variable 'h' should be in `formula` ",
                "to include HLA genotypes, like '", fa, " + h'.")
        }
        if (!is.data.frame(data))
            stop("'data' should be `data.frame`.")

        param <- list(...)
        if (verbose)
        {
            if (is.null(param$family))
            {
                if (is.factor(y))
                    cat("Logistic regression")
                else
                    cat("Linear regression")
            } else
                cat("Regression [", format(param$family)[1L], "]", sep="")
            cat(" (", model, " model) with ", length(y), " individual",
                .plural(length(y)), ":\n", sep="")
        }

        mat <- vector("list", length(allele))
        summ <- NULL
        for (i in seq_along(allele))
        {
            s <- allele[i]
            data$h <- switch(model,
                dominant =
                    as.integer(with(hla$value, (allele1==s) | (allele2==s))),
                additive =
                    with(hla$value, (allele1==s) + (allele2==s)),
                recessive =
                    as.integer(with(hla$value, (allele1==s) & (allele2==s))),
                genotype =
                    as.factor(with(hla$value, (allele1==s) + (allele2==s)))
            )

            a <- try({
                if (is.null(param$family) & is.factor(y))
                {
                    if (!isTRUE(use.prob))
                    {
                        m <- glm(formula, data=data, family=binomial, ...)
                    } else {
                        prob <- hla$value$prob
                        m <- glm(formula, data=data, family=binomial, weights=prob, ...)
                    }
                } else {
                    if (!isTRUE(use.prob))
                    {
                        m <- glm(formula, data=data, ...)
                    } else {
                        prob <- hla$value$prob
                        m <- glm(formula, data=data, weights=prob, ...)
                    }
                }
                NULL
            }, silent=TRUE)
            if (!inherits(a, "try-error"))
            {
                summ <- summary(m)
                z <- summ$coefficients
                if (nrow(z) > 1L)
                {
                    ci <- confint.default(m)
                    v <- cbind(z[-1L,1L], ci[-1L,1L], ci[-1L,2L], z[-1L,4L])
                    v <- c(t(v))
                    nm <- rownames(z)[-1L]
                    names(v) <- c(rbind(paste0(nm, ".est"), paste0(nm, ".2.5%"),
                        paste0(nm, ".97.5%"), paste0(nm, ".pval")))
                    if (is.factor(y) & isTRUE(showOR))
                    {
                        if (model != "genotype")
                            nm <- c("h.est", "h.2.5%", "h.97.5%")
                        else
                            nm <- c("h1.est", "h1.2.5%", "h1.97.5%", "h2.est", "h2.2.5%", "h2.97.5%")
                        j <- match(nm, names(v))
                        nm <- names(v)
                        nm[j] <- paste0(nm[j], "_OR")
                        v[j] <- exp(v[j])
                        names(v) <- nm
                    }
                    mat[[i]] <- v
                }
            }
        }

        if (verbose & !is.null(summ$call))
        {
            s <- gsub("formula = formula", fa, format(summ$call), fixed=TRUE)
            cat("  ", s, "\n", sep="")
        }

        nm <- NULL
        for (i in seq_along(allele))
            nm <- c(nm, names(mat[[i]]))
        nm <- unique(nm)
        if (!is.null(nm))
        {
            for (i in seq_along(allele))
            {
                n <- length(mat[[i]])
                if (n == 0L)
                {
                    mat[[i]] <- rep(NA, length(nm))
                } else {
                    v <- mat[[i]]
                    mat[[i]] <- v[match(nm, names(v))]
                }
            }
            mat <- t(matrix(unlist(mat), nrow=length(nm)))
            colnames(mat) <- nm
            pidx <- c(pidx, ncol(ans) + seq(4L, ncol(mat), 4L))
            ans <- cbind(ans, mat)
        } else
            warning(model.fit, " does not work.", immediate.=TRUE)
    } else {
        if (verbose) cat(model, "model:\n")
        if (isTRUE(use.prob))
        {
            warning(ifelse(is.factor(y),
                "Chi-squared and Fisher's exact tests do not use posterior probabilities.",
                "T test or ANOVA does not use posterior probabilities."),
                immediate.=TRUE)
        }
    }

    if (verbose)
        .assoc_show(ans, pidx, TRUE)

    invisible(ans)
}



# Association tests are applied to HLA protein sequences
hlaAssocTest.hlaAASeqClass <- function(hla, formula, data,
    model=c("dominant", "additive", "recessive", "genotype"),
    model.fit=c("glm"), prob.threshold=NaN, use.prob=FALSE, showOR=FALSE,
    show.all=FALSE, verbose=TRUE, ...)
{
    stopifnot(inherits(hla, "hlaAASeqClass"))
    stopifnot(inherits(formula, "formula"))
    model <- match.arg(model)
    model.fit <- match.arg(model.fit)
    stopifnot(is.logical(use.prob), length(use.prob)==1L)
    stopifnot(is.logical(showOR), length(showOR)==1L)
    stopifnot(is.logical(show.all), length(show.all)==1L)
    if (missing(data))
    {
        data <- environment(formula)
    } else if (is.data.frame(data))
    {
        if (nrow(data) != nrow(hla$value))
            stop("'hla' and 'data' must have the same length.")
    }

    stopifnot(is.numeric(prob.threshold), length(prob.threshold)==1L)
    if (is.finite(prob.threshold))
    {
        if (!is.data.frame(data))
            stop("'data' should be a data.frame, if 'prob.threshold' is used.")
        p <- hla$value$prob
        if (is.null(p))
            stop("No posterior probability in 'hla' ('hla$value$prob' should be available).")
        flag <- (p >= prob.threshold)
        flag[is.na(flag)] <- FALSE
        hla <- hlaAlleleSubset(hla, flag)
        data <- data[flag, ]
        if (verbose)
        {
            m <- sum(!flag)
            cat("Exclude ", m, " individual", .plural(m),
                " from the study due to the call threshold (",
                prob.threshold, ")\n", sep="")
        }
    }


    # formula
    fa <- format(formula)
    s <- unlist(strsplit(fa, "~", fixed=TRUE))
    if (length(s) != 2L)
        stop("Invalid `formula`: ", fa)
    if (s[1L] == "")
        stop("No dependent variable in `formula`: ", fa)

    yv <- trimws(s[1L])
    y <- data[[yv]]
    if (is.null(y))
        stop(sprintf("Dependent variable '%s' does not exist in `data`.", yv))
    else if (length(y) != nrow(hla$value))
        stop(sprintf("'hla' and '%s' must have the same length.", yv))

    if (isTRUE(use.prob))
    {
        if (is.null(hla$value$prob))
            stop("There is no posterior probability.")
    }


    # need proportion (%)
    flag <- is.factor(y)
    if (flag)
    {
        flag <- (nlevels(y) == 2L)
        if (flag) yy <- as.integer(y) - 1L
    } else {
        if (all(y %in% c(0,1), na.rm=TRUE))
        {
            y <- as.factor(y)
            flag <- (nlevels(y) == 2L)
            if (flag) yy <- as.integer(y) - 1L
        }
    }

    # regression
    vars <- attr(terms(formula), "term.labels")
    if (length(vars) > 0L)
    {
        if (!is.element("h", vars))
        {
            stop("Independent variable 'h' should be in `formula` ",
                "to include HLA genotypes, like '", fa, " + h'.")
        }
        if (!is.data.frame(data))
            stop("'data' should be `data.frame`.")
    }


    if (is.factor(y))
    {
        if (length(vars) > 0L)
        {
            param <- list(...)
            if (verbose)
            {
                if (is.null(param$family))
                    cat("Logistic regression")
                else
                    cat("Regression [", format(param$family)[1L], "]", sep="")
                cat(" (", model, " model) with ", length(y), " individual",
                    .plural(length(y)), ":\n", sep="")
            }
        }

        y2 <- rep(y, 2L)
        matseq <- .matrix_sequence(c(hla$value$allele1, hla$value$allele2))
        pos <- 1L - hla$start.position + 1L

        z <- apply(matseq, 1L, FUN=function(x)
        {
            x[x == 42] <- NA  # 42 = '*'
            xx <- as.factor(x)
            xl <- as.integer(levels(xx))
            s <- rawToChar(as.raw(xl))
            a <- try(v <- fisher.test(xx, y2), silent=TRUE)
            pos <<- pos + 1L
            i <- pos + hla$start.position - 2L

            rv <- data.frame(
                pos = pos - 1L,
                num = sum(!is.na(x)),
                ref = substr(hla$reference, i, i),
                poly = paste(unlist(strsplit(s, "", fixed=TRUE)), collapse=","),
                fisher.p = ifelse(inherits(a, "try-error"), NaN, a$p.value),
                stringsAsFactors=FALSE)

            if (length(vars) > 0L)
            {
                if (length(xl) == 2L) xl <- xl[1L]
                a1 <- x[seq.int(1L, length(x)/2)]
                a2 <- x[seq.int(length(x)/2 + 1L, length(x))]

                tv <- NULL
                for (k in xl)
                {
                    if (k == 45L)
                    {   # -, reference, - vs. others
                        data$h <- switch(model,
                            dominant  = as.integer((a1!=k) | (a2!=k)),
                            additive  = (a1!=k) + (a2!=k),
                            recessive = as.integer((a1!=k) & (a2!=k)),
                            genotype = as.factor((a1!=k) + (a2!=k))
                        )
                    } else {
                        data$h <- switch(model,
                            dominant  = as.integer((a1==k) | (a2==k)),
                            additive  = (a1==k) + (a2==k),
                            recessive = as.integer((a1==k) & (a2==k)),
                            genotype = as.factor((a1==k) + (a2==k))
                        )
                    }

                    a <- try({
                        if (!isTRUE(use.prob))
                        {
                            m <- glm(formula, data=data, family=binomial, ...)
                        } else {
                            prob <- hla$value$prob
                            m <- glm(formula, data=data, family=binomial,
                                weights=prob, ...)
                        }
                        NULL
                    }, silent=TRUE)

                    if (!inherits(a, "try-error"))
                    {
                        summ <- summary(m)
                        z <- summ$coefficients
                        if (nrow(z) > 1L)
                        {
                            ci <- confint.default(m)
                            v <- cbind(z[-1L,1L], ci[-1L,1L], ci[-1L,2L],
                                z[-1L,4L])
                            v <- c(t(v))
                            nm <- rownames(z)[-1L]
                            names(v) <- c(rbind(paste0(nm, ".est"),
                                paste0(nm, ".2.5%"), paste0(nm, ".97.5%"),
                                paste0(nm, ".pval")))
                            if (isTRUE(showOR))
                            {
                                if (model != "genotype")
                                    nm <- c("h.est", "h.2.5%", "h.97.5%")
                                else
                                    nm <- c("h1.est", "h1.2.5%", "h1.97.5%",
                                        "h2.est", "h2.2.5%", "h2.97.5%")
                                j <- match(nm, names(v))
                                nm <- names(v)
                                nm[j] <- paste0(nm[j], "_OR")
                                v[j] <- exp(v[j])
                                names(v) <- nm
                            }

                            tv <- rbind(tv, v)
                        }
                    }
                }
                if (!is.null(tv))
                {
                    if (!is.data.frame(tv))
                    {
                        rownames(tv) <- NULL
                        tv <- as.data.frame(tv)
                    }
                    xl <- sapply(as.integer(levels(xx)), function(x)
                        rawToChar(as.raw(x)))
                    tv <- cbind(
                        amino.acid = sapply(seq_along(xl), function(i)
                            if (xl[i] == "-")
                                paste(xl[i], "vs", paste(xl[-i], collapse=","))
                            else
                                paste(paste(xl[-i], collapse=","), "vs", xl[i])
                            ),
                        tv, stringsAsFactors=FALSE)
                }

                n <- 0L
                if (!is.null(tv)) n <- nrow(tv)
                if (n > 0L)
                {
                    if (n > 1L)
                    {
                        rv <- as.data.frame(sapply(rv, function(x) rep(x, n),
                            simplify=FALSE), stringsAsFactors=FALSE)
                    }
                    rv <- cbind(rv, tv)
                }
            }

            rv
        })

        ans <- data.frame(
            pos = unlist(sapply(z, function(x) x$pos)),
            num = unlist(sapply(z, function(x) x$num)),
            ref = unlist(sapply(z, function(x) x$ref)),
            poly = unlist(sapply(z, function(x) x$poly)),
            fisher.p = unlist(sapply(z, function(x) x$fisher.p)),
            stringsAsFactors=FALSE)

        n <- max(lengths(z))
        pidx <- c(5L)
        if (n > 5L)
        {
            for (i in 6L:n)
            {
                ans <- cbind(ans, unlist(sapply(z, function(x)
                    if (i <= ncol(x)) x[,i] else NA
                )))
            }
            names(ans) <- names(z[[match(n, lengths(z))]])
            pidx <- c(pidx, seq.int(7L, n, 4L) + 3L)
        }

        if (verbose)
            .assoc_show(ans, pidx, show.all)

    } else {
        stop(sprintf("Dependent variable '%s' should be factor.", yv))
    }



    invisible(ans)
}
