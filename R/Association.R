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

##########################################################################
# Fit statistical models in assocation tests
#

hlaAssocTest <- function(hla, formula, data,
    model=c("dominant", "additive", "recessive", "genotype"),
    model.fit=c("glm"), showOR=FALSE, ...)
{
    stopifnot(inherits(hla, "hlaAlleleClass"))
    stopifnot(inherits(formula, "formula"))
    model <- match.arg(model)
    model.fit <- match.arg(model.fit)
    stopifnot(is.logical(showOR), length(showOR)==1L)
    if (missing(data)) 
        data <- environment(formula)

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

    # ans
    allele <- with(hla$value, hlaUniqueAllele(c(allele1, allele2)))

    # genotype distribution: dominant, additive, recessive, genotype
    mat <- NULL
    for (i in seq_along(allele))
    {
        s <- allele[i]
        switch(model,
            dominant = {
                    v <- with(hla$value, (allele1==s) | (allele2==s))
                    x <- c(sum(v==FALSE, na.rm=TRUE), sum(v==TRUE, na.rm=TRUE))
                },
            additive = {
                    v <- with(hla$value, (allele1==s) + (allele2==s))
                    x <- c(sum(v==0L, na.rm=TRUE), sum(v==1L, na.rm=TRUE),
                        sum(v==2L, na.rm=TRUE))
                },
            recessive = {
                    v <- with(hla$value, (allele1==s) & (allele2==s))
                    x <- c(sum(v==FALSE, na.rm=TRUE), sum(v==TRUE, na.rm=TRUE))
                },
            genotype = {
                    v <- with(hla$value, (allele1==s) + (allele2==s))
                    x <- c(sum(v==0L, na.rm=TRUE), sum(v==1L, na.rm=TRUE),
                        sum(v==2L, na.rm=TRUE))
                }
        )
        mat <- rbind(mat, x)
    }
    switch(model,
        dominant = { colnames(mat) <- c("[-/-]", "[-/h,h/h]") },
        additive = { colnames(mat) <- c("[-/-]", "[-/h]", "[h/h]") },
        recessive = { colnames(mat) <- c("[-/-,-/h]", "[h/h]") },
        genotype = { colnames(mat) <- c("[-/-]", "[-/h]", "[h/h]") }
    )
    ans <- as.data.frame(mat)
    rownames(ans) <- allele

    # chi-sq tests
    if (is.factor(y))
    {
        w1 <- w2 <- w3 <- rep(NaN, length(allele))
        for (i in seq_along(allele))
        {
            s <- allele[i]
            x <- switch(model,
                dominant = with(hla$value, (allele1==s) | (allele2==s)),
                additive = with(hla$value, (allele1==s) + (allele2==s)),
                recessive = with(hla$value, (allele1==s) & (allele2==s)),
                genotype = {
                    v <- with(hla$value, (allele1==s) + (allele2==s)) + 1L
                    attr(v, "levels") <- c("-/-", "-/h", "h/h")
                    attr(v, "class") <- "factor"
                    v
                }  
            )
            x <- as.factor(x)
            a <- try(v <- suppressWarnings(chisq.test(x, y)), silent=TRUE)
            if (!inherits(a, "try-error"))
            {
                w1[i] <- v$statistic
                w2[i] <- v$p.value
            }
            a <- try(v <- fisher.test(x, y), silent=TRUE)
            if (!inherits(a, "try-error"))
                w3[i] <- a$p.value
        }
        ans$chisq.st <- w1
        ans$chisq.p <- w2
        ans$fisher.p <- w3

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
        mat <- vector("list", length(allele))
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

            if (is.null(param$family))
            {
                a <- try({
                    if (is.factor(y))
                        m <- glm(formula, data=data, family=binomial, ...)
                    else
                        m <- glm(formula, data=data, ...)
                }, silent=TRUE)
            }
            if (!inherits(a, "try-error"))
            {
                z <- summary(m)$coefficients
                if (nrow(z) > 1L)
                {
                    ci <- confint.default(m)
                    v <- cbind(z[-1L,1L], ci[-1L,1L], ci[-1L,2L], z[-1L,4L])
                    v <- c(t(v))
                    nm <- rownames(z)[-1L]
                    names(v) <- c(rbind(paste0(nm, ".est"), paste0(nm, ".25%"),
                        paste0(nm, ".75%"), paste0(nm, ".pval")))
                    if (is.factor(y) & showOR)
                    {
                        if (model != "genotype")
                            nm <- c("h.est", "h.25%", "h.75%")
                        else
                            nm <- c("h1.est", "h1.25%", "h1.75%", "h2.est", "h2.25%", "h2.75%")
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
            ans <- cbind(ans, mat)
        } else
            warning(model.fit, " does not work.", immediate.=TRUE)
    }

    ans
}
