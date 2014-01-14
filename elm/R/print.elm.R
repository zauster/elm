##
## non-verbose:
##
##        Exact linear models
##
## data: Y and X
##
## Tested coefficients:
##        H_0     Estimate     Threshold    Test statistic
## weight > 0         1.23          0.23              0.08
##                             Rejection       chosen Test
##                                    No         Bernoulli
##
##        H_0     Estimate     Threshold    Test statistic
## gender > 0         3.89          12.5              18.4
##                             Rejection       chosen Test
##                                   Yes         Hoeffding
##
## parameter:
##    Y in [0, 10]
##    alpha: 0.05
##    iterations: 1000

## verbose:
##        Exact linear models
##
## data: Y and X
##
## Tested coefficients:
##         H_0     Estimate     Threshold    Test statistic
## weight <= 0         1.23          0.23              0.08
##                              Rejection       chosen Test
##                                     No         Bernoulli
##       betaj         Test     Threshold            TypeII
##        1.42    Bernoulli          0.23              0.50
##        1.42    Hoeffding          3.44              0.54
##        1.42     Cantelli         13.93              0.56
##        ...
##
##
## parameter:
##    Y in [0, 10]
##    alpha: 0.05
##    iterations: 1000


print.elm <- function(x, ..., digits = max(3L, getOption("digits") - 3L))
    {
        cat("\n")
        cat(strwrap(x$method, prefix = "\t"), sep="\n")
        cat("\n")
        cat("data:", x$yname, "and", x$xname, "\n")
        cat("n = ", x$parameter$n, ", m = ", x$parameter$m, "\n", sep = "")

        cat("Tested coefficients:", length(x$coefTests), "\n")
        for(lst in x$coefTests)
            {
                for(i in 2:4)
                    {
                        lst$chosenTest[[i]] <- format(round(lst$chosenTest[[i]],
                                                            digits = digits),
                                              digits = digits)
                    }
                cat("\n")
                cat("----------- Coefficient:", lst$coefname,
                    paste(rep("-", times = 30 - nchar(lst$coefname)),
                          collapse = ""), "\n")
                print(unlist(lst$chosenTest), quote = FALSE, right = TRUE)
                ## print(unlist(lst$chosenTest[1:3]), quote = FALSE, right = TRUE)
                ## print(unlist(lst$chosenTest[4:6]), quote = FALSE, right = TRUE)
                cat("TypeII minimizing beta in the alternative: ", lst$betaj, sep = "")
                cat("\n")

                if(x$verbose == TRUE)
                    {
                        cat("\nOLS Nonstandardized:\n")
                        print(format(lst$OLSNonstandardized$NonstandardizedTests,
                                     digits = digits), quote = FALSE, right = TRUE)
                        print(format(lst$OLSNonstandardized$NonstandardizedTypeII,
                                     digits = digits), quote = FALSE, right = TRUE)

                        cat("\nMM Nonstandardized:\n")
                        print(format(lst$MMNonstandardized$NonstandardizedTests,
                                     digits = digits), quote = FALSE, right = TRUE)
                        print(format(lst$MMNonstandardized$NonstandardizedTypeII,
                                     digits = digits), quote = FALSE, right = TRUE)

                        cat("\nOLS Bernoulli:\n")
                        print(format(lst$OLSBernoulli$BernoulliTest,
                                     digits = digits), quote = FALSE, right = TRUE)

                        cat("\nMM Bernoulli:\n")
                        print(format(lst$MMBernoulli$BernoulliTest,
                                     digits = digits), quote = FALSE, right = TRUE)
                    }
            }



        cat("\nparameters:\n")
        if(!is.null(x$parameter$bounds))
            {
                cat(paste("   ", x$yname, " in ", x$parameter$bounds, "\n", sep = ""))
            }
        cat("   alpha:", x$parameter$alpha)
        if(!is.null(x$parameter$iterations))
            {
                cat("\n   iterations:", x$parameter$iterations)
            }
        cat("\n")

    }

## printCoefs <- function(x) ## x being a list of results from a test of
##     ## a coefficient
##     {

##     }


## print.summary.lm <- function (x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor,
##     signif.stars = getOption("show.signif.stars"), ...)
## {
##     cat("\nCallxxx:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
##         "\n\n", sep = "")
##     resid <- x$residuals
##     df <- x$df
##     rdf <- df[2L]
##     cat(if (!is.null(x$weights) && diff(range(x$weights)))
##         "Weighted ", "Residuals:\n", sep = "")
##     if (rdf > 5L) {
##         nam <- c("Min", "1Q", "Median", "3Q", "Max")
##         rq <- if (length(dim(resid)) == 2L)
##             structure(apply(t(resid), 1L, quantile), dimnames = list(nam,
##                 dimnames(resid)[[2L]]))
##         else {
##             zz <- zapsmall(quantile(resid), digits + 1L)
##             structure(zz, names = nam)
##         }
##         print(rq, digits = digits, ...)
##     }
##     else if (rdf > 0L) {
##         print(resid, digits = digits, ...)
##     }
##     else {
##         cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
##         cat("\n")
##     }
##     if (length(x$aliased) == 0L) {
##         cat("\nNo Coefficients\n")
##     }
##     else {
##         if (nsingular <- df[3L] - df[1L])
##             cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
##                 sep = "")
##         else cat("\nCoefficients:\n")
##         coefs <- x$coefficients
##         if (!is.null(aliased <- x$aliased) && any(aliased)) {
##             cn <- names(aliased)
##             coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn,
##                 colnames(coefs)))
##             coefs[!aliased, ] <- x$coefficients
##         }
##         print(str(coefs))
##         cat("\n")
##         printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
##             na.print = "NA", ...)
##     }
##     cat("\nResidual standard error:", format(signif(x$sigma,
##         digits)), "on", rdf, "degrees of freedom")
##     cat("\n")
##     if (nzchar(mess <- naprint(x$na.action)))
##         cat("  (", mess, ")\n", sep = "")
##     if (!is.null(x$fstatistic)) {
##         cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits))
##         cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared,
##             digits = digits), "\nF-statistic:", formatC(x$fstatistic[1L],
##             digits = digits), "on", x$fstatistic[2L], "and",
##             x$fstatistic[3L], "DF,  p-value:", format.pval(pf(x$fstatistic[1L],
##                 x$fstatistic[2L], x$fstatistic[3L], lower.tail = FALSE),
##                 digits = digits))
##         cat("\n")
##     }
##     correl <- x$correlation
##     if (!is.null(correl)) {
##         p <- NCOL(correl)
##         if (p > 1L) {
##             cat("\nCorrelation of Coefficients:\n")
##             if (is.logical(symbolic.cor) && symbolic.cor) {
##                 print(symnum(correl, abbr.colnames = NULL))
##             }
##             else {
##                 correl <- format(round(correl, 2), nsmall = 2,
##                   digits = digits)
##                 correl[!lower.tri(correl)] <- ""
##                 print(correl[-1, -p, drop = FALSE], quote = FALSE)
##             }
##         }
##     }
##     cat("\n")
##     invisible(x)
## }

## function (x, digits = max(3L, getOption("digits") - 2L), signif.stars = getOption("show.signif.stars"),
##     signif.legend = signif.stars, dig.tst = max(1L, min(5L, digits -
##         1L)), cs.ind = 1:k, tst.ind = k + 1, zap.ind = integer(),
##     P.values = NULL, has.Pvalue = nc >= 4 && substr(colnames(x)[nc],
##         1, 3) == "Pr(", eps.Pvalue = .Machine$double.eps, na.print = "NA",
##     ...)
## {
##     if (is.null(d <- dim(x)) || length(d) != 2L)
##         stop("'x' must be coefficient matrix/data frame")
##     nc <- d[2L]
##     if (is.null(P.values)) {
##         scp <- getOption("show.coef.Pvalues")
##         if (!is.logical(scp) || is.na(scp)) {
##             warning("option \"show.coef.Pvalues\" is invalid: assuming TRUE")
##             scp <- TRUE
##         }
##         P.values <- has.Pvalue && scp
##     }
##     else if (P.values && !has.Pvalue)
##         stop("'P.values' is TRUE, but 'has.Pvalue' is not")
##     if (has.Pvalue && !P.values) {
##         d <- dim(xm <- data.matrix(x[, -nc, drop = FALSE]))
##         nc <- nc - 1
##         has.Pvalue <- FALSE
##     }
##     else xm <- data.matrix(x)
##     k <- nc - has.Pvalue - (if (missing(tst.ind))
##         1
##     else length(tst.ind))
##     if (!missing(cs.ind) && length(cs.ind) > k)
##         stop("wrong k / cs.ind")
##     Cf <- array("", dim = d, dimnames = dimnames(xm))
##     ok <- !(ina <- is.na(xm))
##     for (i in zap.ind) xm[, i] <- zapsmall(xm[, i], digits)
##     if (length(cs.ind)) {
##         acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
##         if (any(ia <- is.finite(acs))) {
##             digmin <- 1 + if (length(acs <- acs[ia & acs != 0]))
##                 floor(log10(range(acs[acs != 0], finite = TRUE)))
##             else 0
##             Cf[, cs.ind] <- format(round(coef.se, max(1L, digits -
##                 digmin)), digits = digits)
##         }
##     }
##     if (length(tst.ind))
##         Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst),
##             digits = digits)
##     if (any(r.ind <- !((1L:nc) %in% c(cs.ind, tst.ind, if (has.Pvalue) nc))))
##         for (i in which(r.ind)) Cf[, i] <- format(xm[, i], digits = digits)
##     ok[, tst.ind] <- FALSE
##     okP <- if (has.Pvalue)
##         ok[, -nc]
##     else ok
##     x1 <- Cf[okP]
##     dec <- getOption("OutDec")
##     if (dec != ".")
##         x1 <- chartr(dec, ".", x1)
##     x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)
##     if (length(not.both.0 <- which(x0 & !is.na(x0)))) {
##         Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(1L,
##             digits - 1L))
##     }
##     if (any(ina))
##         Cf[ina] <- na.print
##     if (P.values) {
##         if (!is.logical(signif.stars) || is.na(signif.stars)) {
##             warning("option \"show.signif.stars\" is invalid: assuming TRUE")
##             signif.stars <- TRUE
##         }
##         if (any(okP <- ok[, nc])) {
##             pv <- as.vector(xm[, nc])
##             Cf[okP, nc] <- format.pval(pv[okP], digits = dig.tst,
##                 eps = eps.Pvalue)
##             signif.stars <- signif.stars && any(pv[okP] < 0.1)
##             if (signif.stars) {
##                 Signif <- symnum(pv, corr = FALSE, na = FALSE,
##                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
##                   symbols = c("***", "**", "*", ".", " "))
##                 Cf <- cbind(Cf, format(Signif))
##             }
##         }
##         else signif.stars <- FALSE
##     }
##     else signif.stars <- FALSE
##     print.default(Cf, quote = FALSE, right = TRUE, na.print = na.print,
##         ...)
##     if (signif.stars && signif.legend)
##         cat("---\nSignif. codes:  ", attr(Signif, "legend"),
##             "\n", sep = "")
##     invisible(x)
## }
