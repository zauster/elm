
minTypeII <- function(betaj, X, ww, XROWS, XCOLS, ej, tau_jB,
                      tauj_2, tau_j, tauj_inf, betabarj, tbarmin,
                      alpha, ds, b,
                      taumm_jB, taujmm_2, taujmm_inf, tbarminmm,
                      dsmm, bmm, root = FALSE)
    {
        ## OLS Nonstandardized
        sigmasqbarOLS <- calcSigmasqbar(X, ww, XROWS, XCOLS, ej,
                                        tau_jB, tauj_2, tau_j,
                                        betabarj, betaj,
                                        type = "typeII")
        ## cat(sigmasqbarOLS)
        ## cat("\nbeta_j: ", betaj)

        OLSNonstandardizedTypeII <- calcTypeIINonstandardized(wb1start = c(0.1, 0.1),
                                                              lowerBE = rep(10^-6, 2),
                                                              sigmasqbar = sigmasqbarOLS["TypeII"],
                                                              betaj = betaj,
                                                              betabarj = betabarj,
                                                              tbarmin = tbarmin,
                                                              tauj_2 = tauj_2,
                                                              tauj_inf = tauj_inf)
        OLSNonstandardizedTypeII <- min(OLSNonstandardizedTypeII)
        names(OLSNonstandardizedTypeII) <- "Nonstandardized (OLS)"

        ## OLS Bernoulli
        OLSBernoulliTypeII <- calcTypeIIBernoulli(betaj = betaj,
                                                  betabarj = betabarj,
                                                  alpha = alpha,
                                                  ds = ds, b = b,
                                                  XROWS = XROWS)$typeII
        names(OLSBernoulliTypeII) <- "Bernoulli (OLS)"

        ## MM Nonstandardized
        sigmasqbarMM <- calcSigmasqbar(X, ww, XROWS, XCOLS, ej,
                                       taumm_jB, taujmm_2, taujmm_inf,
                                       betabarj, betaj,
                                        type = "typeII")
        MMNonstandardizedTypeII <- calcTypeIINonstandardized(wb1start = c(0.1, 0.1),
                                                             lowerBE = rep(10^-6, 2),
                                                             sigmasqbar = sigmasqbarMM["TypeII"],
                                                             betaj = betaj,
                                                             betabarj = betabarj,
                                                             tbarmin = tbarminmm,
                                                             tauj_2 = taujmm_2,
                                                             tauj_inf = taujmm_inf)
        MMNonstandardizedTypeII <- min(MMNonstandardizedTypeII)
        names(MMNonstandardizedTypeII) <- "Nonstandardized (MM)"

        ## MM Bernoulli
        MMBernoulliTypeII <- calcTypeIIBernoulli(betaj = betaj,
                                                 betabarj = betabarj,
                                                 alpha = alpha,
                                                 ds = dsmm, b = bmm,
                                                 XROWS = XROWS)$typeII
        names(MMBernoulliTypeII) <- "Bernoulli (MM)"

        res <- c(OLSNonstandardizedTypeII,
                   OLSBernoulliTypeII,
                   MMNonstandardizedTypeII,
                   MMBernoulliTypeII)

        if(root == TRUE)
            {
                res <- res[which.min(res)] - 0.5
                ## res <- (res[which.min(res)] - 0.5)^2
            }


        return(res)
    }

findMinTypeII <- function(upperbetabound, X, ww, XROWS, XCOLS, ej, tau_jB,
                      tauj_2, tau_j, tauj_inf, betabarj, tbarmin,
                      alpha, ds, b,
                      taumm_jB, taujmm_2, taujmm_inf, tbarminmm,
                      dsmm, bmm, silent, root = TRUE)
    {
        res <- try(uniroot(minTypeII, c(betabarj, upperbetabound), X, ww, XROWS, XCOLS, ej,
                       tau_jB, tauj_2, tau_j,
                       tauj_inf, betabarj, tbarmin,
                       alpha, ds, b,
                       taumm_jB, taujmm_2,
                       taujmm_inf, tbarminmm, dsmm, bmm, root = root,
                       tol = .Machine$double.eps^0.5),
                   silent = silent)

        if(is.error(res))
            {
                res <- findMinTypeII(upperbetabound * 1.1,
                                     X, ww, XROWS, XCOLS, ej, tau_jB,
                      tauj_2, tau_j, tauj_inf, betabarj, tbarmin,
                      alpha, ds, b,
                      taumm_jB, taujmm_2, taujmm_inf, tbarminmm,
                      dsmm, bmm, silent, root = TRUE)
            }

        res
    }

## findUpperLimit <- function(lower, upper, steps = 25, X, ww, XROWS, XCOLS, ej,
##                            tau_jB, tauj_2, tau_j,
##                            betaj, type = "typeI",
##                            zero = 10^-12)
##     {
##         no.error <- TRUE
##         limit <- lower
##         j <- 1
##                 range <- seq(lower, upper, length.out = steps)
##                 print(range)
##         while(no.error == TRUE && j <= steps)
##             {
##                 res <- calcSigmasqbar2(range[j], X = X, ww = ww,
##                                   XROWS = XROWS, XCOLS = XCOLS, ej = ej,
##                                   tau_jB = tau_jB, tauj_2 = tauj_2,
##                                   tau_j = tau_j,
##                                   betabarj = betabarj,
##                                   type = "typeII")
##                 cat("\nrange: ", range[j])
##                 cat("\tres", res)
##                 no.error <- !is.error(res)
##                 j <- j + 1
##                 ## if(sum(res) == steps)
##                 ##     {
##                 ##         lower <- upper
##                 ##         limit <- upper
##                 ##         upper <- upper + range
##                 ##         no.error <- TRUE
##                 ##     } else {
##                 ##         limit <- limit + range[sum(res)]
##                 ##         no.error <- FALSE
##                 ##     }
##             }
##         limit <- range[j - 1]
##     }

findHighestBeta <- function(Y, X, j)
    {
        ## n <- nrow(X)
        ## Y <- seq(min(Y), max(Y), length.out = n)
        ## limits <- rbind(apply(X, 2, min), apply(X, 2, max))
        ## X <- apply(limits, 2, function(mat) seq(mat[1], mat[2], length.out = n))

        ## print(summary(lm(Y ~ X)))
        ## res <- lm(Y ~ X)$coefficients

        Y <- c(min(Y), max(Y))
        X <- c(min(X[, j]), max(X[, j]))

        res <- lm(Y ~ X)$coefficients[2]
        res
    }

is.error <- function(x)
    {
        inherits(x, "try-error")
    }

findMinTypeII2 <- function(upperbetabound, X, ww, XROWS, XCOLS, ej, tau_jB,
                      tauj_2, tau_j, tauj_inf, betabarj, tbarmin,
                      alpha, ds, b,
                      taumm_jB, taujmm_2, taujmm_inf, tbarminmm,
                      dsmm, bmm, root = TRUE)
    {
        res <- optimize(minTypeII, c(betabarj, upperbetabound),
                        X, ww, XROWS, XCOLS, ej,
                       tau_jB, tauj_2, tau_j,
                       tauj_inf, betabarj, tbarmin,
                       alpha, ds, b,
                       taumm_jB, taujmm_2,
                       taujmm_inf, tbarminmm, dsmm, bmm, root = root,
                       tol = .Machine$double.eps^0.5)

        res
    }





## res <- NULL
## for(x in 1:99/100)
##     {
        ## res <- rbind(res, minTypeII(x, X, ww, XROWS, XCOLS, ej,
        ##                             tau_jB, tauj_2, tau_j,
        ##                             tauj_inf, betabarj, tbarmin,
        ##                             alpha, ds, b,
        ##                             taumm_jB, taujmm_2,
        ##                             taujmm_inf, tbarminmm, dsmm, bmm))
##     }
## matplot(res, type = "l")
