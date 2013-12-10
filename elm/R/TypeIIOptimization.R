
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
        names(OLSBernoulliTypeII) <- "    Bernoulli (OLS)"

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
        names(MMBernoulliTypeII) <- "    Bernoulli (MM)"

        res <- c(OLSNonstandardizedTypeII,
                   OLSBernoulliTypeII,
                   MMNonstandardizedTypeII,
                   MMBernoulliTypeII)

        if(root == TRUE)
            {
                res <- res[which.min(res)] - 0.5
            }


        return(res)
    }

findMinTypeII <- function(X, ww, XROWS, XCOLS, ej, tau_jB,
                      tauj_2, tau_j, tauj_inf, betabarj, tbarmin,
                      alpha, ds, b,
                      taumm_jB, taujmm_2, taujmm_inf, tbarminmm,
                      dsmm, bmm, root = TRUE)
    {
        res <- uniroot(minTypeII, c(0.0001, 0.99), X, ww, XROWS, XCOLS, ej,
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
