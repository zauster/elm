
minTypeII <- function(betaj, X, ww, XROWS, XCOLS, ej, Tau_jB,
                      TAUj_2, Tau_j, TAUj_inf, betabarj, tbarmin,
                      alpha, ds, b,
                      Taumm_jB, TAUjmm_2, TAUjmm_inf, tbarminmm,
                      dsmm, bmm, root = FALSE)
    {
        ## OLS Nonstandardized
        sigmasqbarOLS <- calcSigmasqbar(X, ww, XROWS, XCOLS, ej,
                                        Tau_jB, TAUj_2, Tau_j,
                                        betabarj, betaj)
        OLSNonstandardizedTypeII <- calcTypeIINonstandardized(wb1start = c(0.1, 0.1),
                                                              lowerBE = rep(10^-6, 2),
                                                              sigmasqbar = sigmasqbarOLS["TypeII"],
                                                              betaj = betaj,
                                                              betabarj = betabarj,
                                                              tbarmin = tbarmin,
                                                              TAUj_2 = TAUj_2,
                                                              TAUj_inf = TAUj_inf)
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
                                       Taumm_jB, TAUjmm_2, TAUjmm_inf,
                                       betabarj, betaj)
        MMNonstandardizedTypeII <- calcTypeIINonstandardized(wb1start = c(0.1, 0.1),
                                                             lowerBE = rep(10^-6, 2),
                                                             sigmasqbar = sigmasqbarMM["TypeII"],
                                                             betaj = betaj,
                                                             betabarj = betabarj,
                                                             tbarmin = tbarminmm,
                                                             TAUj_2 = TAUjmm_2,
                                                             TAUj_inf = TAUjmm_inf)
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
            }


        return(res)
    }

findMinTypeII <- function(X, ww, XROWS, XCOLS, ej, Tau_jB,
                      TAUj_2, Tau_j, TAUj_inf, betabarj, tbarmin,
                      alpha, ds, b,
                      Taumm_jB, TAUjmm_2, TAUjmm_inf, tbarminmm,
                      dsmm, bmm, root = TRUE)
    {
        res <- uniroot(minTypeII, c(0.0001, 0.99), X, ww, XROWS, XCOLS, ej,
                       Tau_jB, TAUj_2, Tau_j,
                       TAUj_inf, betabarj, tbarmin,
                       alpha, ds, b,
                       Taumm_jB, TAUjmm_2,
                       TAUjmm_inf, tbarminmm, dsmm, bmm, root = root,
                       tol = .Machine$double.eps^0.5)

        res
    }




## res <- NULL
## for(x in 1:99/100)
##     {
        ## res <- rbind(res, minTypeII(x, X, ww, XROWS, XCOLS, ej,
        ##                             Tau_jB, TAUj_2, Tau_j,
        ##                             TAUj_inf, betabarj, tbarmin,
        ##                             alpha, ds, b,
        ##                             Taumm_jB, TAUjmm_2,
        ##                             TAUjmm_inf, tbarminmm, dsmm, bmm))
##     }
## matplot(res, type = "l")
