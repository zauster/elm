
testCoefficient <- function(j, Y, X, ww,
                            betahat, betabarj, alpha,
                            alternative, XROWS, XCOLS,
                            tau, iterations,
                            qq, qqmm, lambda, lambdamm)
{
    COEFNAME <- ifelse(!is.null(colnames(X)[j]),
                       colnames(X)[j],
                       paste("beta_", j, sep = ""))
    nullHypothesis <- paste(COEFNAME, ifelse(alternative == "greater",
                                             "<=", ">="),
                            betabarj)
    betahatj <- betahat[j]
    tau_j <- tau[,j]
    tauj_2 <- sum((tau_j)^2)
    tauj_inf <- max(abs(tau_j))

    ## tauOLS <- list(j = tau[, j],
    ##                sq = sum(tau[, j]^2),
    ##                inf = max(abs(tau[, j])))

    ej <- rep(0, times = XCOLS)
    ej[j] <- 1

    ## Finding upper bound on standard deviation of estimator:
    ## adjusting tau_j st Dmat positive definite,
    ## rescale and then eliminate 0 entries
    ## increases bound so test remains exact

    tauj_inf_min001 <- min(0.01, tauj_inf) ## less computations
    tauj_inf_min1 <- min(1, tauj_inf)

    tau_jB <- sapply(1:XROWS,
                     function(i) ifelse(tau_j[i] >= 0,
                                        max(tau_j[i]/tauj_inf_min001,
                                            qq * tauj_inf_min1),
                                        min(tau_j[i]/tauj_inf_min001,
                                            -qq * tauj_inf_min1)))


    ##
    ## MM Estimate
    ##
    obj <- c(rep(0,XROWS),1)
    mat <- matrix(rbind(cbind(diag(XROWS), rep(-1, XROWS)),
                        cbind(diag(XROWS), rep(1, XROWS)),
                        cbind(t(X), rep(0, XCOLS)),
                        cbind(t(rep(0, XROWS)), 1)),
                  nrow = 2 * XROWS + XCOLS + 1)
    dir <- c(rep("<=", XROWS), rep(">=", XROWS), rep("==", XCOLS), ">=")
    rhs <- c(rep(0, 2 * XROWS), ej, 0)
    XROWSplusone <- XROWS + 1
    sol <- Rglpk_solve_LP(obj, mat, dir, rhs, types = NULL, max = FALSE,
                          bounds = list(lower = list(ind = 1:XROWSplusone,
                                            val = rep(-XROWS,XROWSplusone)),
                              upper = list(ind = 1:XROWSplusone,
                                  val = rep(XROWS,XROWSplusone))),
                          verbose = FALSE)

    taumm_j <- sol$solution
    taumm_j <- taumm_j[1:XROWS]
    betahatjmm = crossprod(taumm_j, Y)
    taujmm_2 = sum((taumm_j)^2)
    taujmm_inf = max(abs(taumm_j))

    ## adjusting tau_j st Dmat positive definite,
    ## rescale and then eliminate 0 entries
    ## increases bound so test remains exact

    taumm_jB <- rep(0,times = XROWS)
    taujmm_2_min <- min(1, taujmm_2) ## less computations

    taumm_jB <- sapply(1:XROWS,
                       function(i) ifelse(taumm_j[i] >= 0,
                                          max(taumm_j[i]/taujmm_2_min,
                                              qqmm * taujmm_2_min),
                                          min(taumm_j[i]/taujmm_2_min,
                                              -qqmm * taujmm_2_min)))

    b <- XROWS * tauj_inf
    d <- sapply(1:XROWS,
                function(i) (1 - lambda) * max(-tau_j[i] * ww,
                                               -tau_j[i] * (ww + 1)) + lambda * (tauj_inf - max(tau_j[i] * ww, tau_j[i] * (ww + 1))))
    ds <- sum(d)

    bmm <- XROWS * taujmm_inf
    dmm <- sapply(1:XROWS, function(i) (1 - lambdamm) * max(-taumm_j[i] * ww,
                                                            -taumm_j[i] * (ww + 1)) + lambdamm * (taujmm_inf - max(taumm_j[i] * ww, taumm_j[i] * (ww + 1))))
    dsmm <- sum(dmm)

    ## ## ## ##
    ## ## ## ## Tests
    ## ## ## ##

    ##
    ## Nonstandardized Test (OLS)
    ##
    sigmasqbarOLS <- calcSigmasqbar(X = X,
                                    ww = ww,
                                    XROWS = XROWS,
                                    XCOLS = XCOLS,
                                    ej = ej,
                                    tau_jB = tau_jB,
                                    tauj_2 = tauj_2,
                                    tau_j = tau_j,
                                    betabarj = betabarj,
                                    betaj = betaj,
                                    type = "typeI")
    OLSNonstandardizedTests <- calcNonstandardizedTest(tauj_2 = tauj_2,
                                                       tauj_inf = tauj_inf,
                                                       alpha = alpha,
                                                       sigmasqbar = sigmasqbarOLS["TypeI"])
    tbarOLS <- min(OLSNonstandardizedTests)


    ##
    ## Nonstandardized Tests (MM)
    ##
    sigmasqbarMM <- calcSigmasqbar(X = X,
                                   ww = ww,
                                   XROWS = XROWS,
                                   XCOLS = XCOLS,
                                   ej = ej,
                                   tau_jB = taumm_jB,
                                   tauj_2 = taujmm_2,
                                   tau_j = taumm_j,
                                   betabarj = betabarj,
                                   betaj = betaj,
                                   type = "typeI")
    MMNonstandardizedTests <- calcNonstandardizedTest(tauj_2 = taujmm_2,
                                                      tauj_inf = taujmm_inf,
                                                      alpha = alpha,
                                                      sigmasqbar = sigmasqbarMM["TypeI"])
    tbarMM <- min(OLSNonstandardizedTests)


    ##
    ## TypeII Optimization
    ##
    ## find betaj that brings typeII to 0.5
    optbetaj <- findMinTypeII(X = X, ww = ww, XROWS = XROWS, XCOLS = XCOLS,
                              ej = ej, tau_jB = tau_jB, tauj_2 = tauj_2,
                              tau_j = tau_j, tauj_inf = tauj_inf,
                              betabarj = betabarj,
                              tbarmin = tbarOLS,
                              alpha = alpha, ds = ds, b = b,
                              taumm_jB = taumm_jB, taujmm_2 = taujmm_2,
                              taujmm_inf = taujmm_inf,
                              tbarminmm = tbarMM,
                              dsmm = dsmm, bmm = bmm)$root

    TypeII <- minTypeII(betaj = optbetaj,
                        X = X, ww = ww, XROWS = XROWS, XCOLS = XCOLS,
                        ej = ej, tau_jB = tau_jB, tauj_2 = tauj_2,
                        tau_j = tau_j, tauj_inf = tauj_inf,
                        betabarj = betabarj,
                        tbarmin = tbarOLS,
                        alpha = alpha, ds = ds, b = b,
                        taumm_jB = taumm_jB, taujmm_2 = taujmm_2,
                        taujmm_inf = taujmm_inf,
                        tbarminmm = tbarMM,
                        dsmm = dsmm, bmm = bmm)

    sigmasqbarOLS <- calcSigmasqbar(X, ww, XROWS, XCOLS, ej,
                                    tau_jB, tauj_2, tau_j,
                                    betabarj, betaj = optbetaj,
                                    type = "typeII")
    OLSNonstandardizedTypeII <- calcTypeIINonstandardized(wb1start = c(0.1, 0.1),
                                                          lowerBE = rep(10^-6, 2),
                                                          sigmasqbar = sigmasqbarOLS["TypeII"],
                                                          betaj = optbetaj,
                                                          betabarj = betabarj,
                                                          tbarmin = tbarOLS,
                                                          tauj_2 = tauj_2,
                                                          tauj_inf = tauj_inf)
    OLSBernoulliTypeII <- calcTypeIIBernoulli(betaj = optbetaj,
                                              betabarj = betabarj,
                                              alpha = alpha,
                                              ds = ds, b = b,
                                              XROWS = XROWS)

    OLSBernoulliTest <- calcBernoulliTest(Y = Y, XROWS = XROWS,
                                          tauj_inf = tauj_inf,
                                          tau_j = tau_j, ww = ww,
                                          Bernoulliparameter = OLSBernoulliTypeII,
                                          d = d, ds = ds, b = b,
                                          betabarj = betabarj,
                                          alternative = alternative,
                                          iterations = iterations)

    sigmasqbarMM <- calcSigmasqbar(X, ww, XROWS, XCOLS, ej,
                                   taumm_jB, taujmm_2, taujmm_inf,
                                   betabarj, betaj = optbetaj,
                                   type = "typeII")
    MMNonstandardizedTypeII <- calcTypeIINonstandardized(wb1start = c(0.1, 0.1),
                                                         lowerBE = rep(10^-6, 2),
                                                         sigmasqbar = sigmasqbarMM["TypeII"],
                                                         betaj = optbetaj,
                                                         betabarj = betabarj,
                                                         tbarmin = tbarMM,
                                                         tauj_2 = taujmm_2,
                                                         tauj_inf = taujmm_inf)

    MMBernoulliTypeII <- calcTypeIIBernoulli(betaj = optbetaj,
                                             betabarj = betabarj,
                                             alpha = alpha,
                                             ds = dsmm, b = bmm,
                                             XROWS = XROWS)
    MMBernoulliTest <- calcBernoulliTest(Y = Y, XROWS = XROWS,
                                         tauj_inf = taujmm_inf,
                                         tau_j = taumm_j, ww = ww,
                                         Bernoulliparameter = MMBernoulliTypeII,
                                         d = dmm, ds = dsmm, b = bmm,
                                         betabarj = betabarj,
                                         alternative = alternative,
                                         iterations = iterations)

    minimizingTest <- which.min(TypeII)
    results <- switch(minimizingTest,
                      {c(tbarOLS,
                         betahatj - betabarj)
                   },
                      {c(OLSBernoulliTest$theta,
                         OLSBernoulliTest$prob_rejection)
                   },
                      {c(tbarMM,
                         betahatj - betabarj)
                   },
                      {c(MMBernoulliTest$theta,
                         MMBernoulliTest$prob_rejection)})

    chosenTest <- list(nullHypothesis,
                       results[1],
                       results[2],
                       betahatj,
                       ifelse(results[1] < results[2], "Yes", "No"),
                       names(TypeII)[minimizingTest])
    if(minimizingTest == 1 | minimizingTest == 3)
        {
            names(chosenTest) <- c("H_0", "tbar", "Test Statistic",
                                   ifelse(minimizingTest == 1,
                                          "(OLS) Estimate",
                                          "(MM) Estimate"),
                                   "Rejection", "chosen Test")
}
else
    {
        names(chosenTest) <- c("H_0", "theta", "Prob rejection",
                               ifelse(minimizingTest == 2,
                                      "(OLS) Estimate",
                                      "(MM) Estimate"),
                               "Rejection", "chosen Test")
        }

    OLSNonstandardized <- list(NonstandardizedTests = OLSNonstandardizedTests,
                               NonstandardizedTypeII = OLSNonstandardizedTypeII)
    OLSBernoulli <- list(BernoulliTest = OLSBernoulliTest,
                         BernoulliTypeII = OLSBernoulliTypeII)
    MMNonstandardized <- list(NonstandardizedTests = MMNonstandardizedTests,
                              NonstandardizedTypeII = MMNonstandardizedTypeII)
    MMBernoulli <- list(BernoulliTest = MMBernoulliTest,
                        BernoulliTypeII = MMBernoulliTypeII)

    res <- list(j = j,
                betabarj = betabarj,
                betaj = optbetaj,
                tbars = c(tbarOLS = tbarOLS,
                    tbarMM = tbarMM),
                chosenTest = chosenTest,
                typeIImatrix = TypeII,
                OLSNonstandardized = OLSNonstandardized,
                OLSBernoulli = OLSBernoulli,
                MMNonstandardized = MMNonstandardized,
                MMBernoulli = MMBernoulli)
    res
}
