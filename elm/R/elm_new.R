
## Given Y=X*beta+error where there are no assumptions imposed on the
## errors, it tests the one sided
## hypothesis H0: betaj <= betabarj against H1: betaj > betabarj
## (coded alternative = "greater")
## where j is index of coefficient.

## It also tests H0: betaj >= betabarj against H1: betaj < betabarj.
## (coded alternative = "less")


elm.new <- function(Y, X, lower = 0, upper = 1,
                alternative = "greater",
                alpha = 0.05, j = 2,
                betabarj = 0, betaj = 0.5,
                lambda = 1, lambdamm = 1,
                iterations = 1000, qq = 0.0001, qqmm = 0.0001)
{
    ## check required packages
    ## require(Rglpk)
    ## require(quadprog)

    ## for prettier output
    YNAME <- deparse(substitute(Y))
    XNAME <- deparse(substitute(X))
    bounds <- paste("[", lower, ", ", upper, "]", sep = "")

    if(min(Y) < lower | max(Y) > upper)
        {
            stop("Some values of Y are not within the range [lower, upper]!")
        }

    if(ncol(X) == 1)
        {
            stop("There is only one covariate, please add a constant.")
        }

    ## check if X has full rank
    if(det(t(X) %*% X) == 0)
        {
            stop("Some columns are not independent. We need to have full rank, please eliminate redundant columns.")
        }

    ## adjustments to deal with both inequalities
    if(alternative == "greater")
        ## if(IE == "<=")
        {
            if(betaj <= betabarj)
                {
                    stop("Please choose a betaj > betabarj.")
                }
        }
    else if(alternative == "less")
        ## else if(IE == ">=")
        {
            lower <- -1 * upper
            upper <- -1 * lower
            Y <- -1 * Y
            betabarj <- -1 * betabarj
            betaj <- -1 * betaj
            if(betaj >= betabarj)
                {
                    stop("Please choose a betaj < betabarj.")
                }
        }
    else
        {
            stop("Please choose either alternative = 'greater' or 'less'.")
        }


    XROWS <- nrow(X)
    XCOLS <- ncol(X)

    Y <- Y /(upper - lower) ## puts Y in [ww,ww + 1] where ww = lower/(upper-lower)
    X <- X /(upper - lower) ## rescales X so that beta remains unchanged
    ww <- lower/(upper - lower)

    ##
    ## OLS estimate
    ##
    betahat <- solve((crossprod(X))) %*% crossprod(X, Y) # TODO:
                                        # numerically stable?
    betahatj <- betahat[j]

    Tau <- X %*% solve(crossprod(X))
    Tau_j <- Tau[,j]
    TAUj_2 <- sum((Tau_j)^2)
    TAUj_inf <- max(abs(Tau_j))

    ## Type 1 OLS
    ej <- rep(0, times = XCOLS)
    ej[j] <- 1

    ## Finding upper bound on standard deviation of estimator:
    ## adjusting tau_j st Dmat positive definite,
    ## rescale and then eliminate 0 entries
    ## increases bound so test remains exact

    ## Tau_jB <- rep(0, times = XROWS)
    TAUj_inf_min001 <- min(0.01, TAUj_inf) ## less computations
    TAUj_inf_min1 <- min(1, TAUj_inf)

    Tau_jB <- sapply(1:XROWS, function(i) ifelse(Tau_j[i] >= 0,
                                                 max(Tau_j[i]/TAUj_inf_min001,
                                                     qq * TAUj_inf_min1),
                                                 min(Tau_j[i]/TAUj_inf_min001,
                                                     -qq * TAUj_inf_min1)))

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

    Taumm_j <- sol$solution
    Taumm_j <- Taumm_j[1:XROWS]
    betahatjmm = crossprod(Taumm_j, Y)
    TAUjmm_2 = sum((Taumm_j)^2)
    TAUjmm_inf = max(abs(Taumm_j))

    ## adjusting tau_j st Dmat positive definite,
    ## rescale and then eliminate 0 entries
    ## increases bound so test remains exact

    Taumm_jB <- rep(0,times = XROWS)
    TAUjmm_2_min <- min(1, TAUjmm_2) ## less computations

    Taumm_jB <- sapply(1:XROWS, function(i) ifelse(Taumm_j[i] >= 0,
                                                   max(Taumm_j[i]/TAUjmm_2_min,
                                                       qqmm * TAUjmm_2_min),
                                                   min(Taumm_j[i]/TAUjmm_2_min,
                                                       -qqmm * TAUjmm_2_min)))

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
                                    Tau_jB = Tau_jB,
                                    TAUj_2 = TAUj_2,
                                    Tau_j = Tau_j,
                                    betabarj = betabarj,
                                    betaj = betaj)

    OLSNonstandardizedTests <- calcNonstandardizedTest(TAUj_2 = TAUj_2,
                                                       TAUj_inf = TAUj_inf,
                                                       alpha = alpha,
                                                       sigmasqbar = sigmasqbarOLS["TypeI"])


    ##
    ## Nonstandardized Tests (MM)
    ##
    sigmasqbarMM <- calcSigmasqbar(X = X,
                                   ww = ww,
                                   XROWS = XROWS,
                                   XCOLS = XCOLS,
                                   ej = ej,
                                   Tau_jB = Taumm_jB,
                                   TAUj_2 = TAUjmm_2,
                                   Tau_j = Taumm_j,
                                   betabarj = betabarj,
                                   betaj = betaj)

    MMNonstandardizedTests <- calcNonstandardizedTest(TAUj_2 = TAUjmm_2,
                                                      TAUj_inf = TAUjmm_inf,
                                                      alpha = alpha,
                                                      sigmasqbar = sigmasqbarMM["TypeI"])


    ##
    ## TypeII Optimization
    ##
    ## find betaj that brings typeII to 0.5
    optbetaj <- findMinTypeII(X = X, ww = ww, XROWS = XROWS, XCOLS = XCOLS,
                                   ej = ej, Tau_jB = Tau_jB, TAUj_2 = TAUj_2,
                                   Tau_j = Tau_j, TAUj_inf = TAUj_inf,
                                   betabarj = betabarj,
                                   tbarmin = min(OLSNonstandardizedTests),
                                   alpha = alpha, ds = ds, b = b,
                                   Taumm_jB = Taumm_jB, TAUjmm_2 = TAUjmm_2,
                                   TAUjmm_inf = TAUjmm_inf,
                                   tbarminmm = min(MMNonstandardizedTests),
                                   dsmm = dsmm, bmm = bmm)$root

    TypeII <- minTypeII(betaj = optbetaj,
                        X = X, ww = ww, XROWS = XROWS, XCOLS = XCOLS,
                        ej = ej, Tau_jB = Tau_jB, TAUj_2 = TAUj_2,
                        Tau_j = Tau_j, TAUj_inf = TAUj_inf,
                        betabarj = betabarj,
                        tbarmin = min(OLSNonstandardizedTests),
                        alpha = alpha, ds = ds, b = b,
                        Taumm_jB = Taumm_jB, TAUjmm_2 = TAUjmm_2,
                        TAUjmm_inf = TAUjmm_inf,
                        tbarminmm = min(MMNonstandardizedTests),
                        dsmm = dsmm, bmm = bmm)

    sigmasqbarOLS <- calcSigmasqbar(X, ww, XROWS, XCOLS, ej,
                                        Tau_jB, TAUj_2, Tau_j,
                                        betabarj, betaj = optbetaj)
    OLSNonstandardizedTypeII <- calcTypeIINonstandardized(wb1start = c(0.1, 0.1),
                                                          lowerBE = rep(10^-6, 2),
                                                          sigmasqbar = sigmasqbarOLS["TypeII"],
                                                          betaj = optbetaj,
                                                          betabarj = betabarj,
                                                          tbarmin = min(OLSNonstandardizedTests),
                                                          TAUj_2 = TAUj_2,
                                                          TAUj_inf = TAUj_inf)
    OLSBernoulliTypeII <- calcTypeIIBernoulli(betaj = optbetaj,
                                              betabarj = betabarj,
                                              alpha = alpha,
                                              ds = ds, b = b,
                                              XROWS = XROWS)
    b <- XROWS * TAUj_inf
    d <- sapply(1:XROWS, function(i) (1 - lambda) * max(-Tau_j[i] * ww,
                                                        -Tau_j[i] * (ww + 1)) + lambda * (TAUj_inf - max(Tau_j[i] * ww, Tau_j[i] * (ww + 1))))
    ds <- sum(d)

    OLSBernoulliTest <- calcBernoulliTest(Y = Y, XROWS = XROWS,
                                          TAUj_inf = TAUj_inf,
                                          Tau_j = Tau_j, ww = ww,
                                          kbar = OLSBernoulliTypeII$kbar,
                                          pbar = OLSBernoulliTypeII$pbar,
                                          alphabar = OLSBernoulliTypeII$alphabar,
                                          d = d,
                                          ds = ds, b = b,
                                          ## betaj = betaj,
                                          betabarj = betabarj,
                                          alternative = alternative,
                                          iterations = iterations,
                                          theta = OLSBernoulliTypeII$theta)

    sigmasqbarMM <- calcSigmasqbar(X, ww, XROWS, XCOLS, ej,
                                       Taumm_jB, TAUjmm_2, TAUjmm_inf,
                                       betabarj, betaj = optbetaj)
    MMNonstandardizedTypeII <- calcTypeIINonstandardized(wb1start = c(0.1, 0.1),
                                                         lowerBE = rep(10^-6, 2),
                                                         sigmasqbar = sigmasqbarMM["TypeII"],
                                                         betaj = optbetaj,
                                                         betabarj = betabarj,
                                                         tbarmin = min(MMNonstandardizedTests),
                                                         TAUj_2 = TAUjmm_2,
                                                         TAUj_inf = TAUjmm_inf)

    MMBernoulliTypeII <- calcTypeIIBernoulli(betaj = optbetaj,
                                             betabarj = betabarj,
                                             alpha = alpha,
                                             ds = dsmm, b = bmm,
                                             XROWS = XROWS)
    bmm <- XROWS * TAUjmm_inf
    dmm <- sapply(1:XROWS, function(i) (1 - lambdamm) * max(-Taumm_j[i] * ww,
                                                        -Taumm_j[i] * (ww + 1)) + lambdamm * (TAUjmm_inf - max(Taumm_j[i] * ww, Taumm_j[i] * (ww + 1))))
    dsmm <- sum(dmm)
    MMBernoulliTest <- calcBernoulliTest(Y = Y, XROWS = XROWS,
                                         TAUj_inf = TAUjmm_inf,
                                         Tau_j = Taumm_j, ww = ww,
                                         kbar = MMBernoulliTypeII$kbar,
                                         pbar = MMBernoulliTypeII$pbar,
                                         alphabar = MMBernoulliTypeII$alphabar,
                                         d = dmm,
                                         ds = dsmm, b = bmm,
                                         ## betaj = betaj,
                                         betabarj = betabarj,
                                         alternative = alternative,
                                         iterations = iterations,
                                         theta = MMBernoulliTypeII$theta)

    ##
    ## end
    ##
    method <- "Exact linear models"

    OLSNonstandardized <- list(NonstandardizedTests = OLSNonstandardizedTests,
                               NonstandardizedTypeII = OLSNonstandardizedTypeII)
    OLSBernoulli <- list(BernoulliTest = OLSBernoulliTest,
                         BernoulliTypeII = OLSBernoulliTypeII)

    MMNonstandardized <- list(Nonstandardized = MMNonstandardizedTests,
                              NonstandardizedTypeII = MMNonstandardizedTypeII)
    MMBernoulli <- list(BernoulliTest = MMBernoulliTest,
                        BernoulliTypeII = MMBernoulliTypeII)
    parameter <- list(n = XROWS,
                      m = XCOLS,
                      alternative = alternative,
                      j = j,
                      optbetaj = optbetaj,
                      alpha = alpha,
                      bounds = bounds,
                      iterations = iterations)

    structure(list(method = method,
                   yname = YNAME,
                   xname = XNAME,
                   betabarj = betabarj,
                   betaj = betaj,
                   betahatjOLS = betahatj,
                   betahatjMM = betahatjmm,
                   parameter = parameter,
                   ## theta = theta,
                   OLSNonstandardized = OLSNonstandardized,
                   OLSBernoulli = OLSBernoulli,
                   MMNonstandardized = MMNonstandardized,
                   MMBernoulli = MMBernoulli),
              class = "elm")
}
