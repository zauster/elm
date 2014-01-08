
## Given Y=X*beta+error where there are no assumptions imposed on the
## errors, it tests the one sided
## hypothesis H0: betaj <= betabarj against H1: betaj > betabarj
## (coded alternative = "greater")
## where j is index of coefficient.

## It also tests H0: betaj >= betabarj against H1: betaj < betabarj.
## (coded alternative = "less")

source("NonstandardizedTest.R")
source("BernoulliTest.R")
source("TypeIIOptimization.R")
source("testCoefficient.R")
source("Sigmasqbar.R")
source("miscfun.R")
source("print.elm.R")


elm <- function(Y, X, lower = 0, upper = 1,
                alternative = "greater",
                alpha = 0.05, coefs = 2,
                betabarj = 0,
                upperbetabound = 1,
                lambda = 1, lambdamm = 1,
                qq = 0.0001, qqmm = 0.0001,
                iterations = 1000,
                silent = FALSE, ## warning during search for optimal beta
                verbose = TRUE, ## results for all tests
                na.action = getOption("na.action"))
{
    ## check required packages
    ## require(Rglpk)
    ## require(quadprog)

    ## for prettier output
    YNAME <- deparse(substitute(Y))
    XNAME <- deparse(substitute(X))

    if(is.data.frame(X) == TRUE)
        {
            X <- data.matrix(X)
        }

    X <- as.matrix(X)
    Y <- as.vector(Y)

    if(any(apply(X, 2, is.constant1)) == FALSE)
        {
            warning("No intercept found, thus included.")
            X <- cbind(1, X)
        }

    bounds <- paste("[", lower, ", ", upper, "]", sep = "")

    if(min(Y) < lower | max(Y) > upper)
        {
            stop("Some values of Y are not within the range [lower, upper]!")
        }

    if(ncol(X) == 1)
        {
            stop("There is only one covariate, please add a constant.")
        }

    if(length(Y) != dim(X)[1])
        {
            stop("X and Y are of unequal length.")
        }

    ## NA options
    switch(na.action,
           na.omit = {complete <- complete.cases(cbind(Y, X))
                      Y <- Y[complete]
                      X <- X[complete, ]},
           na.fail = {na.fail(Y)
                      na.fail(X)},
           na.pass = {},
           na.exclude = {complete <- complete.cases(cbind(Y, X))
                         Y <- Y[complete]
                         X <- X[complete, ]})

    ## check if X has full rank
    if(det(t(X) %*% X) == 0)
        {
            stop("Some columns are not independent. We need to have full rank, please eliminate redundant columns.")
        }

    ## adjustments to deal with both inequalities
    if(alternative == "greater")
        ## if(IE == "<=")
        {
            ## if(betaj <= betabarj)
            ##     {
            ##         stop("Please choose a betaj > betabarj.")
            ##     }
        }
    else if(alternative == "less")
        ## else if(IE == ">=")
        {
            lower <- -1 * upper
            upper <- -1 * lower
            Y <- -1 * Y
            betabarj <- -1 * betabarj
            ## betaj <- -1 * betaj
            ## if(betaj >= betabarj)
            ##     {
            ##         stop("Please choose a betaj < betabarj.")
            ##     }
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
    betahat <- solve(crossprod(X)) %*% crossprod(X, Y) # TODO:
                                        # numerically stable?
    tau <- X %*% solve(crossprod(X))

    coefTests <- list()
    for(j in coefs)
        {
            coefTests[[length(coefTests) + 1L]] <- testCoefficient(j = j, Y = Y, X = X,
                                                                   ww = ww,
                                                                   betahat = betahat,
                                                                   betabarj = betabarj,
                                                                   alpha = alpha,
                                                                   upperbetabound = upperbetabound,
                                                                   alternative = alternative,
                                                                   XROWS = XROWS,
                                                                   XCOLS = XCOLS, tau = tau,
                                                                   iterations = iterations,
                                                                   qq = qq, qqmm = qqmm,
                                                                   lambda = lambda,
                                                                   lambdamm = lambdamm,
                                                                   silent = silent)
        }

    ##
    ## end
    ##
    method <- "Exact linear models"

    parameter <- list(n = XROWS,
                      m = XCOLS,
                      alternative = alternative,
                      j = j,
                      ## optbetaj = optbetaj,
                      alpha = alpha,
                      bounds = bounds,
                      iterations = iterations)

    structure(list(method = method,
                   yname = YNAME,
                   xname = XNAME,
                   verbose = verbose,
                   parameter = parameter,
                   coefTests = coefTests),
              class = "elm")
}
