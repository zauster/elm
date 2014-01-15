
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
                alpha = 0.05,
                coefs = 2, ## coefficients to be tested
                nullvalue = 0, ## the threshold value in the null hypothesis
                upperbetabound = 1,
                lambda = 1, lambdamm = 1,
                qq = 0.0001, qqmm = 0.0001,
                iterations = 1000,
                steppc = 0.1,
                intercept = TRUE,
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

    X <- data.matrix(X)
    Y <- as.vector(Y)

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

    if(any(apply(X, 2, is.constant1)) == FALSE & intercept == TRUE)
        {
            warning("No intercept found, thus included.")
            X <- cbind(1, X)
        }

    XROWS <- nrow(X)
    XCOLS <- ncol(X)
    bounds <- paste("[", lower, ", ", upper, "]", sep = "")

    if(min(Y) < lower | max(Y) > upper)
        {
            stop("Some values of Y are not within the range [lower, upper]!")
        }

    if(length(Y) != XROWS)
        {
            stop("X and Y are of unequal length.")
        }

    m <- length(coefs)
    if(m != length(nullvalue) & length(nullvalue) > 1)
        {
            stop("Number of tested coefficients and null values does not correspond!")
        }

    if(m > 1 & length(nullvalue) == 1)
        {
            nullvalue <- rep(nullvalue, times = m)
        }

    if(!is.null(upperbetabound) && upperbetabound == nullvalue)
        {
            stop("Please set a reasonable value for the upperbetabound.")
        }

    ## check if X has full rank
    if(det(t(X) %*% X) == 0)
        {
            stop("Some columns are not independent. We need to have full rank, please eliminate redundant columns.")
        }

    ## adjustments to deal with both inequalities
    if(alternative == "greater")
        {

        }
    else if(alternative == "less")
        {
            lower.new <- -1 * upper
            upper <- -1 * lower
            lower <- lower.new
            Y <- -1 * Y
            nullvalue <- -1 * nullvalue
        }
    else if(alternative != "two.sided")
        {
            stop("Please choose either alternative = 'greater', 'less' or 'two.sided'.")
        }

    YY <- Y /(upper - lower) ## puts Y in [ww,ww + 1] where ww = lower/(upper-lower)
    XX <- X /(upper - lower) ## rescales X so that beta remains unchanged
    ww <- lower/(upper - lower)

    ##
    ## OLS estimate
    ##
    betahat <- solve(crossprod(XX)) %*% crossprod(XX, YY) # TODO:
                                        # numerically stable?
    tau <- XX %*% solve(crossprod(XX))

    coefTests <- list()
    if(alternative != "two.sided")
        {
            for(j in coefs)
                {
                    if(is.null(upperbetabound))
                        {
                            betabound <- 0.9 * findHighestBeta(YY, XX, j, alternative)
                            ## cat("\nupper: ", betabound)
                        }
                    else
                        {
                            betabound <- upperbetabound[coefs == j]
                        }
                    coefTests[[length(coefTests) + 1L]] <- testCoefficient(j = j, Y = YY, X = XX,
                                                                           ww = ww,
                                                                           betahat = betahat,
                                                                           betabarj = nullvalue[coefs == j],
                                                                           alpha = alpha,
                                                                           upperbetabound = betabound,
                                                                           steppc = steppc,
                                                                           alternative = alternative,
                                                                           XROWS = XROWS,
                                                                           XCOLS = XCOLS, tau = tau,
                                                                           iterations = iterations,
                                                                           qq = qq, qqmm = qqmm,
                                                                           lambda = lambda,
                                                                           lambdamm = lambdamm,
                                                                           silent = silent)

                }
        }
    else
        {
            for(j in coefs)
                {
                    ## cat("\nTwo.sided")
                    ## cat("\nupper: ", upper, "\tlower: ", lower, "\n")

                    if(is.null(upperbetabound))
                        {
                            betabound <- 0.9 * findHighestBeta(YY, XX, j, alternative)
                            ## cat("\nupper: ", betabound)
                        }
                    else
                        {
                            betabound <- upperbetabound[coefs == j]
                        }
                    res.upper <- testCoefficient(j = j, Y = YY, X = XX,
                                                 ww = ww,
                                                 betahat = betahat,
                                                 betabarj = nullvalue[coefs == j],
                                                 alpha = alpha/2,
                                                 upperbetabound = betabound,
                                                 steppc = steppc,
                                                 alternative = alternative,
                                                 XROWS = XROWS,
                                                 XCOLS = XCOLS, tau = tau,
                                                 iterations = iterations,
                                                 qq = qq, qqmm = qqmm,
                                                 lambda = lambda,
                                                 lambdamm = lambdamm,
                                                 silent = silent)
                    lower.new <- -1 * upper
                    upper <- -1 * lower
                    lower <- lower.new
                    ## cat("\nupper: ", upper, "\tlower: ", lower, "\n")
                    Y <- -1 * Y
                    nullvalue <- -1 * nullvalue
                    YY <- Y /(upper - lower)
                    XX <- X /(upper - lower)
                    ww <- lower/(upper - lower)
                    betahat <- solve(crossprod(XX)) %*% crossprod(XX, YY) # TODO:
                                        # numerically stable?
                    tau <- XX %*% solve(crossprod(XX))
                    res.lower <- testCoefficient(j = j, Y = YY, X = XX, ww = ww,
                                                 betahat = betahat,
                                                 betabarj = nullvalue[coefs == j],
                                                 alpha = alpha/2,
                                                 upperbetabound = betabound,
                                                 steppc = steppc,
                                                 alternative = "less",
                                                 XROWS = XROWS,
                                                 XCOLS = XCOLS, tau = tau,
                                                 iterations = iterations,
                                                 qq = qq, qqmm = qqmm,
                                                 lambda = lambda,
                                                 lambdamm = lambdamm,
                                                 silent = silent)

                    if(res.upper$chosenTest$Rejection == res.lower$chosenTest$Rejection)
                        {
                            if(res.upper$chosenTest$'chosen Test' == res.upper$chosenTest$'chosen Test')
                                {
                                    cat("\nSame same\n")
                                    if(res.upper$chosenTest[[3]] > res.lower$chosenTest[[3]])
                                        {
                                            cat("\nUpper")
                                            res <- res.upper
                                        }
                                    else
                                        {
                                            cat("\nLower")
                                            res <- res.lower
                                        }
                                }
                            else
                                {
                                    cat("\nDifferent test for upper and lower bound!")
                                    res <- res.upper
                                }
                        }
                    else if(res.upper$chosenTest$Rejection == TRUE)
                        {
                            cat("\nUpper")
                            res <- res.upper
                        }
                    else
                        {
                            cat("\nLower")
                            res <- res.lower
                        }
                    res$chosenTest$H_0 <- gsub(">", "!", res$chosenTest$H_0)
                    coefTests[[length(coefTests) + 1L]] <- res
                }
        }
    ##
    ## end

    ## guestj7
    ## ecoguest
    method <- "Exact linear models"
    alpha <- ifelse(alternative == "two.sided",
                    paste(alpha/2, "(on both sides)"),
                    alpha)


    parameter <- list(n = XROWS,
                      m = XCOLS,
                      alternative = alternative,
                      j = j,
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
