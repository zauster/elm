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


print.elm <- function(x, ..., digits = 5)
    {
        cat("\n")
        cat(strwrap(x$method, prefix = "\t"), sep="\n")
        cat("\n")
        cat("data:", x$yname, "and", x$xname, "\n")
        cat("n =", x$parameter$n, ", m =", x$parameter$m)
        cat("\n")

        cat("Tested coefficients:")
        ## for(coef in coefficients)


        ## if(x$printdetails == TRUE)
        ##     {
        ##         cat("\nnot yet available\n")
        ##     }

        cat("estimated OLS coefficient of beta_", x$parameter$j, ":",
            "\t", x$betahatjOLS, "\n", sep = "")
        cat("estimated MM coefficient of beta_", x$parameter$j, ":",
            "\t", x$betahatjMM, "\n", sep = "")
        cat("TypeII minimizing betaj:", x$parameter$optbetaj, "\n")

        ## cat(paste("\treject H_0: beta_", x$parameter$j,
        ##           ifelse(x$parameter$alternative == "<=", " <= ", " >= "),
        ##           x$betabarj, "?", sep = ""),
        ##     "\t", ifelse(x$rejection == TRUE, "Yes", "No"), "\n")
        ## cat("\t\tbased on:\t\t", x$chosentest$test, "\n")
        ## cat("\t\tEstimator:\t\t", x$chosentest$estimator, "\n")
        ## cat("\t\ttheta:\t\t\t", x$chosentest$theta, "\n")
        ## cat(paste("TypeII bound for beta_", x$j, ":", sep = ""),
        ##     "\t\t", x$chosentest$typeII, "\n")

        cat("\nOLS Nonstandardized:\n")
        print(x$OLSNonstandardized)
        cat("\n")

        cat("OLS Bernoulli:\n")
        print(x$OLSBernoulli$BernoulliTest)
        cat("Bernoulli TypeII: ", x$OLSBernoulli$BernoulliTypeII$typeII)
        cat("\n")

        cat("\nMM Nonstandardized:\n")
        print(x$MMNonstandardized)
        cat("\n")

        cat("MM Bernoulli:\n")
        print(x$MMBernoulli$BernoulliTest)
        cat("Bernoulli TypeII: ", x$MMBernoulli$BernoulliTypeII$typeII)
        cat("\n")

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

printCoefs <- function(x) ## x being a list of results from a test of
    ## a coefficient
    {

    }
