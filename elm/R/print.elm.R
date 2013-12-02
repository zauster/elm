
##        SUMMARY RESULTS, n= 195 m= 7  Y range [ 60 , 80 ]
## 1 Estimated coefficient of beta 1 :             55.25884
## 2        Reject H0: beta 1 <= 100 ? NO  at alpha =  0.05
## 3                         based on:            Bernoulli
## 4                        Estimator:              Minimax
## 5                            theta:              0.24519
## 6 Type II bound for beta 1 >= 270 :              0.63771

print.elm <- function(x, digits = 5)
    {
        cat("\n")
        cat(strwrap(x$method, prefix = "\t"), sep="\n")
        cat("\n")
        cat("data:", x$yname, "and", x$xname, "\n")
        cat("n =", x$parameter$n, ", m =", x$parameter$m)
        cat("\n")

        ## if(x$printdetails == TRUE)
        ##     {
        ##         cat("\nnot yet available\n")
        ##     }

        cat("estimated coefficient of beta_", x$parameter$j, ":",
            "\t", x$betahatj, "\n", sep = "")
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
