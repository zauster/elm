##
##
##         Confidence intervals for exact linear models
##
## data: Y and X
## n = 195, m = 7
## Number of tested coefficients: 1
##
## ----------- Coefficient: distance ----------------------
##             estimate           lower             upper
##  distance     0.1534         -0.1023            0.3211
##

print.elm <- function(x, ..., verbose = NULL,
                      digits = max(3L, getOption("digits") - 3L))
    {
        cat("\n")
        cat(strwrap(x$method, prefix = "\t"), sep="\n")
        cat("\n")
        cat("data:", x$yname, "and", x$xname, "\n")
        cat("n = ", x$parameter$n, ", m = ", x$parameter$m, "\n", sep = "")
        cat("Number of tested coefficients:", length(x$coefTests), "\n")

        if(!is.null(verbose) & is.logical(verbose))
            {
                x$verbose <- verbose
            }

        for(lst in x$coefTests)
            {
                for(i in 2:4)
                    {
                        lst$chosenTest[[i]] <- format(round(lst$chosenTest[[i]],
                                                            digits = digits),
                                              digits = digits)
                    }
                lst$chosenTest$Rejection <- ifelse(lst$chosenTest$Rejection == TRUE,
                                                   "Yes", "No")
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
                        cat("\nOLS based: \t(estimate = ",
                            format(lst$betahatj[1], digits = digits), ")\n", sep = "")
                        cat("Nonstandardized Type I: \n")
                        print(format(lst$OLSNonstandardized$NonstandardizedTests,
                                     digits = digits), quote = FALSE, right = TRUE)
                        cat("Nonstandardized Type II: \n")
                        print(format(lst$OLSNonstandardized$NonstandardizedTypeII,
                                     digits = digits), quote = FALSE, right = TRUE)

                        cat("\nBernoulli:\n")
                        print(format(lst$OLSBernoulli$BernoulliTest,
                                     digits = digits), quote = FALSE, right = TRUE)

                        cat("\nMM based: \t(estimate = ",
                            format(lst$betahatj[2], digits = digits), ")\n", sep = "")
                        cat("Nonstandardized Type I: \n")
                        print(format(lst$MMNonstandardized$NonstandardizedTests,
                                     digits = digits), quote = FALSE, right = TRUE)
                        cat("Nonstandardized Type II: \n")
                        print(format(lst$MMNonstandardized$NonstandardizedTypeII,
                                     digits = digits), quote = FALSE, right = TRUE)

                        cat("\nBernoulli:\n")
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


##      Confidence intervals for exact linear models

## data: y and x
## n = 876, m = 21
## Number of tested coefficients: 1

##                                                  lower       upper
##     safi_lr04      Nonstandardized Test          0.066       0.234
##                                                  lower       upper
##     safi_lr04      Nonstandardized Test          0.066       0.234

print.elmCI <- function(x, ..., verbose = NULL,
                      digits = max(3L, getOption("digits") - 3L))
    {
        cat("\n")
        cat(strwrap(x$method, prefix = "\t"), sep="\n")
        cat("\n")
        cat("data:", x$yname, "and", x$xname, "\n")
        cat("n = ", x$n, ", m = ", x$m, "\n", sep = "")
        cat("Number of tested coefficients:", length(x$CI), "\n")

        if(!is.null(verbose) & is.logical(verbose))
            {
                x$verbose <- verbose
            }

        CIs <- matrix(numeric(0), ncol = 4)
        for(lst in x$CI)
            {
                ## cat("\n")
                lst[[2]] <- format(round(lst[[2]], digits = digits), digits = digits)
                lst <- unlist(lst)
                names(lst) <- NULL
                ## print(lst, quote = FALSE, right = TRUE)
                CIs <- rbind(CIs, lst)
            }

        rownames(CIs) <- rep(c(""), nrow(CIs))
        colnames(CIs) <- c("Coefficient", "lower", "upper", "used Test")

        cat("\n")
        print(CIs, quote = FALSE, right = TRUE)
        ## cat("\nparameters:\n")
        ## if(!is.null(x$parameter$bounds))
        ##     {
        ##         cat(paste("   ", x$yname, " in ", x$parameter$bounds, "\n", sep = ""))
        ##     }
        ## cat("   alpha:", x$parameter$alpha)
        ## if(!is.null(x$parameter$iterations))
        ##     {
        ##         cat("\n   iterations:", x$parameter$iterations)
        ##     }
        cat("\n")

    }
