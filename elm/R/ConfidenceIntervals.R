
## steps:

## take elm-object
## see what test
## if nonstandardized -> nonstandardized
## if bernoulli -> bernoulli
## iterate over betabarj to find the bound

elmCI <- function(elm, alpha = 0.05, coefs = 2)
    {
        for(lst in elm$coefTests)
            {
                switch(lst$chosenTest$'chosen Test',
                       'Nonstandardized (OLS)' = {
                           cat("1\n")
                           res <- findNonstandardizedCI(upperbetabound = lst$upperbetabound,
                                                 nullvalue = lst$nullvalue,
                                                 betahatj = lst$betahatj[1],
                                                 elm = lst$model,
                                                 alpha = alpha)$root
                           cat("\nCI bound: ", res, "\n")
                       },
                       'Nonstandardized (MM)' = {
                           cat("2")
                       },
                       'Bernoulli (OLS)' = {
                           cat("3")
                       },
                       'Bernoulli (MM)' = {
                           cat("4")
                       })
            }
    }

findNonstandardizedCI <- function(upperbetabound, nullvalue, betahatj,
                                  elm, alpha, iter = 0, step = 0.05 * upperbetabound)
    {
        betainterval <- c(-upperbetabound, nullvalue)
        ## cat(betainterval)

        res <- try(uniroot(calcNonstandardizedCI,
                           interval = betainterval,
                           betahatj = betahatj,
                           elm = elm,
                           alpha = alpha),
                   silent = TRUE)

        if(is.error(res) == TRUE)
            {
                if(grepl("uniroot", res[1]) == TRUE)
                    {
                        iter <- iter + 1
                        res <- findNonstandardizedCI(upperbetabound = upperbetabound + step,
                                                     nullvalue = nullvalue,
                                                     betahatj = betahatj,
                                                     elm = elm, alpha = alpha,
                                                     iter = iter)
                    }
                else
                    {
                        print(res)
                        stop("Could not find CI bound.")
                    }
            }
        res
    }

calcNonstandardizedCI <- function(nullvalue, betahatj,
                                  elm, alpha)
    {
        sigmasqbar <- calcSigmasqbar(X = elm$X,
                                     ww = elm$ww,
                                     XROWS = elm$XROWS,
                                     XCOLS = elm$XCOLS,
                                     ej = elm$ej,
                                     tau_jB = elm$tau_jB,
                                     tauj_2 = elm$tauj_2,
                                     tau_j = elm$tau_j,
                                     Dmat = elm$Dmat,
                                     dvec = elm$dvec,
                                     Amat = elm$Amat,
                                     betabarj = nullvalue,
                                     betaj = 0, ## not needed anyway
                                     type = "typeI")
        NonstandardizedTests <- calcNonstandardizedTest(tauj_2 = elm$tauj_2,
                                                        tauj_inf = elm$tauj_inf,
                                                        alpha = alpha,
                                                        sigmasqbar = sigmasqbar["TypeI"])
        ## cat("\nnon: ", NonstandardizedTests)
        tbar <- min(NonstandardizedTests)
        ## cat("\ntbar: ", tbar)
        ## cat("\nbetahat: ", betahatj)

        res <- tbar - (betahatj - nullvalue)
        res
    }

findBernoulliCI <- function()
    {

    }

calcBernoulliCI <- function(nullvalue, )
