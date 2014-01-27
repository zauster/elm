
## steps:

## take elm-object
## see what test
## if nonstandardized -> nonstandardized
## if bernoulli -> bernoulli
## iterate over betabarj to find the bound

## TODO
## set bounds correctly: if null rejected, right, else left
## silence when CI
## let user select which test to use, if not, use from nullvalue = 0

elmCI <- function(elm, alpha = 0.05, coefs = 2)
    {
        for(lst in elm$coefTests)
            {
                switch(lst$chosenTest$'chosen Test',
                       'Nonstandardized (OLS)' = {
                           cat("1\n")
                           res <- findNonstandardizedCI(upperbetabound = lst$upperbetabound,
                                                        ## nullvalue =
                                                        ## lst$nullvalue,
                                                        nullvalue = 0,
                                                        betahatj = lst$betahatj[1],
                                                        elm = lst$model,
                                                        alpha = alpha)$root
                           cat("\nCI bound: ", res, "\n")
                       },
                       'Nonstandardized (MM)' = {
                           cat("2")
                       },
                       '       Bernoulli (OLS)' = {
                           cat("3")
                       },
                       '       Bernoulli (MM)' = {
                           cat("4")
                           res <- findBernoulliCI(upperbetabound = lst$upperbetabound,
                                                  nullvalue = 0, #lst$nullvalue,
                                                  betahatj = lst$betahatj[2],
                                                  elm = lst$model,
                                                  alpha = alpha)$root
                           cat("\nCI bound: ", res, "\n")
                       })
            }
    }

findBernoulliCI <- function(upperbetabound, nullvalue, betahatj,
                            elm, iter = 0,
                            step = 0.05 * upperbetabound,
                            alternative, alpha)
    {
        betainterval <- c(-upperbetabound, nullvalue)
        cat("\nbetainterval 1: ", betainterval)

        res <- try(uniroot(calcBernoulliCI,
                           interval = betainterval,
                           betahatj = betahatj,
                           upperbetabound = upperbetabound,
                           elm = elm,
                           alternative = alternative,
                           alpha = alpha),
                   silent = TRUE)

        if(is.error(res) == TRUE)
            {
                if(grepl("uniroot", res[1]) == TRUE)
                    {
                        iter <- iter + 1
                        ifelse(iter > 40, stop("ITERITERITERITER"),
                               upperbetabound)
                        res <- findBernoulliCI(upperbetabound = upperbetabound + step,
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

calcBernoulliCI <- function(nullvalue, betahatj,
                            upperbetabound, elm,
                            alternative, alpha = 0.05,
                            iterations = 1000,
                            zero = 10^-6,
                            max.iter = 5000,
                            lambda = 1,
                            root = TRUE)
    {
        betainterval <- c(nullvalue, upperbetabound)
        cat("\nbetainterval 2: ", betainterval)
        optbetaj <- uniroot(calcTypeIIBernoulli,
                            interval = betainterval,
                            betabarj = nullvalue,
                            alpha = alpha,
                            ds = elm$ds, b = elm$b,
                            XROWS = elm$XROWS,
                            root = root)$root
        BernoulliTypeII <- calcTypeIIBernoulli(betaj = optbetaj,
                                               betabarj = nullvalue,
                                               alpha = alpha,
                                               ds = elm$ds, b = elm$b,
                                               XROWS = elm$XROWS,
                                               root = FALSE)

        res <- calcBernoulliTest(Y = elm$Y, XROWS = elm$XROWS,
                                 tauj_inf = elm$tauj_inf,
                                 tau_j = elm$tau_j,
                                 ww = elm$ww,
                                 Bernoulliparameter = BernoulliTypeII,
                                 d = elm$d, ds = elm$ds,
                                 b = elm$b, a = 0, # betaj,
                                 betabarj = nullvalue,
                                 alternative = alternative,
                                 iterations = iterations,
                                 zero = 10^-6,
                                 max.iter = 5000,
                                 root = root,
                                 lambda = 1)
        res
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
