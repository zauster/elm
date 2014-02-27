
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
## make BernoulliTest output a named vector instead of list [LOW]

elmCI <- function(elm, alpha = 0.025, coefs = 2,
                  dispWarnings = FALSE,
                  useTest = NULL)
    {
        alternative <- elm$parameter$alternative
        ## if(alternative == "less")
        ##     {
        ##         stop("Not implemented yet. Please test with 'alternative = \"greater\"' and try again.")
        ##     }

        CIs <- list()
        for(lst in elm$coefTests)
            {
                ## print(lst)
                rejection <- lst$chosenTest$Rejection

                ## if useTest is not NULL, set the test to find the CI manually
                if(!is.null(useTest))
                    {
                        lst$chosenTest$'chosen Test' <- switch(useTest,
                                                               {"Nonstandardized (OLS)"},
                                                               {"      Bernoulli (OLS)"},
                                                               {"Nonstandardized (MM)"},
                                                               {"      Bernoulli (MM)"})
                    }

                switch(lst$chosenTest$'chosen Test',
                       'Nonstandardized (OLS)' = {
                           ## cat("1\n")
                           res <- findNonstandardizedCI(upperbetabound = lst$upperbetabound,
                                                        nullvalue = 0,
                                                        rejection = rejection,
                                                        alternative = alternative,
                                                        betahatj = lst$betahatj[1],
                                                        elm = lst$modelOLS,
                                                        alpha = alpha)$root
                           ## cat("\nres: ", res, "\n")
                           confint <- c(res, lst$betahatj["OLS"] + (lst$betahatj["OLS"] - res))
                           ## cat("\nCI bound:\n", confint)
                           ## return(confint)
                       },
                       '      Bernoulli (OLS)' = {
                           ## cat("2\n")
                           res <- findBernoulliCI(upperbetabound = lst$upperbetabound,
                                                  nullvalue = 0, #lst$nullvalue,
                                                  betahatj = lst$betahatj[1],
                                                  elm = lst$modelOLS,
                                                  alpha = alpha,
                                                  dispWarnings = dispWarnings)$root
                           ## res <- ifelse(alternative == "greater",
                           ##               res, -res)
                           confint <- c(-res, res)
                       },
                       'Nonstandardized (MM)' = {
                           ## cat("3\n")
                           res <- findNonstandardizedCI(upperbetabound = lst$upperbetabound,
                                                        nullvalue = 0,
                                                        rejection = rejection,
                                                        betahatj = lst$betahatj[2],
                                                        elm = lst$modelMM,
                                                        alpha = alpha)$root
                           ## cat("\nres: ", res, "\n")
                           confint <- c(res, lst$betahatj["MM"] + (lst$betahatj["MM"] - res))
                           ## cat("\nCI bound:\n", confint)
                       },
                       '      Bernoulli (MM)' = {
                           ## cat("4\n")   #
                           res <- findBernoulliCI(upperbetabound = lst$upperbetabound,
                                                  nullvalue = 0, #lst$nullvalue,
                                                  betahatj = lst$betahatj[2],
                                                  elm = lst$modelMM,
                                                  alpha = alpha,
                                                  dispWarnings = dispWarnings)$root
                           ## cat("res: ", res, "\n")
                           confint <- ifelse(alternative == "greater",
                                         res, -res)
                           ## cat("\nCI bound:\n")
                           ## return(res)
                       })

                CIs[[length(CIs) + 1L]] <- list(coefname = lst$coefname,
                                                confint = confint,
                                                chosenTest = lst$chosenTest$'chosen Test')
            }

        method <- "Confidence Intervals for exact linear models"
        structure(list(CI = CIs,
                       method = method,
                       yname = elm$yname,
                       xname = elm$xname,
                       n = elm$parameter$n,
                       m = elm$parameter$m,
                       alpha = alpha),
                  class = "elmCI")
    }


findBernoulliCI <- function(upperbetabound, nullvalue, betahatj,
                            elm, iter = 0,
                            step = 0.05 * upperbetabound,
                            alternative, alpha,
                            dispWarnings = FALSE)
    {
        betainterval <- c(-upperbetabound, nullvalue)
        ## cat("\nbetainterval 1: ", betainterval)

        res <- try(uniroot(calcBernoulliCI,
                           interval = betainterval,
                           betahatj = betahatj,
                           upperbetabound = upperbetabound,
                           elm = elm,
                           alternative = alternative,
                           alpha = alpha,
                           dispWarnings = dispWarnings),
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
                                               elm = elm,
                                               iter = iter,
                                               alternative = alternative,
                                               alpha = alpha,
                                               dispWarnings = dispWarnings)
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
                            root = TRUE,
                            dispWarnings = FALSE)
    {
        betainterval <- c(nullvalue, upperbetabound)
        ## cat("\nbetainterval 2: ", betainterval)
        optbetaj <- uniroot(calcTypeIIBernoulli,
                            interval = betainterval,
                            betabarj = nullvalue,
                            alpha = alpha,
                            ds = elm$ds, b = elm$b,
                            XROWS = elm$XROWS,
                            root = root)$root
        ## cat("\noptbetaj: ", optbetaj)
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
                                 lambda = 1,
                                 dispWarnings = dispWarnings)
        res
    }


findNonstandardizedCI <- function(upperbetabound, nullvalue, rejection, alternative, betahatj,
                                  elm, alpha, iter = 0, step = 0.05 * upperbetabound)
    {
        ## if(rejection == TRUE | (alternative == "less"))
        ## if(alternative == "greater")
        ##     {
        ##         betainterval <- c(-nullvalue, upperbetabound)
        ##     }
        ## else
        ##     {
                betainterval <- c(-upperbetabound, upperbetabound)
        ##     }
        ## cat("\nbetainterval: ", betainterval)

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
                                                     rejection = rejection,
                                                     alternative = alternative,
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
        ## res <- tbar - (nullvalue - betahatj)
        res
    }
