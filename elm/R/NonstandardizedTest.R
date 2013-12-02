
## Nonstandardized test statistics
calcNonstandardizedTest <- function(TAUj_2, TAUj_inf, alpha,
                                    sigmasqbar)
    {
        test <- rep(1, 5)
        names(test) <- c("Berry-Esseen", "Cantelli", "Bhattacharyya",
                          "Hoeffding", "Pinelis")

        ##
        ## Hoeffding
        ##
        test[4] <- (-0.5 * TAUj_2 * log(alpha))^(0.5)

        ##
        ## Pinelis
        ##
        test[5] <- sqrt(TAUj_2) * 0.5 * qnorm(1 - alpha/(factorial(5) * (exp(1)/5)^5 ))

        ##
        ## Cantelli
        ##
        test[2] <- ((sigmasqbar) * (1 - alpha)/alpha)^(0.5)

        ##
        ## Bhattacharyya
        ##
        if(test[2] * (test[2] - TAUj_inf)/sigmasqbar > 1)
            {
                test[3] <- sqrt(sigmasqbar * (1 + sqrt(3 * (1 - alpha)/alpha)))

                if(test[3]^2 * TAUj_inf < sigmasqbar * (TAUj_inf + 3 * test[3]))
                    {

                        r <- uniroot(function(x) ((3 * sigmasqbar - TAUj_inf^2) * sigmasqbar)/((3 * sigmasqbar - TAUj_inf^2) * (sigmasqbar + x^2) + (x^2 - x * TAUj_inf - sigmasqbar)^2) - alpha, c(0, test[3]))

                        test[3] <- r$root
                    }

            }
        else
            {
                test[3] <- Inf
            }

        ##
        ## Berry-Esseen
        ##
        wb1start <- c(0.1, 0.1)
        tbartemp <- test[4] ## start with Hoeffding bound
        test[1] <- Inf
        grid <- 0.00001
        lowerBE <- c(0.000001, 0.000001)

        BEtype1 <- optim(wb1start, minBerryEsseenbound, gr = NULL,
                         sigmasqbar = sigmasqbar,
                         TAUj_inf = TAUj_inf,
                         tbar = tbartemp, ## start with Hoeffding bound
                         method = "L-BFGS-B",
                         lower = lowerBE, upper = c(10, 10))
        ## control = list(fnscale = 1))
        BE <- (BEtype1$value/1000 < alpha)
        while(BE)
            {
                BEtype1 <- optim(wb1start, minBerryEsseenbound, gr = NULL,
                                 sigmasqbar = sigmasqbar,
                                 TAUj_inf = TAUj_inf,
                                 tbar = tbartemp, ## start with Hoeffding bound
                                 method = "L-BFGS-B",
                                 lower = lowerBE, upper = c(10, 10))

                if(BEtype1$value/1000 > alpha)
                    {
                        (BE <- FALSE)
                    }
                else
                    {
                        test[1] <- tbartemp
                        tbartemp <- tbartemp - grid
                        wb1start <- BEtype1$par
                    }
            }
        return(test)
    }

minBerryEsseenbound <- function(wb, sigmasqbar, TAUj_inf, tbar)
    {
        ## if(sigmasqbar < 2 * (wb[1])^2)
        ##     {
        ##         (R <- TAUj_inf * sigmasqbar/(sigmasqbar + (wb[1])^2)^(3/2))
        ##     }
        ## else
        ##     {
        ##         (R <- 2 * TAUj_inf/((27^0.5) * wb[1]))
        ##     }
        R <- ifelse(sigmasqbar < 2 * wb[1]^2,
                    TAUj_inf * sigmasqbar/(sigmasqbar + (wb[1])^2)^(3/2),
                    2 * TAUj_inf/((27^0.5) * wb[1]))
        res <- 1000 * (1 - pnorm((tbar - wb[2])/(sigmasqbar + wb[1]^2)^0.5) + 0.56 * R)/pnorm(wb[2]/wb[1])
        res
    }

calcBhattacharyya <- function(Bhattbar, sigmasqbar, TAUj_inf)
    {
        if(((Bhattbar^2)/sigmasqbar) - Bhattbar * TAUj_inf/sigmasqbar - 1 <= 0)
            {
                gOLStemp <- 1
            } else {
                if(sigmasqbar < (Bhattbar^2) * TAUj_inf/(TAUj_inf + 3 * Bhattbar))
                    {
                        gOLStemp <- (3 * sigmasqbar^2)/(4 * (sigmasqbar^2) - 2 * sigmasqbar * Bhattbar^2 + Bhattbar^4)
                    } else {
                        gOLStemp <- ((3 * sigmasqbar - TAUj_inf^2) * sigmasqbar)/((3 * sigmasqbar - TAUj_inf^2) * (sigmasqbar + Bhattbar^2) + (Bhattbar^2 - Bhattbar * TAUj_inf - sigmasqbar)^2)
                    }
            }
        gOLStemp
    }



## TypeII of the Nonstandardized test
calcTypeIINonstandardized <- function(wb1start, lowerBE,
                                      sigmasqbar,
                                      betaj,
                                      betabarj,
                                      tbarmin,
                                      TAUj_2,
                                      TAUj_inf)
    ## should work for OLS and MM (I see no differences)
    {
        typeII <- rep(1, 5)
        names(typeII) <- c("Berry-Esseen", "Cantelli", "Bhattacharyya",
                           "Hoeffding", "Pinelis")

        ## Berry-Esseen
        typeII[1] <- optim(wb1start, minBerryEsseenbound, gr = NULL,
                           sigmasqbar = sigmasqbar,
                           TAUj_inf = TAUj_inf,
                           tbar = betaj - betabarj - tbarmin,
                           method = "L-BFGS-B",
                           lower = lowerBE, upper = c(10,10),
                           control = list(fnscale = 1))$value/1000

        ## Cantelli
        typeII[2] <- (sigmasqbar)/(sigmasqbar + (betaj - betabarj - tbarmin)^2)

        ## Bhattacharyya
        typeII[3] <- calcBhattacharyya(betaj - betabarj - tbarmin,
                                       sigmasqbar = sigmasqbar,
                                       TAUj_inf = TAUj_inf)

        ## Hoeffding
        typeII[4] <- exp(-2 * (betaj - betabarj - tbarmin)^2/TAUj_2)

        ## Pinelis
        typeII[5] <- factorial(5) * (exp(1)/5)^5 * (1 - pnorm(2 * (betaj - betabarj - tbarmin)/sqrt(TAUj_2)))

        typeII[typeII > 1] <- 1

        return(typeII)
    }

## findOptimalBetajTypeIINonstandardized <- function(betaj)
##     {
##         typeII <- calcTypeIINonstandardized(wb1start, lowerBE,
##                                         sigmasqbar,
##                                         betaj = betaj, betabarj,
##                                         tbarmin, TAUj_2)
##         return(typeII[which.min(typeII)] - 0.5)
##     }
