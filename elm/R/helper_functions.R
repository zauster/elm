
## TypeII of the Nonstandardized test
TypeIINonstandardized <- function(wb1start = c(0.1, 0.1),
                                  ## hOLS = hOLS,
                                  lowerBE = rep(10^-6, 2),
                                  sigmasqbar_betajOLSvalue = 2533,
                                  betaj11 = -270,
                                  betabarj0 = -100,
                                  tbarmin = 134.7,
                                  ## gOLS = gOLS,
                                  TAUj_2 = 12858)
    ## should work for OLS and MM (I see no differences)
    {
        typeII <- rep(1, 5)
        names(typeII) <- c("Berry-Esseen", "Cantelli", "Bhattacharyya",
                           "Hoeffding", "Pinelis")

        ## Berry-Esseen
        typeII[1] <- optim(wb1start, hOLS, gr = NULL,
                           method = "L-BFGS-B",
                           lower = lowerBE, upper = c(10,10),
                           control = list(fnscale = 1))$value/1000

        ## Cantelli
        typeII[2] <- (sigmasqbar_betajOLSvalue)/(sigmasqbar_betajOLSvalue + (betaj11 - betabarj0 - tbarmin)^2)

        ## Bhattacharyya
        typeII[3] <- gOLS(betaj11 - betabarj0 - tbarmin)

        ## Hoeffding
        typeII[4] <- exp(-2 * (betaj11 - betabarj0 - tbarmin)^2/TAUj_2)

        ## Pinelis
        typeII[5] <- factorial(5) * (exp(1)/5)^5 * (1 - pnorm(2 * (betaj11 - betabarj0 - tbarmin)/sqrt(TAUj_2)))

        typeII[typeII > 1] <- 1

        return(typeII)
    }

findOptimalBetajTypeIINonstandardized <- function(betaj)
    {
        typeII <- TypeIINonstandardized(wb1start, lowerBE,
                                        sigmasqbar_betajOLSvalue,
                                        betaj11 = betaj, betabarj0,
                                        tbarmin, TAUj_2)
        return(typeII[which.min(typeII)] - 0.5)
    }
