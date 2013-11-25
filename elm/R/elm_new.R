
elm <- function(Y, X, some parameter)
    {
        for(coefficients)
            {
                res <- elmTest("OLS")
                res <- elmTest("mm")

                ## select estimator
            }

    }

elmTest <- function(estimator)
    {
        ##
        ## calc Tests
        ##
        ## Nonstandardized:
        ## Cantelli
        ## Hoeffding
        ## Bhattarayya
        ## Berry-Esseen
        ## Pintelis
        calcNonstandardized
        ##
        ## Bernoulli:
        calcBernoulli

        ## select Test
        return(list(variablename,
                    coefficient,
                    testValue,
                    testName,
                    testEstimator))
    }
