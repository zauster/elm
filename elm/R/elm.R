

library(Rlab)
data(golf)
Y <- matrix(golf$score)
X <- cbind(1, as.matrix(golf[, -1]))
## printdetails <- T
## IE <- "<="
## alpha <- 0.05
## j <- 1
## betabarj <- 100
## betaj1 <- betabarj + 170
## w1Y <- 60
## w2Y <- 80
## lambda <- 1
## lambdamm <- 1
## monte <- 1000
## qq <- 0.0001 ##OLS (default=0.0001)
## qqmm <- 0.0001 ##MM (default=0.0001)
elm(Y, X, 60, 80, j = 1, betabarj = 100, betaj1 = 284)

## step example
n <- 400
h <- 0.5
YY <- sample(c(0, 1), size = n, replace = TRUE)
XX <- cbind(1, runif(n = n) < h)
## alpha <- 0.05
## Given Y=X*beta+error where there are no assumptions imposed on the
## errors, it tests the one sided
## hypothesis H0: betaj<=betabarj against H1: betaj>betabarj where j
## is index of coefficient.
## It also tests H0: betaj>=betabarj against H1: betaj<betabarj.

## IE <- "<="
## w1Y <- 0
## w2Y <- 1
## j <- 2
## betabarj = 0
## betaj1 = .13
## lambdamm <- 1
## monte <- 1000
## qq <- 0.0001 ##OLS (default=0.0001)
## qqmm <- 0.0001 ##MM (default=0.0001)
elm(YY, XX, 0, 1, j = 2, betabarj = 0, betaj1 = .13)

elm <- function(Y, X, w1Y, w2Y, IE = "<=", alpha = 0.05, j = 2,
                betabarj = 0, betaj1 = betabarj + 1.1,
                lambda = 1, lambdamm = 1,
                monte = 1000, qq = 0.0001, qqmm = 0.0001,
                printdetails = TRUE)
    {

        ## require(Rlab) ## not needed anymore, substituted rbinom for
        ## rbern
        require(Rglpk)
        require(quadprog)
        YNAME <- deparse(substitute(Y))
        XNAME <- deparse(substitute(X))

        if(min(Y) < w1Y | max(Y) > w2Y)
            {
                stop("your values of Y are not within the range [w1Y, w2Y]!")
            }

        if(ncol(X)==1)
            {
                stop("you need at least two covariates, add a constant")
            }

        ## check if X has full rank
        if(det(t(X) %*% X) == 0)
            {
                stop("columns are not independent, need to have full rank, eliminate some")
            }

        ## adjustments to deal with both inequalities
        if(IE == "<=")
            {
                w11Y <- w1Y
                w22Y <- w2Y
                IEc <- 1
                IEa <- ">="
                betabarj0 <- betabarj
                betaj11 <- betaj1
                if(betaj1 <= betabarj)
                    {
                        stop("you need to choose betaj1 > betabarj")
                    }
            } else if(IE == ">=") {
                w11Y <- -w2Y
                w22Y <- -w1Y
                Y <- -1 * Y
                IEc <- -1
                IEa <- "<="
                betabarj0 <- -1 * betabarj
                betaj11 <- betabarj0 + betabarj - betaj1
                if(betaj1 >= betabarj)
                    {
                        stop("you need to choose betaj1 < betabarj")
                    }
            } else {
                stop("you need to choose IE either as <= or >=")
            }


        XROWS <- nrow(X)
        XCOLS <- ncol(X)
        YROWS <- XROWS

        Y <- Y /(w22Y - w11Y) ## puts Y in [ww,ww + 1] where ww = w11Y/(w22Y-w11Y)
        X <- X /(w22Y - w11Y) ## rescales X so that beta remains unchanged
        ww <- w11Y /(w22Y - w11Y)

        Betahat <- solve((crossprod(X))) %*% crossprod(X, Y) # TODO:
                                        # numerically stable?
        Betahatj <- Betahat[j]
        cj <- j

        ## Summaryinitial <- c("alpha ", "coefficient index (j in the paper)",
        ##                     "n","m","betahatj", alpha, j, XROWS,
        ##                     XCOLS, round(Betahatj,4))
        ## dim(Summaryinitial) <- c(5,2)

        ## mydatainitial <- as.data.frame(Summaryinitial)
        ## names(mydatainitial) <- c(" ", "Data input")
        ## if(printdetails)
        ##     {
                ## print(mydatainitial)
        ##         cat("stop 1")
        ##     }

        cat("\n NON-STANDARDIZED TEST \n")

        ## OLS estimate
        Tau <- X %*% solve(crossprod(X))
        Tau_j <- Tau[,j]
        TAUj_2 <- sum((Tau_j)^2)
        TAUj_inf <- max(abs(Tau_j))

        ## Type 1 OLS
        ej <- rep(0, times = XCOLS)
        ej[j] <- 1

        ## Finding upper bound on standard deviation of estimator:
        ## adjusting tau_j st Dmat positive definite,
        ## rescale and then eliminate 0 entries
        ## increases bound so test remains exact

        Tau_jB <- rep(0, times = XROWS)
        i <- 1
        TAUj_inf_min001 <- min(0.01, TAUj_inf) ## less computations
        TAUj_inf_min1 <- min(1, TAUj_inf)
        while (i <= XROWS)
            {
                if(Tau_j[i] >= 0)
                    {
                        Tau_jB[i] <- max(Tau_j[i]/TAUj_inf_min001,
                                         qq * TAUj_inf_min1)
                    } else {
                        Tau_jB[i] <- min(Tau_j[i]/TAUj_inf_min001,
                                         -qq * TAUj_inf_min1)
                    }
                i <- i + 1
            }

        DmatOLS <- 2 * t(X) %*% ((Tau_jB^2) * X)
        dvecOLS <- (1 + 2 * ww) * colMeans((Tau_jB^2) * X) * XROWS
        Amat <- t(matrix(rbind(-t(ej), X, -X), ncol = XCOLS))
        bvec0 <- as.vector(c(-betabarj0, rep(ww, times = XROWS),
                             rep(-(ww + 1), times = XROWS)))

        zero <- 10^-12
        if(det(DmatOLS) <= zero)
            {
                print(paste("Dmat not positive definite for OLS, Det=",
                            det(DmatOLS)))
                sigmasqbar_beta0jOLSvalue <- TAUj_2/4 - (1/XROWS) * (min(0, betabarj0-(1/2) * sum(Tau_j)))^2

            } else {
                sigmasqbar_beta0jOLSQP <- solve.QP(Dmat = DmatOLS,
                                                   dvec = dvecOLS, Amat,
                                                   bvec = bvec0, meq = 0,
                                                   factorized = FALSE)
                zsol0OLS <- sigmasqbar_beta0jOLSQP$solution
                sigmasqbar_beta0jOLSvalue <- sum((Tau_j)^2 * ((X %*% zsol0OLS - ww) * (1 + ww - X %*% zsol0OLS)))
            }

        ##
        ##
        ##
        ##
        ##
        ##
        ## Hoeffding, Pinelis and Cantelli OLS
        ##
        ##
        ##
        ##
        ##
        ##


        Hoeftbar <- (-0.5 * TAUj_2 * log(alpha))^(0.5)
        Pintbar <- sqrt(TAUj_2) * 0.5 * qnorm(1 - alpha/(factorial(5) * (exp(1)/5)^5 ) )
        Cantellitbar <- ((sigmasqbar_beta0jOLSvalue) * (1 - alpha)/alpha)^(0.5)

        ##
        ##
        ##
        ##
        ## Bhattacharyya OLS
        ##
        ##
        ##
        ##
        ##

        if(Cantellitbar * (Cantellitbar - TAUj_inf)/sigmasqbar_beta0jOLSvalue > 1)
            {
                Bhattbar <- sqrt(sigmasqbar_beta0jOLSvalue * (1 + sqrt(3 * (1 - alpha) / alpha)))

                if(Bhattbar^2 * TAUj_inf < sigmasqbar_beta0jOLSvalue * ( TAUj_inf + 3 * Bhattbar))
                    {

                        r <- uniroot(function(x) ((3 * sigmasqbar_beta0jOLSvalue - TAUj_inf^2) * sigmasqbar_beta0jOLSvalue)/((3 * sigmasqbar_beta0jOLSvalue - TAUj_inf^2) * (sigmasqbar_beta0jOLSvalue + x^2) + (x^2 - x * TAUj_inf - sigmasqbar_beta0jOLSvalue)^2) - alpha, c(0, Bhattbar))

                        Bhattbar <- r$root
                    }

            } else {
                Bhattbar <- Inf
            }

        ##
        ##
        ##
        ##
        ##
        ## Berry-Esseen OLS
        ##
        ##
        ##
        ##

        h <- function(wb1)
            {
                if(sigmasqbar_beta0jOLSvalue < 2 * (wb1[1])^2)
                    {
                        (R <- TAUj_inf * sigmasqbar_beta0jOLSvalue/(sigmasqbar_beta0jOLSvalue + (wb1[1])^2)^(3/2))
                    } else {
                        (R <- 2 * TAUj_inf/((27^0.5) * wb1[1]))
                    }
                htemp <- 1000 * (1 - pnorm((tbartemp - wb1[2])/(sigmasqbar_beta0jOLSvalue + wb1[1]^2)^0.5) + 0.56 * R)/pnorm(wb1[2]/wb1[1])
                htemp
            }

        wb1start <- c(0.1, 0.1)
        tbartemp <- Hoeftbar
        BEtbar <- Inf
        grid <- 0.00001
        lowerBE <- c(0.000001, 0.000001)

        BEtype1 <- optim(wb1start, h, gr = NULL,
                         method = "L-BFGS-B",
                         lower = lowerBE, upper = c(10, 10),
                         control = list(fnscale = 1))

        BE <- (BEtype1$value/1000 < alpha)
        while(BE)
            {
                BEtype1 <- optim(wb1start, h, gr = NULL,
                                 method = "L-BFGS-B",
                                 lower = lowerBE, upper = c(10,10),
                                 control = list(fnscale = 1))

                if(BEtype1$value/1000 > alpha)
                    {
                        (BE <- FALSE)
                    } else {
                        BEtbar <- tbartemp
                        ## print(BEtbar)
                        tbartemp <- tbartemp - grid
                        wb1start <- BEtype1$par
                    }
            }

        ##
        ##
        ##
        ##
        ##
        ## Final results for Non-Standardized test OLS cutoffs
        ##
        ##
        ##
        ##
        ##

        tbarmin <- min(c(Cantellitbar, Bhattbar, Hoeftbar, Pintbar, BEtbar))

        Summary1 <- c("betabarj", "Betahatj", "Cantellitbar", "Bhattbar",
                      "Hoeftbar", "Pintbar", "BEtbar",
                      round(betabarj, digits = 5),
                      round(IEc * Betahatj, digits = 5),
                      round(Cantellitbar, digits = 5),
                      round(Bhattbar, digits = 5),
                      round(Hoeftbar, digits = 5),
                      round(Pintbar, digits = 5),
                      round(BEtbar, digits = 5))
        dim(Summary1) <- c(7,2)

        mydata <- as.data.frame(Summary1)
        names(mydata) <- c("statistic name (OLS version)", "value")
        if(printdetails)
            {
                cat("stop 1\n")
                print(mydata)
            }
### xxxxxx
        if(Betahatj - betabarj0 >= tbarmin)
            {
                if(printdetails)
                    {
                        cat("Nonstandardized test >>> REJECTS H_0\n\n")
                    }
                RejectNonstandardizedOLS <- "YES"
            } else {
                if(printdetails)
                    {
                        cat("Non-standardized test does not reject H_0 \n\n")
                    }
                RejectNonstandardizedOLS <- "NO"
            }


        ## Type 2 OLS
        if(det(DmatOLS) <= zero)
            {
                sigmasqbar_betajOLSvalue <- TAUj_2/4 - (1/XROWS) * (betaj11-(1/2) * sum(Tau_j))^2

            } else {
                bvec1 <- as.vector(c(-betaj11,
                                     rep(ww, times = XROWS),
                                     rep(-(ww + 1), times = XROWS)))
                c2 <- 1/max(DmatOLS) ## NEW

                sigmasqbar_betajOLSQP <- solve.QP(Dmat = c2 * DmatOLS,
                                                  dvec = c2 * dvecOLS,
                                                  Amat = Amat, bvec = bvec1,
                                                  meq = 1,
                                                  factorized = FALSE)

                ## if there is an error message
                ## Fehler in solve.QP(Dmat = c2 * DmatOLS, dvec = c2 * dvecOLS, Amat = Amat,  :
                ##  constraints are inconsistent, no solution!
                ## then this often results to numerical problems when solving (2) in
                ## (Gossner and Schlag, 2013)
                ## but can also come if betaj1 is not close enough to 0

                zsolOLS <- sigmasqbar_betajOLSQP$solution

                sigmasqbar_betajOLSvalue <- sum((Tau_j)^2 * ((X %*% zsolOLS-ww) * (1 + ww-X %*% zsolOLS)))
            }

        ## gOLS <- function(Bhattbar)
        ##     {

        ##         if(((Bhattbar^2)/sigmasqbar_betajOLSvalue) - Bhattbar * TAUj_inf/sigmasqbar_betajOLSvalue - 1 <= 0) gOLStemp <- 1
        ##         else
        ##         if(sigmasqbar_betajOLSvalue < (Bhattbar^2) * TAUj_inf/(TAUj_inf + 3 * Bhattbar))
        ##             gOLStemp <- (3 * sigmasqbar_betajOLSvalue^2)/(4 * (sigmasqbar_betajOLSvalue^2)-2 * sigmasqbar_betajOLSvalue * Bhattbar^2 + Bhattbar^4)
        ##         else gOLStemp <- ((3 * sigmasqbar_betajOLSvalue - TAUj_inf^2) * sigmasqbar_betajOLSvalue)/((3 * sigmasqbar_betajOLSvalue - TAUj_inf^2) * (sigmasqbar_betajOLSvalue + Bhattbar^2) + (Bhattbar^2 - Bhattbar * TAUj_inf - sigmasqbar_beta0jOLSvalue)^2)
        ##         gOLStemp
        ##     }
        gOLS <- function(Bhattbar)
            {
                if(((Bhattbar^2)/sigmasqbar_betajOLSvalue) - Bhattbar * TAUj_inf/sigmasqbar_betajOLSvalue - 1 <= 0)
                    {
                        gOLStemp <- 1
                    } else {
                        if(sigmasqbar_betajOLSvalue < (Bhattbar^2) * TAUj_inf/(TAUj_inf + 3 * Bhattbar))
                            {
                                gOLStemp <- (3 * sigmasqbar_betajOLSvalue^2)/(4 * (sigmasqbar_betajOLSvalue^2) - 2 * sigmasqbar_betajOLSvalue * Bhattbar^2 + Bhattbar^4)
                            } else {
                                gOLStemp <- ((3 * sigmasqbar_betajOLSvalue - TAUj_inf^2) * sigmasqbar_betajOLSvalue)/((3 * sigmasqbar_betajOLSvalue - TAUj_inf^2) * (sigmasqbar_betajOLSvalue + Bhattbar^2) + (Bhattbar^2 - Bhattbar * TAUj_inf - sigmasqbar_beta0jOLSvalue)^2)
                            }
                    }
                gOLStemp
            }


        hOLS <- function(wb1)
            {
                if(sigmasqbar_betajOLSvalue < 2 * (wb1[1])^2)
                    {
                        (R <- TAUj_inf * sigmasqbar_betajOLSvalue/(sigmasqbar_betajOLSvalue + (wb1[1])^2)^(3/2))
                    } else {
                        (R <- 2 * TAUj_inf/((27^0.5) * wb1[1]))
                    }
                htempOLS <- 1000 * (1 - pnorm(((betaj11 - betabarj0 - tbarmin) - wb1[2])/(sigmasqbar_betajOLSvalue + wb1[1]^2)^0.5) + 0.56 * R)/pnorm(wb1[2]/wb1[1])
                htempOLS
            }


        TYPEII_BEfuncOLS <- optim(wb1start, hOLS, gr = NULL,
                                  method = "L-BFGS-B",
                                  lower = lowerBE, upper = c(10,10),
                                  control = list(fnscale = 1))


        TYPEIIOLS_C <- min(1, (sigmasqbar_betajOLSvalue)/(sigmasqbar_betajOLSvalue + (betaj11 - betabarj0 - tbarmin)^2))

        TYPEIIOLS_Y <- min(1, gOLS(betaj11 - betabarj0 - tbarmin))

        TYPEIIOLS_BE <- min(1, TYPEII_BEfuncOLS$value/1000)

        TYPEIIOLS_H <- min(1, exp(-2 * (betaj11 - betabarj0 - tbarmin)^2/TAUj_2))

        TYPEIIOLS_P <- min(1, factorial(5) * (exp(1)/5)^5 * (1 - pnorm(2 * (betaj11 - betabarj0 - tbarmin)/sqrt(TAUj_2))))

        if(betaj11 <= betabarj0 + tbarmin)
            {
                TYPEIIOLS_C <- 1
                TYPEIIOLS_Y <- 1
                TYPEIIOLS_BE <- 1
                TYPEIIOLS_H <- 1
                TYPEIIOLS_P <- 1
            }



        if(tbarmin == Cantellitbar)
            {
                CutoffNonstandardizedOLS <- "C "
                ModelNonstandardizedOLS <- "Non-Standardized test"
            }

        if(min(TYPEIIOLS_C, TYPEIIOLS_Y, TYPEIIOLS_BE,
               TYPEIIOLS_H, TYPEIIOLS_P) == TYPEIIOLS_C)
            {
                TYPEIINonstandardizedOLS <- TYPEIIOLS_C
                TYPEIINonstandardizedOLSname <- "C "
            }

        if(tbarmin == Bhattbar)
            {
                CutoffNonstandardizedOLS <- "Bh "
                ModelNonstandardizedOLS <- "Non-Standardized test"
            }

        if(min(TYPEIIOLS_C, TYPEIIOLS_Y, TYPEIIOLS_BE,
               TYPEIIOLS_H, TYPEIIOLS_P) == TYPEIIOLS_Y)
            {
                TYPEIINonstandardizedOLS <- TYPEIIOLS_Y
                TYPEIINonstandardizedOLSname <- "Bh "
            }

        if(tbarmin == Hoeftbar )
            {
                CutoffNonstandardizedOLS <- "H "
                ModelNonstandardizedOLS <- "Non-Standardized test"
            }

        if(min(TYPEIIOLS_C, TYPEIIOLS_Y, TYPEIIOLS_BE,
               TYPEIIOLS_H, TYPEIIOLS_P) == TYPEIIOLS_H)
            {
                TYPEIINonstandardizedOLS <- TYPEIIOLS_H
                TYPEIINonstandardizedOLSname <- "H "
            }

        if(tbarmin == Pintbar )
            {
                CutoffNonstandardizedOLS <- "P "
                ModelNonstandardizedOLS <- "Non-Standardized test"
            }

        if(min(TYPEIIOLS_C, TYPEIIOLS_Y, TYPEIIOLS_BE,
               TYPEIIOLS_H, TYPEIIOLS_P) == TYPEIIOLS_P)
            {
                TYPEIINonstandardizedOLS <- TYPEIIOLS_P
                TYPEIINonstandardizedOLSname <- "P "
            }

        if(tbarmin == BEtbar)
            {
                CutoffNonstandardizedOLS <- "BE "
                ModelNonstandardizedOLS <- "Non-Standardized test"
            }

        if(min(TYPEIIOLS_C, TYPEIIOLS_Y, TYPEIIOLS_BE,
               TYPEIIOLS_H, TYPEIIOLS_P) == TYPEIIOLS_BE)
            {
                TYPEIINonstandardizedOLS <- TYPEIIOLS_BE
                TYPEIINonstandardizedOLSname <- "BE "
            }

        Summary2OLS <- c("beta_j under alternative", "TYPE II bound Cantelli",
                         "TYPE II bound Bhattacharyya", "TYPE II bound Hoeffding",
                         "TYPE II bound Pinelis", "TYPE II bound Berry-Esseen",
                         round(betaj1, digits = 5),
                         round(TYPEIIOLS_C, digits = 5),
                         round(TYPEIIOLS_Y, digits = 5),
                         round(TYPEIIOLS_H, digits = 5),
                         round(TYPEIIOLS_P, digits = 5),
                         round(TYPEIIOLS_BE, digits = 5))
        dim(Summary2OLS) <- c(6,2)

        mydata2OLS <- as.data.frame(Summary2OLS)
        names(mydata2OLS) <- c("Type II related variables for OLS", "values")
        if(printdetails)
            {
                cat("stop 2\n")
                print(mydata2OLS)
            }


        ## -------------------- MM

        obj <- c(rep(0,XROWS),1)
        mat <- matrix(rbind(cbind(diag(XROWS), rep(-1,XROWS)),
                            cbind(diag(XROWS), rep(1,XROWS)),
                            cbind(t(X), rep(0,XCOLS)),
                            cbind(t(rep(0,XROWS)), 1)),
                      nrow = 2 * XROWS + XCOLS + 1)
        dir <- c(rep("<=", XROWS), rep(">=", XROWS), rep("==", XCOLS), ">=")
        rhs <- c(rep(0, 2 * XROWS), ej, 0)
        XROWSplusone <- XROWS + 1
        sol <- Rglpk_solve_LP(obj, mat, dir, rhs, types = NULL, max = FALSE,
                              bounds = list(lower = list(ind = 1:XROWSplusone,
                                                val = rep(-XROWS,XROWSplusone)),
                                  upper = list(ind = 1:XROWSplusone,
                                      val = rep(XROWS,XROWSplusone))),
                              verbose = FALSE)

        Taumm_j <- sol$solution
        Taumm_j <- Taumm_j[1:XROWS]
        Betahatmmj = crossprod(Taumm_j, Y)
        TAUjmm_2 = sum((Taumm_j)^2)
        TAUjmm_inf = max(abs(Taumm_j))


        ##Type 1 MM


        ## adjusting tau_j st Dmat positive definite,
        ## rescale and then eliminate 0 entries
        ## increases bound so test remains exact

        Taumm_jB <- rep(0,times = XROWS)
        i <- 1
        TAUjmm_2_min <- min(1, TAUjmm_2) ## less computations
        while (i <= XROWS)
            {
                if(Taumm_j[i] >= 0)
                    {
                        Taumm_jB[i] <- max(Taumm_j[i]/TAUjmm_2_min,
                                           qqmm * TAUjmm_2_min)
                    } else {
                        Taumm_jB[i] <- min(Taumm_j[i]/TAUjmm_2_min,
                                           -qqmm * TAUjmm_2_min)
                    }
                i <- i + 1
            }

        Dmatmm <- 2 * t(X) %*% ((Taumm_jB^2) * X)
        dvecmm <- (1 + 2 * ww) * colMeans((Taumm_jB^2) * X) * XROWS

        if(det(Dmatmm) <= zero)
            {
                print(paste("Dmat not positive definite for MM, Det=",
                            det(Dmatmm)))
                sigmasqbar_beta0jmmvalue <- TAUjmm_2/4 - (1/XROWS) * (min(0, betabarj0 - (1/2) * sum(Taumm_j)))^2
            } else {
                sigmasqbar_beta0jmmQP <- solve.QP(Dmat = Dmatmm, dvec = dvecmm,
                                                  Amat,
                                                  bvec = bvec0, meq = 0,
                                                  factorized = FALSE)
                zsol0mm <- sigmasqbar_beta0jmmQP$solution
                sigmasqbar_beta0jmmvalue <- sum((Taumm_j)^2 * ((X %*% zsol0mm - ww) * (1 + ww - X %*% zsol0mm)))
            }

        ##
        ##
        ##
        ##
        ##
        ##
        ## Hoeffding, Pinelis and Cantelli MM
        ##
        ##
        ##
        ##
        ##
        ##

        Hoeftbarmm <- (-0.5 * TAUjmm_2 * log(alpha))^0.5
        Pintbarmm <- sqrt(TAUjmm_2) * 0.5 * qnorm(1 - alpha/(factorial(5) * (exp(1)/5)^5))
        Cantellitbarmm <- ((sigmasqbar_beta0jmmvalue) * (1 - alpha)/alpha)^(0.5)

        ##
        ##
        ##
        ##
        ## Bhattacharyya MM
        ##
        ##
        ##
        ##
        ##

        if(Cantellitbarmm * (Cantellitbarmm - TAUjmm_inf)/sigmasqbar_beta0jmmvalue > 1)
            {
                Bhattbarmm <- sqrt(sigmasqbar_beta0jmmvalue * (1 + sqrt(3 * (1 - alpha) / alpha)))

                if(Bhattbarmm^2 * TAUjmm_inf < sigmasqbar_beta0jmmvalue * (TAUjmm_inf + 3  *  Bhattbarmm))
                    {
                        r <- uniroot(function(x) ((3 * sigmasqbar_beta0jmmvalue-TAUj_inf^2) * sigmasqbar_beta0jmmvalue)/((3 * sigmasqbar_beta0jmmvalue-TAUj_inf^2) * (sigmasqbar_beta0jmmvalue + x^2) + (x^2-x * TAUj_inf-sigmasqbar_beta0jmmvalue)^2) - alpha, c(0,Bhattbarmm))

                        Bhattbarmm <- r$root
                    }

            } else {
                Bhattbarmm <- Inf
            }

        ##
        ##
        ##
        ##
        ##
        ## Berry-Esseen MM
        ##
        ##
        ##
        ##

        hmm <- function(wb1)
            {
                if(sigmasqbar_beta0jmmvalue < 2 * (wb1[1])^2)
                    {
                        (R <- TAUj_inf * sigmasqbar_beta0jmmvalue/(sigmasqbar_beta0jmmvalue + (wb1[1])^2)^(3/2))
                    } else {
                        (R <- 2 * TAUj_inf/((27^0.5) * wb1[1]))
                    }
                hmmtemp = 1000 * (1-pnorm((tbartempmm-wb1[2])/(sigmasqbar_beta0jmmvalue + wb1[1]^2)^0.5) + 0.56 * R)/pnorm(wb1[2]/wb1[1])
                hmmtemp
            }

        h <- function(wb1)
            {
                if(sigmasqbar_beta0jOLSvalue < 2 * (wb1[1])^2)
                    {
                        (R <- TAUj_inf * sigmasqbar_beta0jOLSvalue/(sigmasqbar_beta0jOLSvalue + (wb1[1])^2)^(3/2))
                    } else {
                        (R <- 2 * TAUj_inf/((27^0.5) * wb1[1]))
                    }
                htemp <- 1000 * (1 - pnorm((tbartemp - wb1[2])/(sigmasqbar_beta0jOLSvalue + wb1[1]^2)^0.5) + 0.56 * R)/pnorm(wb1[2]/wb1[1])
                htemp
            }


        ## xxx: same as above?
        BEtype1 <- optim(wb1start, h, gr = NULL,
                         method = "L-BFGS-B",
                         lower = lowerBE, upper = c(10,10),
                         control = list(fnscale = 1))

        BE <- (BEtype1$value/1000 < alpha)

        while(BE)
            {
                BEtype1 <- optim(wb1start, h, gr = NULL,
                                 method = "L-BFGS-B",
                                 lower = lowerBE, upper = c(10,10),
                                 control = list(fnscale = 1))

                if(BEtype1$value/1000 > alpha)
                    {
                        (BE <- FALSE)
                    } else {
                        BEtbar <- tbartemp
                        ## print(BEtbar)
                        tbartemp <- tbartemp - grid
                        wb1start <- BEtype1$par
                    }
            }



        wb1mmstart <- c(0.1, 0.1)
        tbartempmm <- Hoeftbarmm
        BEtbarmm <- Inf

        BEtype1mm <- optim(wb1mmstart, hmm, gr = NULL,
                           method = "L-BFGS-B",
                           lower = lowerBE, upper = c(10,10),
                           control = list(fnscale = 1))

        BE <- (BEtype1mm$value/1000 < alpha)

        while(BE)
            {
                BEtype1mm <- optim(wb1mmstart, hmm, gr = NULL,
                                   method = "L-BFGS-B",
                                   lower = lowerBE, upper = c(10,10),
                                   control = list(fnscale = 1))

                if(BEtype1mm$value/1000 > alpha)
                    {
                        (BE <- FALSE)
                    } else {
                        BEtbarmm <- tbartempmm
                        ## print(BEtbarmm)
                        tbartempmm <- tbartempmm - grid
                        wb1mmstart <- BEtype1$par
                    }
            }

        ##
        ##
        ##
        ##
        ##
        ## Final results for non-standardized test cutoff MM
        ##
        ##
        ##
        ##
        ##

        tbarminmm <- min(c(Cantellitbarmm, Bhattbarmm, Hoeftbarmm,
                           Pintbarmm, BEtbarmm))

        if(printdetails)
            {
                cat("\nstop 3\n")
                cat("NON-STANDARDIZED TEST \n")
            }

        Summarymm1 <- c("betabarj", "Betahatj", "Cantellitbar", "Bhattbar",
                        "Hoeftbar", "Pintbar", "BEtbar",
                        round(betabarj, digits = 5),
                        round(IEc * Betahatmmj, digits = 5),
                        round(Cantellitbarmm, digits = 5),
                        round(Bhattbarmm, digits = 5),
                        round(Hoeftbarmm, digits = 5),
                        round(Pintbarmm, digits = 5),
                        round(BEtbarmm, digits = 5))

        dim(Summarymm1) <- c(7,2)

        mydata <- as.data.frame(Summarymm1)
        names(mydata) <- c("statistic name (Minimax version)", "value")
        if(printdetails)
            {
                print(mydata)
            }

        if(Betahatmmj - betabarj0 >= tbarminmm)
            {
                if(printdetails)
                    {
                        cat("Nonstandardized test >>> REJECTS H_0\n\n")
                    }
                RejectNonstandardizedmm <- "YES"
            } else {
                if(printdetails)
                    {
                        cat("Non-standardized test does not reject H_0\n\n")
                    }
                RejectNonstandardizedmm <- "NO"
            }

        ## Type 2 MM
        if(det(Dmatmm) <= zero)
            {
                sigmasqbar_betajmmvalue <- TAUjmm_2/4-(1/XROWS) * (betaj11-(1/2) * sum(Taumm_j))^2
            } else {
                bvec1 <- as.vector(c(-betaj11, rep(ww, times = XROWS),
                                     rep(-(ww + 1), times = XROWS)))

                sigmasqbar_betajmmQP <- solve.QP(Dmat = c2 * Dmatmm,
                                                 dvec = c2 * dvecmm,
                                                 Amat = Amat, bvec = bvec1,
                                                 meq = 1,
                                                 factorized = FALSE)
                zsolmm <- sigmasqbar_betajmmQP$solution
                sigmasqbar_betajmmvalue <- sum((Taumm_j)^2 * ((X %*% zsolmm-ww) * (1 + ww-X %*% zsolmm)))
            }

        gmm <- function(Bhattbar)
            {
                if(((Bhattbar^2)/sigmasqbar_betajmmvalue) - Bhattbar * TAUjmm_inf/sigmasqbar_betajmmvalue - 1 <=0)
                    {
                        gmmtemp <- 1
                    } else {
                        if(sigmasqbar_betajmmvalue < (Bhattbar^2) * TAUjmm_inf/(TAUjmm_inf + 2 * Bhattbar))
                            {
                                gmmtemp <- (3 * sigmasqbar_betajmmvalue^2)/(4 * (sigmasqbar_betajmmvalue^2) - 2 * sigmasqbar_betajmmvalue * Bhattbar^2 + Bhattbar^4)
                            } else {
                                gmmtemp <- ((3 * sigmasqbar_betajmmvalue - TAUjmm_inf^2) * sigmasqbar_betajmmvalue)/((3 * sigmasqbar_betajmmvalue - TAUjmm_inf^2) * (sigmasqbar_betajmmvalue + Bhattbar^2) + (Bhattbar^2 - Bhattbar * TAUjmm_inf - sigmasqbar_beta0jmmvalue)^2)
                            }
                    }
                gmmtemp
            }

        gOLS <- function(Bhattbar)
            {
                if(((Bhattbar^2)/sigmasqbar_betajOLSvalue) - Bhattbar * TAUj_inf/sigmasqbar_betajOLSvalue - 1 <=0)
                    {
                        gOLStemp <- 1
                    } else {
                        if(sigmasqbar_betajOLSvalue < (Bhattbar^2) * TAUj_inf/(TAUj_inf + 3 * Bhattbar))
                            {
                                gOLStemp <- (3 * sigmasqbar_betajOLSvalue^2)/(4 * (sigmasqbar_betajOLSvalue^2) - 2 * sigmasqbar_betajOLSvalue * Bhattbar^2 + Bhattbar^4)
                            } else {
                                gOLStemp <- ((3 * sigmasqbar_betajOLSvalue - TAUj_inf^2) * sigmasqbar_betajOLSvalue)/((3 * sigmasqbar_betajOLSvalue - TAUj_inf^2) * (sigmasqbar_betajOLSvalue + Bhattbar^2) + (Bhattbar^2 - Bhattbar * TAUj_inf - sigmasqbar_beta0jOLSvalue)^2)
                            }
                    }
                gOLStemp
            }

        hmm <- function(wb1)
            {
                if(sigmasqbar_betajmmvalue < 2 * (wb1[1])^2)
                    {
                        (R <- TAUj_inf * sigmasqbar_betajmmvalue/(sigmasqbar_betajmmvalue + (wb1[1])^2)^(3/2))
                    } else {
                        (R <- 2 * TAUj_inf/((27^0.5) * wb1[1]))
                    }
                htempmm <- 1000 * (1 - pnorm(((betaj11 - betabarj0 - tbarminmm) - wb1[2])/(sigmasqbar_betajmmvalue + wb1[1]^2)^0.5) + 0.56 * R)/pnorm(wb1[2]/wb1[1])
                htempmm
            }


        TYPEII_BEfuncmm <- optim(wb1mmstart, hmm, gr = NULL,
                                 method = "L-BFGS-B",
                                 lower = lowerBE, upper = c(10,10),
                                 control = list(fnscale = 1))


        TYPEIImm_C <- min(1, (sigmasqbar_betajmmvalue)/(sigmasqbar_betajmmvalue + (betaj11 - betabarj0 - tbarminmm)^2))

        TYPEIImm_Y <- min(1, gmm(betaj11 - betabarj0 - tbarminmm))

        TYPEIImm_BE <- min(1, TYPEII_BEfuncmm$value/1000)

        TYPEIImm_H <- min(1, exp(-2 * (betaj11 - betabarj0 - tbarminmm)^2/TAUjmm_2))

        TYPEIImm_P <- min(1, factorial(5) * (exp(1)/5)^5 * (1 - pnorm(2 * (betaj11 - betabarj0 - tbarminmm)/sqrt(TAUjmm_2))))

        if(betaj11 <= betabarj0 + tbarminmm)
            {
                TYPEIImm_C <- 1
                TYPEIImm_Y <- 1
                TYPEIImm_BE <- 1
                TYPEIImm_H <- 1
                TYPEIImm_P <- 1
            }

        if(tbarminmm == Cantellitbarmm)
            {
                CutoffNonstandardizedmm <- "C"
                ModelNonstandardizedmm <- "Non-Standardized test"
            }

        if(min(TYPEIImm_C, TYPEIImm_Y, TYPEIImm_BE,
               TYPEIImm_H, TYPEIImm_P) == TYPEIImm_C)
            {
                TYPEIINonstandardizedmm <- TYPEIImm_C
                TYPEIINonstandardizedmmname <- "C"
            }

        if(tbarminmm == Bhattbarmm )
            {
                CutoffNonstandardizedmm <- "Bh"
                ModelNonstandardizedmm <- "Non-Standardized test"
            }

        if(min(TYPEIImm_C, TYPEIImm_Y, TYPEIImm_BE,
               TYPEIImm_H, TYPEIImm_P) == TYPEIImm_Y)
            {
                TYPEIINonstandardizedmm <- TYPEIImm_Y
                TYPEIINonstandardizedmmname <- "Bh"
            }

        if(tbarminmm == Hoeftbarmm )
            {
                CutoffNonstandardizedmm <- "H"
                ModelNonstandardizedmm <- "Non-Standardized test"
            }

        if(min(TYPEIImm_C, TYPEIImm_Y, TYPEIImm_BE,
               TYPEIImm_H, TYPEIImm_P) == TYPEIImm_H)
            {
                TYPEIINonstandardizedmm <- TYPEIImm_H
                TYPEIINonstandardizedmmname <- "H"
            }

        if(tbarminmm == Pintbarmm )
            {
                CutoffNonstandardizedmm <- "P"
                ModelNonstandardizedmm <- "Non-Standardized test"
            }

        if(min(TYPEIImm_C, TYPEIImm_Y, TYPEIImm_BE,
               TYPEIImm_H, TYPEIImm_P) == TYPEIImm_P)
            {
                TYPEIINonstandardizedmm <- TYPEIImm_P
                TYPEIINonstandardizedmmname <- "P"
            }

        if(tbarminmm == BEtbarmm)
            {
                CutoffNonstandardizedmm <- "BE"
                ModelNonstandardizedmm <- "Non-Standardized test"
            }

        if(min(TYPEIImm_C, TYPEIImm_Y, TYPEIImm_BE,
               TYPEIImm_H, TYPEIImm_P) == TYPEIImm_BE)
            {
                TYPEIINonstandardizedmm <- TYPEIImm_BE
                TYPEIINonstandardizedmmname <- "BE"
            }

        Summary2mm <- c("beta_j under alternative", "TYPE II bound Cantelli",
                        "TYPE II bound Bhattacharyya", "TYPE II bound Hoeffding",
                        "TYPE II bound Pinelis", "TYPE II bound Berry-Esseen",
                        round(betaj1, digits = 5),
                        round(TYPEIImm_C, digits = 5),
                        round(TYPEIImm_Y, digits = 5),
                        round(TYPEIImm_H, digits = 5),
                        round(TYPEIImm_P, digits = 5),
                        round(TYPEIImm_BE, digits = 5))
        dim(Summary2mm) <- c(6,2)

        mydata2mm <- as.data.frame(Summary2mm)
        names(mydata2mm) <- c("Type II related variables for Minimax", "values")

        if(printdetails)
            {
                cat("\nstop 4\n")
                print(mydata2mm)
                cat("\n")
            }
        ##if(betaj11<=betabarj0 + tbarminmm){ cat("Doesn't apply (MM).\n\n") }


        ## ------------
        ##
        ##
        ##
        ##
        ##   Non-Randomised Bernoulli Test (OLS)
        ##
        ##
        ##
        ##
        ##
        ##

        a <- 0
        b <- XROWS * TAUj_inf
        i <- 1
        d <- rep(0, times = XROWS)
        while(i <= XROWS)
            {
                d[i] <- (1 - lambda)  *  max(-Tau_j[i] * ww,
                                             -Tau_j[i] * (ww + 1)) + lambda  *  (TAUj_inf - max(Tau_j[i] * ww, Tau_j[i] * (ww + 1)))
                i <- i + 1
            }
        ds <- sum(d)

        if(((betaj11 + ds - a)/(b - a)>1)& (IE == "<="))
            {
                stop(paste("   betaj1 is too large, has to be <= ", round(b-ds,digits = 5)))
            }
        if(((betaj11 + ds - a)/(b - a)>1)& (IE == ">="))
            {
                stop(paste("   betaj1 is too small, has to be >= ", round(ds-b,digits = 5)))
            }

        Z <- XROWS * (Tau_j * Y + d)
        p1 <- (Z - a)/(b - a)

        ## not clear why this was here:
        ## i <- 1
        ## while(i <= XROWS){
        ## if(p1[i]>1) {p1[i] <- 1}
        ## i <- i + 1}
        ## so changed it into this:

        if(max(p1) > 1)
            {
                stop(paste(" largest p1 = ", round(max(p1), 4), ", has to be < 1 "))
            }

        p0 <- 1 - p1

        ## W <- rep(0, times = monte * XROWS)
        ## j <- 1
        ## while(j <= monte)
        ##     {
                ## print(paste("j:", j))
        ##         for(i in 1:XROWS)
        ##             {
                        ## print(paste("i:", i))
                        ## print(paste("i + (j - 1):", i + (j - 1)* XROWS))
        ##                 W[i + (j - 1) * XROWS] <- rbinom(1, size = 1, p1[i])
        ##                 ## W[i + (j-1) * XROWS] <- rbern(1, p1[i])
        ##             }
        ##         j <- j + 1
        ##     }
        ## Wbars <- rep(0, times = monte)
        ## i <- 1
        ## while(i <= monte)
        ##     {
        ##         temp1 <-  (i - 1) * XROWS + 1
        ##         temp2 <- i * XROWS
                ## print(paste("from", temp1, "to", temp2))
        ##         Wbars[i] <- mean(W[temp1:temp2])
        ##         i <- i + 1
        ##     }

        ## much faster than commented code above:
        W <- sapply(1:XROWS, function(x) rbinom(n = monte, size = 1, p1[x]))
        Wbars <- rowMeans(W)

        pbar = (betabarj0 + ds - a)/(b - a)

        ## Bkp is the probability of getting at least k successes
        Bkp <- function(k, p)
            {
                pbinom(k - 1, XROWS, p, lower.tail = FALSE)
            }

        ## in this software we do not fix theta and then find kbar but instead
        ## we search among a set of kk and choose theta such that lamda=0
        ## which means we choose theta such that Bkp(kk,pbar) = theta * alpha
        ## this means that we only need to investigate the tail at kk
        ## which means that we need for iid to be worst case that
        ## kk >= pbar * XROWS + 1

        ## we first look for smallest k st Bkp(k,pbar) <= alpha and
        ## k >= pbar * XROWS + 1
        k1 <- floor(pbar * XROWS + 1)
        while ((Bkp(k1, pbar) > alpha) & (k1 <= XROWS))
            {
                k1 <- k1 + 1
            }
        ## so k1 is this smallest k
        ## if(Bkp(floor(pbar * XROWS),pbar)<=alpha) print (' *  *  *  *  *  OLS Hoeffding binding')

        if(k1 <= XROWS)
            {
                TYPEIIbernoulli <- function(betaj)
                    {
                        if(kbar/XROWS <= pbar)
                            {
                                TYPEIIbernoullitemp <- 1
                            } else {
                                TYPEIIbernoullitemp <- (1 - Bkp(kbar, (betaj + ds - a)/(b - a)))/(1 - theta)
                            }
                        TYPEIIbernoullitemp
                    }

                TYPEII_B <- Inf
                kb <- XROWS
                theta <- Inf
                kk <- k1

                while(Bkp(kk, pbar)/alpha > 0.05 )
                    {
                        alphabar <- Bkp(kk, pbar)
                        theta <- alphabar/alpha
                        kbar <- kk
                        ## print(paste("kbar: ", kbar))
                        ## print(paste("theta", theta))
                        ## print(paste("typeII: ", TYPEIIbernoulli(betaj11)))

                        if(TYPEIIbernoulli(betaj11) < TYPEII_B)
                            {
                                TYPEII_B <- TYPEIIbernoulli(betaj11)
                                kb <- kbar ##value of kbar to remember
                                ##for below
                                ## print(paste("kb: ", kb))

                            }
                        kk <- kk + 1
                        ## print(paste("k1: ", kk))
                    }

                kbar <- kb
                alphabar <- Bkp(kbar,pbar)
                theta <- alphabar/alpha


                ## r_alphaprimeWbar <- function(Wbar)
                ##     {
                ##         if(XROWS * Wbar >= kbar)
                ##             {
                ##                 r_alphaprimeWbartemp <- 1
                ##             }
                ##         else if(XROWS * Wbar == kbar - 1)
                ##             {
                ##                 r_alphaprimeWbartemp <- (alphabar - Bkp(kbar, pbar))/(Bkp(kbar - 1, pbar) - Bkp(kbar, pbar))
                ##             }
                ##         else
                ##             {
                ##                 r_alphaprimeWbartemp <- 0
                ##             }
                ##         r_alphaprimeWbartemp
                ##     }
                ## r_alphaprimeWbars <- rep(0, times = monte)
                ## i <- 1
                ## while(i <= monte)
                ##     {
                ##         r_alphaprimeWbars[i] <- r_alphaprimeWbar(Wbars[i])
                ##         i <- i + 1
                ##     }

                ## vectorized, slightly faster
                r_alphaprimeWbar1 <- function(Wbar)
                    {
                        res <- ifelse(XROWS * Wbar >= kbar, 1,
                                      ifelse(XROWS * Wbar != kbar - 1, 0,
                                             (alphabar - Bkp(kbar, pbar))/(Bkp(kbar - 1, pbar) - Bkp(kbar, pbar))))
                    }
                r_alphaprimeWbars <- r_alphaprimeWbar1(Wbars)


                if(printdetails)
                    {
                        cat("\nstop 5\n")
                        cat("NON-RANDOMISED BERNOULLI TEST (using OLS estimator) \n")
                    }
                if(lambda < 1)
                    {
                        cat(paste("lambda =", lambda," \n"))
                    }

                if(TYPEII_B < 1)
                    {
                        if(mean(r_alphaprimeWbars) >= theta)
                            {
                                if(printdetails)
                                    {
                                        cat(">>> REJECTS H0  (rejection probability of randomized test =",
                                            mean(r_alphaprimeWbars),
                                            "), error=",
                                            exp(-2 * monte * (mean(r_alphaprimeWbars) - theta)^2))
                                    }
                                RejectBernoulliOLS <- "YES"
                            } else {
                                if(printdetails)
                                    {
                                        cat("does not reject H0, rejection probability of randomized test =",
                                            mean(r_alphaprimeWbars),
                                            "), error=",
                                            exp(-2 * monte * (mean(r_alphaprimeWbars) - theta)^2))
                                        cat("\n")
                                    }
                                RejectBernoulliOLS <- "NO"
                            }
                    }

                SummaryBernoulli <- c("theta", "TYPE II bound Bernoulli",
                                      round(theta, digits = 5),
                                      round(min(1, TYPEII_B), digits = 5))
                dim(SummaryBernoulli) <- c(2, 2)


                mydataBernoulli <- as.data.frame(SummaryBernoulli)
                names(mydataBernoulli) <- c("Type II related variables",
                                            "values")
                if(printdetails)
                    {
                        print(mydataBernoulli)
                    }

            } else {
                TYPEII_B <- Inf
                if(printdetails)
                    {
                        cat("stop 6\n")
                        cat("NON-RANDOMISED BERNOULLI TEST (using OLS estimator) \n")
                    }
                if(lambda < 1)
                    {
                        cat(paste("lambda =", lambda," \n"))
                    }
                cat(" .... need to reduce lambda to make Bernoulli applicable under OLS")
                cat("\n")
            }

        ##
        ##
        ##
        ##
        ##
        ##
        ##   Non-Randomised Bernoulli Test (MM)
        ##
        ##
        ##
        ##
        ##
        ##

        a <- 0
        bmm <- XROWS * TAUjmm_inf
        i <- 1
        dmm <- rep(0, times = XROWS)
        while(i <= XROWS)
            {
                dmm[i] <- (1-lambdamm)  *  max(-Taumm_j[i] * ww,
                                               -Taumm_j[i] * (ww + 1)) + lambdamm * (TAUjmm_inf - max(Taumm_j[i] * ww, Taumm_j[i] * (ww + 1)))
                i <- i + 1
            }
        dsmm <- sum(dmm)

        if(((betaj11 + dsmm - a)/(b - a) > 1) & (IE == "<="))
            {
                stop(paste("   betaj1 is too large, has to be <= ",
                           round(b - dsmm,digits = 5)))
            }
        if(((betaj11 + dsmm - a)/(b - a) > 1)&(IE == ">=")) {
            stop(paste("   betaj1 is too small, has to be >= ",
                       round(dsmm - b,digits = 5)))
        }


        Zmm <- XROWS * (Taumm_j * Y + dmm)
        p1mm <- (Zmm - a)/(bmm - a)

        ## i <- 1
        ## while(i <= XROWS)
        ##     {
        ##         if(p1mm[i] > 1)
        ##             {
        ##                 p1mm[i] <- 1
        ##             }
        ##         i <- i + 1
        ##     }
        p1mm[p1mm > 1] <- 1
        p0mm <- 1 - p1mm

        ## Wmm <- rep(0, times = monte * XROWS)
        ## j <- 1
        ## cat("\nwhile 1\n")
        ## while(j <= monte)
        ##     {
        ##         for(i in 1:XROWS)
        ##             {
        ##                 Wmm[i + (j - 1) * XROWS] <- rbinom(1, size = 1, p1mm[i])
        ##                 ## Wmm[i + (j-1) * XROWS] <- rbern(1, p1mm[i])
        ##             }
        ##         j <- j + 1
        ##     }

        ## Wbarsmm <- rep(0, times = monte)
        ## i <- 1
        ## cat("while 2\n")
        ## while(i <= monte)
        ##     {
        ##         temp1mm <- (i-1) * XROWS + 1
        ##         temp2mm <- i * XROWS
        ##         Wbarsmm[i] <- mean(Wmm[temp1mm:temp2mm])
        ##         i <- i + 1
        ##     }
        ## much faster than commented code above:
        Wmm <- sapply(1:XROWS, function(x) rbinom(n = monte, size = 1, p1mm[x]))
        Wbarsmm <- rowMeans(Wmm)

        pbarmm <- (betabarj0 + dsmm - a)/(bmm - a)

        k1mm <- floor(pbarmm * XROWS + 1)
        while ((Bkp(k1mm, pbarmm) > alpha) & (k1mm <= XROWS))
            {
                k1mm <- k1mm + 1
            }
        if(Bkp(floor(pbarmm * XROWS), pbarmm) <= alpha)
            {
                print (' *  *  *  *  *  MM Hoeffding binding')
            }

        if(k1mm <= XROWS)
            {
                TYPEIIbernoullimm <- function(betaj)
                    {
                        if(kbarmm/XROWS <= pbarmm)
                            {
                                TYPEIIbernoullitempmm <- 1
                            } else {
                                TYPEIIbernoullitempmm <- (1 - Bkp(kbarmm, (betaj + dsmm - a)/(bmm - a)))/(1 - thetamm)
                            }
                        TYPEIIbernoullitempmm
                    }

                TYPEII_Bmm <- Inf
                kbmm <- XROWS
                thetamm <- Inf
                kkmm <- k1mm

                while(Bkp(kkmm, pbarmm)/alpha > 0.05)
                    {

                        alphabarmm <- Bkp(kkmm, pbarmm)
                        thetamm <- alphabarmm/alpha
                        kbarmm <- kkmm

                        ## print(thetamm)
                        ## print(TYPEIIbernoullimm(betaj11))

                        if(TYPEIIbernoullimm(betaj11) < TYPEII_Bmm)
                            {
                                TYPEII_Bmm <- TYPEIIbernoullimm(betaj11)
                                kbmm <- kbarmm
                            }
                        kkmm <- kkmm + 1
                    }

                kbarmm <- kbmm
                alphabarmm <- Bkp(kbarmm,pbarmm)
                thetamm <- alphabarmm/alpha

                r_alphaprimeWbarmm <- function(Wbar)
                    {
                        if(XROWS * Wbar >= kbarmm)
                            {
                                r_alphaprimeWbartempmm <- 1
                            } else if(XROWS * Wbar == kbarmm-1) {
                                r_alphaprimeWbartempmm <- (alphabarmm - Bkp(kbarmm, pbarmm))/(Bkp(kbarmm - 1, pbarmm) - Bkp(kbarmm, pbarmm))
                            } else {
                                r_alphaprimeWbartempmm <- 0
                            }
                        r_alphaprimeWbartempmm
                    }

                r_alphaprimeWbarsmm <- rep(0, times = monte)
                i <- 1
                while(i <= monte)
                    {
                        r_alphaprimeWbarsmm[i] <- r_alphaprimeWbarmm(Wbarsmm[i])
                        i <- i + 1
                    }

                if(printdetails)
                    {
                        cat("\n")
                        cat("stop 7\n")
                        cat("NON-RANDOMISED BERNOULLI TEST (using minimax estimator)\n")
                    }
                if(lambdamm < 1)
                    {
                        cat(paste("lambdamm =", lambdamm," \n"))
                    }

                if(TYPEII_Bmm < 1)
                    {
                        if(mean(r_alphaprimeWbarsmm) >= thetamm)
                            {
                                if(printdetails)
                                    {
                                        cat(">>> REJECTS H0  (rejection probability of randomized test =",
                                            mean(r_alphaprimeWbarsmm),
                                            "), error=",
                                            exp(-2 * monte * (mean(r_alphaprimeWbarsmm) - thetamm)^2))
                                    }
                                RejectBernoullimm <- "YES"
                            } else {
                                if(printdetails)
                                    {
                                        cat("does not reject H0, rejection probability of randomized test =",
                                            mean(r_alphaprimeWbarsmm),
                                            "), error=",
                                            exp(-2 * monte * (mean(r_alphaprimeWbarsmm) - thetamm)^2))
                                        cat("\n")
                                    }
                                RejectBernoullimm <- "NO"
                            }
                    }

                SummaryBernoullimm <- c("theta", "TYPE II bound Bernoulli",
                                        round(thetamm, digits = 5),
                                        round(min(1,TYPEII_Bmm), digits = 5))
                dim(SummaryBernoullimm) <- c(2,2)


                mydataBernoullimm <- as.data.frame(SummaryBernoullimm)
                names(mydataBernoullimm) <- c("Type II related variables",
                                              "values")
                if(printdetails)
                    {
                        print(mydataBernoullimm)
                    }

            } else {
                TYPEII_Bmm <- Inf
                if(printdetails)
                    {
                        cat("\n")
                        cat("NON-RANDOMISED BERNOULLI TEST (using minimax estimator)\n")
                    }
                if(lambdamm < 1)
                    {
                        cat(paste("lambdamm =", lambdamm," \n"))
                    }
                cat(" .... need to reduce lambdamm to make Bernoulli applicable under MM")
                cat("\n")
            }

####################################################################################################USER OUTPUT

        ## changed: Error output

        if(printdetails)
            {
                cat("\n")
            }

        if(det(DmatOLS) <= zero)
            {
                print("Upper bound on sigma^2bar used as Dmat not positive definite (OLS)")
            }
        if(det(Dmatmm) <= zero)
            {
                print("Upper bound on sigma^2bar used as Dmat not positive definite (MM)")
            }
        if(k1 > XROWS)
            {
                print(" ... need to decrease lambda to apply Bernoulli under OLS")
            }
        if(k1mm > XROWS)
            {
                print(" ... need to decrease lambdamm to apply Bernoulli under MM")
            }

        ##Selecting the test
        cat("\n")
        if(min(TYPEIINonstandardizedOLS,
               TYPEIINonstandardizedmm,
               TYPEII_B, TYPEII_Bmm) == 1)
            {
                if(IE == "<=")
                    {
                        cat("betaj1 needs to be made larger")
                    } else {
                        cat("betaj1 needs to be made smaller")
                    }
                cat("\n")
            } else {
                if(min(TYPEIINonstandardizedOLS,
                       TYPEIINonstandardizedmm,
                       TYPEII_B, TYPEII_Bmm) == TYPEIINonstandardizedOLS)
                    {
                        ## UserSummary <- c(paste("Estimated coefficient of beta",
                        ##                        cj, ":"),
                        ##                  paste("Reject H0: beta", cj, IE,
                        ##                        betabarj,"?"),
                        ##                  "based on:", "Estimator:",
                        ##                  paste("using",
                        ##                        CutoffNonstandardizedOLS,
                        ##                        "cutoff (t-bar):"),
                        ##                  paste("Type II bound (",
                        ##                        TYPEIINonstandardizedOLSname,
                        ##                        ") for beta", cj, IEa,
                        ##                        round(betaj1,digits = 5)
                        ##                        , ":"), round(IEc * Betahatj,
                        ##                                      digits = 5),
                        ##                  paste(RejectNonstandardizedOLS,
                        ##                        " at alpha = ", alpha),
                        ##                  ModelNonstandardizedOLS, "OLS",
                        ##                  round(tbarmin, digits = 5),
                        ##                  round(TYPEIINonstandardizedOLS,
                        ##                        digits = 5))
                        betahatj <- IEc * Betahatj
                        chosentest <- ModelNonstandardizedOLS
                        estimator <- "OLS"
                        typeII <- TYPEIINonstandardizedOLS
                        theta <- tbarmin
                        rejection <- RejectNonstandardizedOLS
                    }


                if(min(TYPEIINonstandardizedOLS,
                       TYPEIINonstandardizedmm,
                       TYPEII_B, TYPEII_Bmm) == TYPEIINonstandardizedmm)
                    {
                        ## UserSummary <- c(paste("Estimated coefficient of beta",
                        ##                        cj, ":"),
                        ##                  paste("Reject H0: beta", cj, IE,
                        ##                        betabarj,"?"),
                        ##                  "based on:", "Estimator:",
                        ##                  paste("using",
                        ##                        CutoffNonstandardizedmm,
                        ##                        "cutoff (t-bar):"),
                        ##                  paste("Type II bound (",
                        ##                        TYPEIINonstandardizedmmname,
                        ##                        ") for beta", cj, IEa,
                        ##                        round(betaj1,digits = 5),
                        ##                        ":"),
                        ##                  round(IEc * Betahatj, digits = 5),
                        ##                  paste(RejectNonstandardizedmm,
                        ##                        " at alpha = ", alpha),
                        ##                  ModelNonstandardizedmm,
                        ##                  "Minimax",
                        ##                  round(tbarmin, digits = 5),
                        ##                  round(TYPEIINonstandardizedmm,
                        ##                        digits = 5))
                        betahatj <- IEc * Betahatj
                        chosentest <- ModelNonstandardized
                        estimator <- "Minimax"
                        typeII <- TYPEIINonstandardized
                        theta <- tbarmin
                        rejection <- RejectNonstandardizedmm

                    }

                if(min(TYPEIINonstandardizedOLS,
                       TYPEIINonstandardizedmm,
                       TYPEII_B, TYPEII_Bmm) == TYPEII_B)
                    {

                        ## UserSummary <- c(paste("Estimated coefficient of beta",
                        ##                        cj, ":"),
                        ##                  paste("Reject H0: beta", cj, IE,
                        ##                        betabarj, "?"),
                        ##                  "based on:", "Estimator:", "theta:",
                        ##                  paste("Type II bound for beta", cj,
                        ##                        IEa, round(betaj1,digits = 5),
                        ##                        ":"),
                        ##                  round(IEc * Betahatj, digits = 5),
                        ##                  paste(RejectBernoulliOLS,
                        ##                        " at alpha = ", alpha),
                        ##                  "Bernoulli" , "OLS",
                        ##                  round(theta, digits = 5),
                        ##                  round(TYPEII_B, digits = 5))
                        betahatj <- IEc * Betahatj
                        chosentest <- "Bernoulli"
                        estimator <- "OLS"
                        typeII <- TYPEII_B
                        theta <- theta
                        rejection <- RejectBernoulliOLS
                    }

                if(min(TYPEIINonstandardizedOLS,
                       TYPEIINonstandardizedmm,
                       TYPEII_B, TYPEII_Bmm) == TYPEII_Bmm)
                    {
                        ## UserSummary <- c(paste("Estimated coefficient of beta",
                        ##                        cj, ":"),
                        ##                  paste("Reject H0: beta", cj, IE,
                        ##                        betabarj,"?"),
                        ##                  "based on:", "Estimator:", "theta:",
                        ##                  paste("Type II bound for beta", cj,
                        ##                        IEa, round(betaj1, digits = 5),
                        ##                        ":"),
                        ##                  round(IEc * Betahatj, digits = 5),
                        ##                  paste(RejectBernoullimm,
                        ##                        " at alpha = ", alpha),
                        ##                  "Bernoulli" , "Minimax",
                        ##                  round(thetamm, digits = 5),
                        ##                  round(TYPEII_Bmm, digits = 5))
                        betahatj <- IEc * Betahatj
                        chosentest <- "Bernoulli"
                        estimator <- "Minimax"
                        typeII <- TYPEII_Bmm
                        theta <- thetamm
                        rejection <- RejectBernoullimm
                    }

                ## dim(UserSummary) <- c(6,2)
                ## mydataUserSummary <- as.data.frame(UserSummary)
                ## names(mydataUserSummary) <- c(paste("SUMMARY RESULTS, n=",
                ##                                     XROWS, "m=", XCOLS),
                ##                               paste(" Y range [", w1Y,
                ##                                     ",", w2Y, "]"))
                ## print(mydataUserSummary)
            }
        if(min(TYPEIINonstandardizedOLS,
               TYPEIINonstandardizedmm,
               TYPEII_B, TYPEII_Bmm) < 0.49)
            {
                if(IE == "<=")
                    {
                        cat("...... decrease betaj1 to get type II closer to 0.5 \n")
                    } else {
                        cat("...... increase betaj1 to get type II closer to 0.5 \n")
                    }
            }
        if(min(TYPEIINonstandardizedOLS,
               TYPEIINonstandardizedmm,
               TYPEII_B, TYPEII_Bmm) > 0.51)
            {
                if(IE == "<=")
                    {
                        cat("...... increase betaj1 (if possible) to get type II closer to 0.5 \n")
                    } else {
                        cat("...... decrease betaj1 (if possible) to get type II closer to 0.5 \n")
                    }
            }

        method <- "Exact linear models"
        ## rejection <- ifelse
        bounds <- paste("[", w1Y, ", ", w2Y, "]", sep = "")

        OLSTestsTypeI <- list() ## TypeI errors
        OLSTestsTypeII <- list() ## TypeII errors
        OLSTests <- list(OLSTestsTypeI, OLSTestsTypeII)

        MMTestsTypeI <- list()
        MMTestsTypeII <- list()
        MMTests <- list(MMTestsTypeI, MMTestsTypeII)

        ## a copy of the results of the best test
        chosentest <- list(test = "Name of test",
                           estimator = "OLS or MM",
                           typeI = 0.55,
                           typeII = 0.44,
                           theta = 0.33)
        OLSNrBernoulli <- list() ## Nonrandomized Bernoulli Test OLS Parameter
        MMNrBernoulli <- list() ## Nonrandomized Bernoulli Test MM
        ## Parameter

        ## parameter <- list(n = XROWS,
        ##                   m = XCOLS,
        ##                   j = cj,
        ##                   alternative = IE,
        ##                   alpha = alpha,
        ##                   bounds = bounds,
        ##                   iterations = monte)

        ## structure(list(method = method,
        ##                yname = YNAME,
        ##                xname = XNAME,
        ##                parameter = parameter,
        ##                coef = cj,
        ##                betabarj = betabarj0,
        ##                ## betahatj = betahatj,
        ##                chosentest = chosentest,
        ##                estimator = estimator,
        ##                typeII = typeII,
        ##                theta = theta,
        ##                rejection = rejection,
        ##                bounds = bounds,
        ##                theta = theta),
        ##           class = "elm")

    }
