
## Nonrandomized Bernoulli Test
calcBernoulliTest <- function(Y, XROWS, TAUj_inf, Tau_j,
                              ww, kbar, pbar, alphabar,
                              d, ds, b, a = 0, # betaj,
                              betabarj, alternative, iterations,
                              theta, lambda = 1)
    {
        ## if(((betaj + ds - a)/(b - a) > 1) & (alternative == "greater"))
        ##     {
        ##         stop(paste("betaj is too large, has to be <= ", round(b - ds, digits = 5)))
        ##     }
        ## if(((betaj + ds - a)/(b - a) > 1) & (alternative == "less"))
        ##     {
        ##         stop(paste("betaj is too small, has to be >= ", round(ds - b, digits = 5)))
        ##     }

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
                stop(paste("largest p1 = ", round(max(p1), 4), ", has to be < 1 "))
            }

        ## p0 <- 1 - p1

        W <- sapply(1:XROWS, function(x) rbinom(n = iterations, size = 1, p1[x]))
        Wbars <- rowMeans(W)
        rej <- mean(r_alphaprimeWbar1(Wbars, XROWS, kbar, pbar, alphabar))
        error <- exp(-2 * iterations * (rej - theta)^2)

        return(c(probability_rejection = rej,
                    error = error,
                    theta = theta))
        ## to return:
        ## theta
        ## rej
        ## error
    }

calcTypeIIBernoulli <- function(betaj, betabarj, alpha, ds, b, XROWS, a = 0)
    {
        pbar = (betabarj + ds - a)/(b - a)

        ## in this software we do not fix theta and then find kbar but instead
        ## we search among a set of kk and choose theta such that lamda=0
        ## which means we choose theta such that Bkp(kk,pbar) = theta * alpha
        ## this means that we only need to investigate the tail at kk
        ## which means that we need for iid to be worst case that
        ## kk >= pbar * XROWS + 1

        ## we first look for smallest k st Bkp(k,pbar) <= alpha and
        ## k >= pbar * XROWS + 1

        k1 <- floor(pbar * XROWS + 1)
        k1 <- k1 + sum(Bkp(k1:XROWS, pbar, XROWS) > alpha)
        ## above is the same as below, but vectorized
        ## while((Bkp(k1, pbar) > alpha) & (k1 <= XROWS))
        ##     {
        ##         k1 <- k1 + 1
        ##     }
        ## so k1 is this smallest k

        if(k1 <= XROWS)
            {
                TYPEII_B <- Inf
                kbar <- XROWS
                theta <- Inf
                ## kk <- k1

                while(Bkp(k1, pbar, XROWS)/alpha > 0.05)
                    {
                        alphabar <- Bkp(k1, pbar, XROWS)
                        theta <- alphabar/alpha
                        ## kbar <- k1
                        ## print(paste("theta: ", theta))
                        typeIItemp <- TYPEIIbernoulliNew(betaj, k1, pbar, ds, b, XROWS, theta)
                        ## print(paste("typeII: ", typeIItemp))

                        if(typeIItemp < TYPEII_B)
                            {
                                TYPEII_B <- typeIItemp
                                kbar <- k1 ##value of kbar to remember
                                ##for below
                                ## print(paste("kb: ", kbar))

                            }
                        k1 <- k1 + 1
                        ## print(paste("k1: ", k1))
                    }

                ## kbar <- kb
                alphabar <- Bkp(kbar, pbar, XROWS)
                theta <- alphabar/alpha
                ## print(c(kbar, alphabar, theta))
            }
        else
            {
                TYPEII_B <- Inf
                theta <- Inf
            }

        return(list(typeII = ifelse(TYPEII_B < 1, TYPEII_B, 1),
                    theta = theta,
                    kbar = kbar,
                    pbar = pbar,
                    alphabar = alphabar))
    }

## calcTypeIIBernoulliNew <- function(betaj, ds, b, kbar, theta)
##     {
##         pbetaj <- (betaj + ds)/b
##         (1 - calcLambda(theta, alpha, kbar, pbetaj) * Bkp(kbar - 1, pbetaj) - (1 - calcLambda(theta, alpha, kbar, pbetaj)) * Bkp(kbar, pbetaj))/(1 - theta)
##     }

## calcLambda <- function(theta, alpha, kbar, pbar)
##     {
##         min((theta * alpha - Bkp(kbar, pbar))/(Bkp(kbar - 1, pbar) - Bkp(kbar, pbar)), 1)
##     }

## Bkp is the probability of getting at least k successes
Bkp <- function(k, p, n)
    {
        pbinom(k - 1, n, p, lower.tail = FALSE)
    }

TYPEIIbernoulliNew <- function(betaj, kbar, pbar, ds, b, XROWS, theta, a = 0)
    {
        ## ifelse is slightly faster than if/else
        ifelse(kbar/XROWS <= pbar, 1,
               (1 - Bkp(kbar, (betaj + ds - a)/(b - a), XROWS))/(1 - theta))
    }

## vectorized, slightly faster
r_alphaprimeWbar1 <- function(Wbar, XROWS, kbar, pbar, alphabar)
    {
        res <- ifelse(XROWS * Wbar >= kbar, 1,
                      ifelse(XROWS * Wbar != kbar - 1, 0,
                             (alphabar - Bkp(kbar, pbar, XROWS))/(Bkp(kbar - 1, pbar, XROWS) - Bkp(kbar, pbar, XROWS))))
    }