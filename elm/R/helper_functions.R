## calculate sigmasqbar-variants (one for the typeI, one for typeII)
calcSigmasqbar <- function(X, ww, XROWS, XCOLS, ej,
                           tau_jB,
                           tauj_2, tau_j,
                           betabarj, betaj,
                           type = "typeI",
                           zero = 10^-12)
    {
        sigmasqbar <- vector(length = 2)
        names(sigmasqbar) <- c("TypeI", "TypeII")

        Dmat <- 2 * t(X) %*% ((tau_jB^2) * X)
        dvec <- (1 + 2 * ww) * colMeans((tau_jB^2) * X) * XROWS
        Amat <- t(matrix(rbind(-t(ej), X, -X), ncol = XCOLS))

        if(det(Dmat) <= zero)
            {
                ## typeI
                warning("Dmat is nearly singular!")
                if(type == "typeI")
                    {
                        term1 <- min(0, betabarj - (1/2) * sum(tau_j))
                        sigmasqbar[1] <- tauj_2/4 - (1/XROWS) * (term1)^2
                    }
                ## typeII
                else if (type == "typeII")
                    {
                        term1 <- betaj - (1/2) * sum(tau_j)
                        sigmasqbar[2] <- tauj_2/4 - (1/XROWS) * (term1)^2
                    }
            } else {
                ## typeI
                if(type == "typeI")
                    {
                        bvec <- as.vector(c(-betabarj,
                                            rep(ww, times = XROWS),
                                            rep(-(ww + 1), times = XROWS)))
                        meq <- 0
                        zsol <- solve.QP(Dmat = Dmat,
                                         dvec = dvec,
                                         Amat,
                                         bvec = bvec, meq = meq,
                                         factorized = FALSE)$solution
                        sigmasqbar[1] <- sum((tau_j)^2 * ((X %*% zsol - ww) * (1 + ww - X %*% zsol)))
                    }
                else if (type == "typeII")
                    {
                        ## typeII
                        c2 <- 1/max(Dmat) ## NEW
                        bvec <- as.vector(c(-betaj,
                                            rep(ww, times = XROWS),
                                            rep(-(ww + 1), times = XROWS)))
                        meq <- 1
                        zsol <- solve.QP(Dmat = c2 * Dmat,
                                         dvec = c2 * dvec,
                                         Amat,
                                         bvec = bvec, meq = meq,
                                         factorized = FALSE)$solution
                        ## if there is an error message
                        ## Fehler in solve.QP(Dmat = c2 * DmatOLS, dvec = c2 * dvecOLS, Amat = Amat,  :
                        ##  constraints are inconsistent, no solution!
                        ## then this often results to numerical problems when solving (2) in
                        ## (Gossner and Schlag, 2013)
                        ## but can also come if betaj is not close enough to
                        ## 0
                        sigmasqbar[2] <- sum((tau_j)^2 * ((X %*% zsol - ww) * (1 + ww - X %*% zsol)))
                    }
            }
        sigmasqbar
    }
