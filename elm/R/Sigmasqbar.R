## calculate sigmasqbar-variants (one for the typeI, one for typeII)
calcSigmasqbar <- function(X, ww, XROWS, XCOLS, ej,
                           tau_jB,
                           tauj_2, tau_j,
                           Dmat, dvec, Amat,
                           betabarj, betaj,
                           type = "typeI",
                           silent = TRUE,
                           zero = 10^-12)
    {
        sigmasqbar <- vector(length = 2)
        names(sigmasqbar) <- c("TypeI", "TypeII")

        ## sigmas <<- sigmas + 1 ## to count all calls of calcSigma
        ## Dmat <- 2 * t(X) %*% ((tau_jB^2) * X)
        ## dvec <- (1 + 2 * ww) * colMeans((tau_jB^2) * X) * XROWS
        ## Amat <- t(matrix(rbind(-1 * t(ej), X, -X), ncol = XCOLS))

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
                else if(type == "typeII")
                    {
                        term1 <- betaj - (1/2) * sum(tau_j)
                        sigmasqbar[2] <- tauj_2/4 - (1/XROWS) * (term1)^2
                    }
            } else {
                ## typeI
                if(type == "typeI")
                    {
                        bvec <- as.vector(c(-1 * betabarj,
                                            rep(ww, times = XROWS),
                                            rep(-(ww + 1), times = XROWS)))
                        meq <- 0
                        zsol <- try(solve.QP(Dmat = Dmat,
                                             dvec = dvec,
                                             Amat,
                                             bvec = bvec, meq = meq,
                                             factorized = FALSE)$solution,
                                    silent = silent)
                        if(is.error(zsol) == TRUE)
                            {
                                if(grepl("constraints are inconsistent", zsol[1]) == TRUE)
                                    {
                                        warning("Could not find the optimal value for sigmasqbar. Will use upper bound instead.")
                                        term1 <- min(0, betabarj - (1/2) * sum(tau_j))
                                        sigmasqbar[1] <- tauj_2/4 - (1/XROWS) * (term1)^2
                                    }
                                else
                                    {
                                        warning("WARNING: Unusual error in finding optimal sigmasqbar, will use upper bound instead. Are the parameters given sensible?")
                                        term1 <- min(0, betabarj - (1/2) * sum(tau_j))
                                        sigmasqbar[1] <- tauj_2/4 - (1/XROWS) * (term1)^2
                                    }
                            }
                        else
                            {
                                sigmasqbar[1] <- sum((tau_j)^2 * ((X %*% zsol - ww) * (1 + ww - X %*% zsol)))
                            }
                    }
                else if(type == "typeII")
                    {
                        ## typeII
                        ## browser()
                        c2 <- 1/max(Dmat) ## NEW
                        ## cat("\nc2: ", c2)
                        bvec <- as.vector(c(-1 * betaj,
                                            rep(ww, times = XROWS),
                                            rep(-(ww + 1), times = XROWS)))
                        ## cat("\nbvec: ", bvec)
                        meq <- 1
                        zsol <- try(solve.QP(Dmat = c2 * Dmat,
                                             dvec = c2 * dvec,
                                             Amat,
                                             bvec = bvec, meq = meq,
                                             factorized = FALSE)$solution,
                                    silent = silent)
                        if(is.error(zsol) == TRUE)
                            {
                                if(grepl("constraints are inconsistent", zsol[1]) == TRUE)
                                    {
                                        warning("Could not find the optimal value for sigmasqbar. Will use upper bound instead.")
                                        term1 <- betaj - (1/2) * sum(tau_j)
                                        sigmasqbar[2] <- tauj_2/4 - (1/XROWS) * (term1)^2
                                    }
                                else
                                    {
                                        warning("WARNING: Unusual error in finding optimal sigmasqbar, will use upper bound instead. Are the parameters given sensible?")
                                        term1 <- betaj - (1/2) * sum(tau_j)
                                        sigmasqbar[2] <- tauj_2/4 - (1/XROWS) * (term1)^2
                                    }
                            }
                        else
                            {
                                sigmasqbar[2] <- sum((tau_j)^2 * ((X %*% zsol - ww) * (1 + ww - X %*% zsol)))
                            }

                        ## print(str(zsol))
                        ## cat("\nzsol: ", zsol)
                        ## if there is an error message
                        ## Fehler in solve.QP(Dmat = c2 * DmatOLS, dvec = c2 * dvecOLS, Amat = Amat,  :
                        ##  constraints are inconsistent, no solution!
                        ## then this often results to numerical problems when solving (2) in
                        ## (Gossner and Schlag, 2013)
                        ## but can also come if betaj is not close enough to
                        ## 0
                    }
            }
        sigmasqbar
    }
