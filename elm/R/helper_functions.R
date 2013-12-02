## calculate sigmasqbar-variants (one for the typeI, one for typeII)
calcSigmasqbar <- function(X, ww, XROWS, XCOLS, ej,
                           Tau_jB,
                           TAUj_2, Tau_j,
                           betabarj, betaj1,
                           zero = 10^-12)
    {
        sigmasqbar <- vector(length = 2)
        names(sigmasqbar) <- c("TypeI", "TypeII")

        Dmat <- 2 * t(X) %*% ((Tau_jB^2) * X)
        dvec <- (1 + 2 * ww) * colMeans((Tau_jB^2) * X) * XROWS
        Amat <- t(matrix(rbind(-t(ej), X, -X), ncol = XCOLS))

        if(det(Dmat) <= zero)
            {
                ## typeI
                term1 <- min(0, betabarj - (1/2) * sum(Tau_j))
                sigmasqbar[1] <- TAUj_2/4 - (1/XROWS) * (term1)^2

                ## typeII
                term1 <- betaj1 - (1/2) * sum(Tau_j)
                sigmasqbar[2] <- TAUj_2/4 - (1/XROWS) * (term1)^2
            } else {
                ## typeI
                bvec <- as.vector(c(-betabarj,
                                    rep(ww, times = XROWS),
                                    rep(-(ww + 1), times = XROWS)))
                meq <- 0
                zsol <- solve.QP(Dmat = Dmat,
                                 dvec = dvec,
                                 Amat,
                                 bvec = bvec, meq = meq,
                                 factorized = FALSE)$solution
                sigmasqbar[1] <- sum((Tau_j)^2 * ((X %*% zsol - ww) * (1 + ww - X %*% zsol)))

                ## typeII
                c2 <- 1/max(Dmat) ## NEW
                bvec <- as.vector(c(-betaj1,
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
                ## but can also come if betaj1 is not close enough to
                ## 0
                sigmasqbar[2] <- sum((Tau_j)^2 * ((X %*% zsol - ww) * (1 + ww - X %*% zsol)))
            }
        sigmasqbar
    }
