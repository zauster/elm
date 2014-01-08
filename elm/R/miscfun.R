
is.constant1 <- function(x)
    {
        length(unique(x)) == 1
    }

is.constant2 <- function(x)
    {
        sd(x, na.rm = TRUE) == 0
    }

is.constant3 <- function(x)
    {
        all(is.na(x)) | all(x[1L] == x)
    }

calcMaxBetaj <- function(beta, Amat, bvec, ej)
    {
        bvec <- c(-beta[2], bvec)
        sum(Amat %*% beta - bvec)^2
    }

is.error <- function(x) inherits(x, "try-error")
