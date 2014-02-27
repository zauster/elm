
## toto (OLS Nonstandardized)
toto <- read.csv2("toto1N exam.csv")

y <- toto[, 1]
x <- data.matrix(cbind(1, toto[, 2:21]))
x <- -data.matrix(cbind(1, toto[, 2:21]))


elm(y, x, coefs = 2, alpha = 0.025)
elmo(y, x, betaj1 = 0.1405107, alpha = 0.025)

elm(y, xi, coefs = 2, upperbetabound = NULL, intercept = F)
elmo(y, xi, coefs = 2, betaj1 = 0.1511815)

## CI testing
safi2g <- elm(y, x, coefs = 2:6, alpha = 0.025)
safi2l <- elm(y, x, coefs = 2, alpha = 0.025, alternative = "less")
safi3g <- elm(y, x, coefs = 3, alpha = 0.025)
safi3l <- elm(y, x, coefs = 3, alpha = 0.025, alternative = "less")
safi4g <- elm(y, x, coefs = 4, alpha = 0.025)
safi4l <- elm(y, x, coefs = 4, alpha = 0.025, alternative = "less")
safi5g <- elm(y, x, coefs = 5, alpha = 0.025)
safi5l <- elm(y, x, coefs = 5, alpha = 0.025, alternative = "less")
safi6g <- elm(y, x, coefs = 6, alpha = 0.025)
safi6l <- elm(y, x, coefs = 6, alpha = 0.025, alternative = "less")

elmCI(safi2g)
elmCI(safi2l)
elmCI(safi3g)
elmCI(safi3l)
elmCI(safi4l)
elmCI(safi5g)
elmCI(safi6g)

## OLS Bernoulli
yy <- runif(40) > 0.4
xx1 <- runif(40) > 0.4
xx2 <- cbind(1, xx1)
elm(yy, xx2, upperbetabound = NULL)
elmo(yy, xx2, betaj1 = 0.466101)

## MM Bernoulli
ymm <- runif(600)
xmm <- cbind(1, rnorm(600))
elm(ymm, xmm, upperbetabound = NULL)
elmo(ymm, xmm, betaj1 = 0.06718927)

## MM Nonstandardized
## (Problem with the Bernoulli Test (1 > 1))
yym <- runif(200) > .1
xxm <- cbind(1, rnorm(200), rbeta(200, 2, 4))
res <- elm(yym, xxm, coefs = 2:3, upperbetabound = NULL)
elmo(yym, xxm, betaj1 = 0.1026914)


fun <- function(x, y)
    {
        res <- match.call()

        print(x*y)

        return(res)
    }
