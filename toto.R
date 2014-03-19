
## toto (OLS Nonstandardized)
toto <- read.csv2("toto1N exam.csv")

y <- toto[, 1]
x <- data.matrix(cbind(1, toto[, 2:21]))
xx <- data.matrix(toto[, 2:22])


elm(y, x, coefs = 2, alpha = 0.025)
elmo(y, x, betaj1 = 0.1405107, alpha = 0.025)

elm(y, xi, coefs = 2, upperbetabound = NULL, intercept = F)
elmo(y, xi, coefs = 2, betaj1 = 0.1511815)

## CI testing
safi2g <- elm(y, x, coefs = 2:6, alpha = 0.025, verbose = F)
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


## Schlag Test
X <- as.matrix(read.table("T4A1X.dat"))
Y <- as.matrix(read.table("T4A1Y.dat"))
elm(Y,X,intercept=T,alpha=0.025,nullvalue=0.2,coefs=1)
safi2o <- elm(Y, X, alpha = 0.025, coefs = 1:5, intercept = F)
elmCI(safi2o)
elmo(Y, X, 0, 1, alpha = 0.025, coefs = 1, betabarj = 0.01, betaj1 = 0.1903)



> elmCI(safi2g)

        Confidence Intervals for exact linear models

data: y and x
n = 876, m = 21
Number of tested coefficients: 5

           Coefficient   lower   upper             used Test
             safi_lr04  0.0066  0.2217 Nonstandardized (OLS)
                 sbsk1 -0.0671  0.1861 Nonstandardized (OLS)
            sbsk1_demo -0.2103  0.1576 Nonstandardized (OLS)
                  demo -0.9504  0.9631 Nonstandardized (OLS)
 house_ever_anyfert_05  0.2748  0.4627 Nonstandardized (OLS)

elm(y, x, alpha = 0.025, coefs = 2, nullvalue = 0.0066)
elmo(y, x, alpha = 0.025, coefs = 2, betabarj = 0.007, betaj1 = 0.1470979)
elm(y, x, alpha = 0.025, coefs = 2, nullvalue = 0.2217, alternative = "less")
elmo(y, x, alpha = 0.025, coefs = 2, betabarj = 0.222, IE = ">=", betaj1 = 0.0814)
