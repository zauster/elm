\name{elm}
\alias{elm}
\title{
Exact linear models
}
\description{
  Uses exact tests for the coefficients of linear regressions.
}
\usage{
elm(Y, X, lower = 0, upper = 1,
    alternative = "greater",
    alpha = 0.05,
    coefs = 2, 
    nullvalue = 0, 
    upperbetabound = 1,
    lambda = 1, lambdamm = 1,
    qq = 0.0001, qqmm = 0.0001,
    iterations = 1000,
    steppc = 0.1,
    silent = FALSE,
    verbose = TRUE,
    na.action = getOption("na.action"))
}

\arguments{
  \item{Y}{
    dependent variable, as matrix.
}
  \item{X}{
    independent variable, as matrix.
}
  \item{lower, upper}{
    the theoretical lower and upper bounds on the data outcomes known ex-ante before gathering the data.
}
  \item{alternative}{
    the hypothesis to be tested, "less" or "greater" (default).
}
  \item{alpha}{
    the type I error.
}
  \item{coefs}{
    index of the coefficient to be tested
}
  \item{nullvalue}{
    the critical value for the null hypothesis
}
  \item{upperbetabound}{
    the upper bound of beta in the set of the alternative
  hypothesis. The program tries to find a beta in [\code{nullvalue},
  \code{upperbetabound}] which brings the typeII error to 0.5. If
  \code{upperbetabound} is set to \code{NULL}, it will try to guess a
  usable \code{upperbetabound} and slowly increase (\code{steppc}
  controls the increases) it until it finds an
  optimal beta. This could be, however, computationally expensive.
}
\item{steppc}{
  Controls the size of the steps taken in finding the optimal beta. The
  stepwise increase is \code{upperbetabound} * \code{steppc}. Default is
  0.1. 
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{lambdamm}{
%%     ~~Describe \code{lambdamm} here~~
}
  \item{iterations}{
    number of iterations
}
  \item{qq}{
%%     ~~Describe \code{qq} here~~
}
  \item{qqmm}{
%%     ~~Describe \code{qqmm} here~~
}
\item{silent}{
  Should warnings during the procedure be displayed? Default is FALSE.
}
\item{verbose}{
  If FALSE, it prints only essential summary of the test. Default is TRUE.
}
\item{na.action}{
  How to cope with missing values. Uses system-default as default value.
}
}
\details{
  This function computes several exact tests for the coefficient of a
  linear regression. For an explanation as to how the tests are
  constructed, please refer to the paper mentioned below.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Olivier Gossner, Karl H. Schlag, "Finite-sample exact tests for linear regressions with bounded dependent variables", Journal of Econometrics, Volume 177, Issue 1, November 2013, Pages 75-84, ISSN 0304-4076, \url{http://dx.doi.org/10.1016/j.jeconom.2013.06.003}.
}
\author{
Karl Schlag, Olivier Gossner, Gareth Liu-Evans and Oliver Reiter
}
%% \note{
%% %%  ~~further notes~~
%% }

\seealso{
  \url{http://homepage.univie.ac.at/karl.schlag/research/statistics.html}
}
\examples{
## step example
n <- 40
h <- 0.5
Y <- sample(c(0, 1), size = n, replace = TRUE)
X <- cbind(1, runif(n = n) < h)
elm(Y, X, 0, 1, coefs = 2, nullvalue = 0, upperbetabound = NULL)
}
\keyword{linear regression}
\keyword{exact method}

