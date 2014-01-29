\name{elmCI}
\alias{elmCI}
\title{
Confidence intervals for exact linear models
}
\description{
  Computes the bounds of confidence intervals from a run of \code{elm}
}
\usage{
elmCI(elm, alpha = 0.05, coefs = 2,
      dispWarnings = FALSE,
      useTest = NULL)
}

\arguments{
  \item{elm}{
    test results from a previous \code{elm}-run.
}
\item{alpha}{
  the confidence level of the interval is 1 - \code{alpha}
}
\item{coefs}{
  the index of the coefficients to be tested
}
  \item{dispWarnings}{
    indicator whether warnings should be displayed during the run of
  \code{elmCI}. Default is \code{FALSE}, since the root-finding
  algorithms sometimes searches out-of-bounds, which can result in a
  large amount of warnings.
}
  \item{useTest}{
    By default (\code{NULL}), the test that minimizes the typeII error in the run of
  \code{elm} is used to calculate the bound of the confidence
  interval. If you have doubts about the location of the bound, you can
  use this switch to coerce the program to use one of the four tests
  available. "1" is the Nonstandardized OLS, "2" the Bernoulli OLS, "3"
  the Nonstandardized MM and "4" corresponds to the Bernoulli MM
  test. After using this option, it is recommended to rerun \code{elm}
  with the bound specified by the user and compare it to the bound by
  the original test. See the examples for an illustration of this point (TODO).
}
}
\details{
  This program computes the bound of a confidence interval, given the
  results of a run of \code{elm}. At this stage, it is only possible to
  derive the bound corresponding to the alternative chosen in the run of
  \code{elm}. That means, if \code{alternative = "greater"} was
  specified, the program will look for the lower bound of the regression
  coefficient. For the upper bound, the results of a run of \code{elm}
  with \code{alternative = "less"} will be needed. Two sided confidence
  interval are as of now not implemented. Please run the program at
  \code{alpha/2} with both alternatives.
}
%% \value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%% }
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
h <- 0.3
Y <- sample(c(0, 1), size = n, replace = TRUE)
X <- cbind(1, runif(n = n) < h)

## first lower bound:
result.lower <- elm(Y, X, 0, 1, coefs = 2, nullvalue = 0, upperbetabound =
  NULL, alternative = "greater")
## see what the results look like
result.lower
## compute lower bound for the coefficient
elmCI(result.lower, alpha = 0.025)

## now upper bound:
result.upper <- elm(Y, X, 0, 1, coefs = 2, nullvalue = 0, upperbetabound =
  NULL, alternative = "less")
## see what the results look like
result.upper
## compute upper bound for the coefficient
elmCI(result.upper, alpha = 0.025)

}
\keyword{linear regression}
\keyword{exact method}

