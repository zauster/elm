\name{elm}
\alias{elm}
\title{
Exact linear models
}
\description{
  Uses exact tests for the coefficients of linear regressions.
}
\usage{
elm(Y, X, lower = 0, upper = 1, alternative = "greater",
    alpha = 0.05, j = 2, betabarj = 0, betaj = 0.5,
    lambda = 1, lambdamm = 1, iterations = 1000,
    qq = 0.0001, qqmm = 0.0001)
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
  \item{j}{
    index of the coefficient to be tested
}
  \item{betabarj}{
    the critical value for the null hypothesis
}
  \item{betaj}{
    value of the coefficient under the alternative. to bring the type II
  error to 0.5
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
elm(Y, X, 0, 1, j = 2, betabarj = 0, betaj = .43)
}
\keyword{linear regression}
\keyword{exact method}

