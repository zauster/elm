\name{elm}
\alias{elm}
\title{
Exact linear models
}
\description{
  Uses exact tests for the coefficients of linear regressions.
}
\usage{
elm(Y, X, w1Y, w2Y, IE = "<=", alpha = 0.05, j = 2,
    betabarj = 0, betaj1 = betabarj + 1.1,
    lambda = 1, lambdamm = 1, monte = 1000,
    qq = 1e-04, qqmm = 1e-04, printdetails = TRUE)
}

\arguments{
  \item{Y}{
    dependent variable, as matrix.
}
  \item{X}{
    independent variable, as matrix.
}
  \item{w1Y, w2Y}{
    the theoretical lower and upper bounds on the data outcomes known ex-ante before gathering the data.
}
  \item{IE}{
    the hypothesis to be tested, ">=" or "<=" (default).
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
  \item{betaj1}{
    value of the coefficient under the alternative. to bring the type II
  error to 0.5
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{lambdamm}{
%%     ~~Describe \code{lambdamm} here~~
}
  \item{monte}{
    number of iterations
}
  \item{qq}{
%%     ~~Describe \code{qq} here~~
}
  \item{qqmm}{
%%     ~~Describe \code{qqmm} here~~
}
  \item{printdetails}{
    logical, should details of the tests be printed out?
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}

