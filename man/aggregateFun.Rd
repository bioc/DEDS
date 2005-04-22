\name{aggregateFun}
\alias{aggregateFun}

\title{Aggregate Statistical Functions for DEDS}
\description{
  This function takes statistic functions and  creates a
  function that takes a matrix as a single argument. The statisic
  functions are bound in the environment of the returned function and
  are applied sequentially to the argument of the returned function.
}
\usage{
aggregateFun(...)
}

\arguments{
  \item{\dots}{Functions of various statistics, could be in a list. }
}
\details{
  The function takes several statistics functions or a list of these
  functions and returns a function (F) with bindings to the input
  statistics functions. F takes a data matrix as its single argument,
  and apply the bound statistical functions sequentially to the data
  matrix.
}
\value{
  It returns a function that takes a matrix as its single argument. The
  function returns a matrix of statistics, with \eqn{m} rows
  correspondign to variables (hypotheses) and \eqn{n} columns
  corresponding to specified statistics.
}


\author{
  Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.
}

\seealso{\code{\link{comp.t}}, \code{\link{comp.FC}}, \code{\link{comp.F}},
  \code{\link{comp.SAM}}, \code{\link{comp.modt}}, \code{\link{comp.modF}},
  \code{\link{comp.B}}}

\examples{
X <- matrix(rnorm(100, 0, 1), nc=10)
L <- rep(0:1,c(5,5))
t.fun <- comp.t(L)
fc.fun <- comp.FC(L)
sam.fun <- comp.SAM(L)
ffun <- aggregateFun(list(t.fun, fc.fun, sam.fun))
stats <- ffun(X)

}

\keyword{manip}






