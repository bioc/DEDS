\name{comp.F}
\alias{comp.F}
\title{Computing F-statistic for Differential Expression}
\description{
  \code{comp.F} returns a function of one argument with bindings for
  \code{L}. This function accepts a microarray data matrix as its single
  argment, when evaluated, computes F statistics for each row of the
  matrix.
}

\usage{
comp.F(L = NULL)
}
\arguments{
  \item{L}{A vector of integers corresponding to observation (column)
    class labels. For \eqn{k} classes, the labels must be integers
    between 0 and \eqn{k-1}.
  }
}

\value{
  \code{comp.F} returns a function with bindings for \code{L},
  which calculates and returns of vector of F statistics  for
  each row in the data matrix.
}

\author{Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.
  }

\seealso{\code{\link{comp.FC}}, \code{\link{comp.t}}}
\examples{

X <- matrix(rnorm(1000,0,0.5), nc=10)
L <- rep(0:1,c(5,5))

# genes 1-10 are differentially expressed
X[1:10,6:10]<-X[1:10,6:10]+1

# two sample test, unequal variance
F <- comp.F(L)
f.X <- F(X)
}

\keyword{univar}




