\name{comp.t}
\alias{comp.t}
 
\title{Computing One and Two Sample t-statistic for Differential Expression }
\description{
  \code{comp.t} returns a function of one argument with bindings for
  \code{L}, \code{mu}, \code{var.equal}. This function accepts a
  microarray data matrix as its single argment, when evaluated, computes
  t statistics for each row of the matrix.
}

\usage{
comp.t(L = NULL, mu = 0, var.equal = FALSE)
}
\arguments{
  \item{L}{A vector of integers corresponding to observation (column)
    class labels. For \eqn{k} classes, the labels must be integers
    between 0 and \eqn{k-1}.
  }
  \item{mu}{A number indicating the true value of the mean (or
    difference in means if you are performing a two sample test).
   }
   \item{var.equal}{a logical variable indicating whether to treat the
     two variances as being equal. If \code{TRUE} then the pooled
     variance is used to estimate the variance otherwise the Welch
     statistic will be calculated.
   }
}
\details{
  The function returned by \code{comp.t} calculates t statistics for
  each row of the microarary data matrix, given specific class labels. 
}
\value{
  \code{comp.t} returns a function with bindings for \code{L},
  \code{mu}, \code{var.equal}, which calculates and returns of vector of
  t statistics for each row in the data matrix.
}

\author{Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.
  }

\seealso{\code{\link{comp.FC}}, \code{\link{comp.F}}}
\examples{
X <- matrix(rnorm(1000,0,0.5), nc=10)
L <- rep(0:1,c(5,5))

# genes 1-10 are differentially expressed
X[1:10,6:10]<-X[1:10,6:10]+1

# two sample test, unequal variance
t.fun <- comp.t(L)
t.X <- t.fun(X)

# two sample test, equal variance
t.fun <- comp.t(L, var.equal=TRUE)
t.X <- t.fun(X)
}

\keyword{univar}



