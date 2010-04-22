\name{comp.FC}
\alias{comp.FC}

\title{Computing Fold Change for Differential Expression}
\description{
  \code{comp.FC} returns a function of one argument with bindings for
  \code{L}, \code{is.log} and \code{FUN}. This function accepts a
  microarray data matrix as its single argument, when evaluated, computes
  fold change for each row of the matrix.
}

\usage{
comp.FC(L = NULL, is.log = TRUE, FUN = mean)
}
\arguments{
  \item{L}{A vector of integers corresponding to observation (column)
    class labels. For \eqn{k} classes, the labels must be integers
    between 0 and \eqn{k-1}.}
  \item{is.log}{A logical variable indicating whether the data has been
    logged.}
  \item{FUN}{The summary statistics function used to calculate fold
    change, the default is set as \code{\link{mean}}, the user can also
    use \code{\link{median}}}.
}

\details{
  The function returned by \code{comp.FC} calculates fold change for each row
  of the matrix, given specific class labels. If \code{is.log=TRUE},
  fold change is calculated by subtraction; if \code{is.log=FALSE},
  fold change is calculated by division. 
}
\value{
  \code{comp.FC} returns a function with bindings for \code{L},
  \code{is.log} and \code{FUN}, which calculates and returns a vector of
  fold changes for each row in the data matrix.
}


\author{Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.
  }

\seealso{\code{\link{comp.t}},\code{\link{comp.F}}}
  
\examples{
X <- matrix(rnorm(1000,0,0.5), nc=10)
L <- rep(0:1,c(5,5))

# genes 1-10 are differentially expressed
X[1:10,6:10]<-X[1:10,6:10]+1

f <- comp.FC(L=L)
f.X <- f(X)
}

\keyword{univar}



