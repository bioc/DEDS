\name{comp.B}
\alias{comp.B}
\title{Computing B-statistics for Differential Expression}
\description{
  \code{comp.B} returns a function of one argument with bindings for
  \code{L} and \code{proportion}. This function accepts a microarray
  data matrix as its single argment, when evaluated, computes lod-odds
  of differential expression by emprical Bayes shrinkage of the standard
  error toward a common value. The lod-odds are sometimes called B
  statistics. 
}
\usage{
comp.B(L = NULL, proportion = 0.01)
}

\arguments{
  \item{L}{A vector of integers corresponding to observation (column)
    class labels. For \eqn{k} classes, the labels must be integers
    between 0 and \eqn{k-1}.}
  \item{proportion}{A numeric variable specifying the proportion of
    differential expression.}
}
\details{
  The function returned by \code{comp.B} calculates B statistics for
  each row of the microarray data matrix, with bindings for \code{L} and
  \code{proportion}. It interfaces to a C function. \code{\link{comp.stat}}
  is another function that wrapps around the same C function that could
  be used for computing B statistics (see examples below).
}
\value{
  \code{comp.B} returns a function (F) with the bindings for
  \code{L} and \code{proportion} . The function F when supplied with
  a microarray data matrix and evaluated will return a numeric vector of
  B  statistics for each row of the matrix.
}
\references{
  Lonnstedt, I. and Speed, T. P. (2002). Replicated microarray
  data. \emph{Statistica Sinica} 12, 31-46.
  
   Smyth, G. K. (2003). Linear models and empirical Bayes methods for
  assessing differential expression in microarray
  experiments. http://www.statsci.org/smyth/pubs/ebayes.pdf
}

\author{
  Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jeany@maths.usyd.edu.au}.
}

\seealso{\code{\link{comp.modt}},\code{\link{comp.stat}}.}
\examples{
X <- matrix(rnorm(1000,0,0.5), nc=10)
L <- rep(0:1,c(5,5))

# genes 1-10 are differentially expressed
X[1:10,6:10]<-X[1:10,6:10]+1

# compute B statistics, proportion set as 0.01
B.fun <- comp.B(L)
B.X <- B.fun(X)

# compute B statistics, proportion set as 0.1
B.fun <- comp.B(L, proportion=0.1)
B.X <- B.fun(X)

# Another way of computing B statistics
B.X<- comp.stat(X, L, "B")
}

\keyword{univar}





