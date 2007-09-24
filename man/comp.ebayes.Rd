\name{comp.ebayes}
\alias{comp.ebayes}

\title{Computing Empirical Bayes Statistics for Differential Expression}

\description{
  \code{comp.ebayes} returns a function of one argument with bindings for
  \code{L} and \code{proportion}. This function accepts a microarray
  data matrix as its single argment, when evaluated, computes lod-odds
  (B statistics) and moderated t statistics of differential expression
  by emprical Bayes shrinkage of the standard error toward a common value. 
}

\usage{
comp.ebayes(L = NULL, proportion = 0.01)
}

\arguments{
  \item{L}{A vector of integers corresponding to observation (column)
    class labels. For \eqn{k} classes, the labels must be integers
    between 0 and \eqn{k-1}.}
  \item{proportion}{A numeric variable specifying the proportion of
    differential expression.}
}
\details{
   The function returned by \code{comp.ebayes} calculates B statistics
   and moderated t statistics for each row of the microarray data
   matrix, with bindings for \code{L} and \code{proportion}. It
   interfaces to a C function. 
}
\value{
  \code{comp.ebayes} returns a function (F) with the bindings for
  \code{L} and \code{proportion} . The function F when supplied with
  a microarray data matrix and evaluated will return a matrix of two
  columns:
  \item{t}{Moderated t statistics}
  \item{B}{B statistics (log-odds) of differential expression}}
}

\references{
  Lonnstedt, I. and Speed, T. P. (2002). Replicated microarray
  data. \emph{Statistica Sinica} \bold{12}, 31-46.
  
  Smyth, G. K. (2003). Linear models and empirical Bayes methods for
  assessing differential expression in microarray
  experiments. http://www.statsci.org/smyth/pubs/ebayes.pdf
}

\author{
  Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.
}

\seealso{\code{\link{comp.modt}},\code{\link{comp.B}}.}
\examples{
X <- matrix(rnorm(1000,0,0.5), nc=10)
L <- rep(0:1,c(5,5))

# genes 1-10 are differentially expressed
X[1:10,6:10]<-X[1:10,6:10]+1

# compute B and moderated t statistics, proportion set as 0.01
ebayes.fun <- comp.ebayes(L)
ebayes.X <- ebayes.fun(X)

# compute B and moderated t statistics, proportion set as 0.1
ebayes.fun <- comp.ebayes(L, proportion=0.1)
ebayes.X <- ebayes.fun(X)

}
\keyword{htest}







