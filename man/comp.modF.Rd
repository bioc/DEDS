\name{comp.modF}
\alias{comp.modF}

\title{Computing Moderated t-statistics for Differential Expression}
\description{
  \code{comp.modF} returns a function of one argument with bindings for
  \code{L}. The function accepts a microarray data matrix as its single
  argument, when evaluated, computes moderated F-statistics by emprical
  Bayes shrinkage of the standard error toward a common value.
}

\usage{
comp.modF(L = NULL)
}

\arguments{
  \item{L}{A vector of integers corresponding to observation (column)
    class labels. For \eqn{k} classes, the labels must be integers
    between 0 and \eqn{k-1}.}
}
\details{
  The function returned by \code{comp.modF} computes moderated F
  statistics for the assessment of differential expression. It
  interfaces to a C function. \code{\link{comp.stat}} is another
  function that wrapps around the C function that could be used for
  computing moderated F statistics. For details of moderated statistics,
  see Smyth (2003). 
}
\value{
   \code{comp.modF} returns a function (F) with the bindings for
   \code{L}. The function F when supplied with a microarray data matrix
   and evaluated will return a numeric vector of moderated F statistics
   for each row of the matrix.
 }

\references{
  Lönnstedt, I. and Speed, T. P. (2002). Replicated microarray
  data. \emph{Statistica Sinica} \bold{12}, 31-46.
  
  Smyth, G. K. (2003). Linear models and empirical Bayes methods for
  assessing differential expression in microarray
  experiments. http://www.statsci.org/smyth/pubs/ebayes.pdf
}

\author{
  Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.
}

\seealso{\code{\link{comp.FC}}, \code{\link{comp.modt}}, \code{\link{comp.stat}}}

\examples{
X <- matrix(rnorm(1000,0,0.5), nc=10)
L <- rep(0:1,c(5,5))

# genes 1-10 are differentially expressed
X[1:10,6:10]<-X[1:10,6:10]+1
fmod <- comp.modF(L)
fmod.X <- fmod(X)

# Another way of computing moderated F statistics
fmod.X <- comp.stat(X, L, "modf")
}

\keyword{univar}



