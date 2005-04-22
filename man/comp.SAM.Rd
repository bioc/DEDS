\name{comp.SAM}
\alias{comp.SAM}

\title{ Computing SAM Statistics for Differential Expression}
\description{
  \code{comp.SAM} returns a function of one argument. This function has a
  environment with bindings for a series of arguments (see below). It
  accepts a microarray data matrix as its single argument, when
  evaluated, computes SAM statistics for each row of the matrix.
}

\usage{
comp.SAM(L = NULL, prob = 0.5, B = 200, stat.only = TRUE, verbose = FALSE,
deltas, s.step=0.01, alpha.step=0.01, plot.it=FALSE)
}

\arguments{
  \item{L}{A vector of integers corresponding to observation (column)
    class labels. For \eqn{k} classes, the labels must be integers
    between 0 and \eqn{k-1}.}
  \item{prob}{A numeric variable used to set the fudge factor
    \eqn{s_0} in terms of the percentile of the standard deviations of the
    genes. If set as \code{NULL}, \eqn{s_0} is calculated using the
    algorithm by Tusher et al. (see reference).}
  \item{B}{The number of permutations. For a complete enumeration,
      \code{B} should be 0 (zero) or any number not less than the total
      number of permutations.}
  \item{stat.only}{A logical variable, if \code{TRUE}, only statistics
    are calculated and returned; if \code{FALSE}, false discovery rates
    (FDRs) for a set of \eqn{\delta}{delta}(\code{deltas}) are
    calculated and returned.}
  \item{verbose}{A logical variable, if \code{TRUE}, informative mesages
    are printed during the computation process.}
  \item{deltas}{A vector of values for the threshold \eqn{\delta}{delta}; see
    Tusher et al.}
  \item{s.step}{A numeric variable specifying the size of the moving window
    acorss the gene-wise standard deviations for the selection of the
    fudge factor \eqn{s_0}.}
  \item{alpha.step}{A numeric variable specifying the increment of
    a percentile sequence between 0 and 1, from which the fudge factor
    will be chosen to minimize the coefficient of variation of statistics.  }
  \item{plot.it}{A logical variable, if \code{TRUE}, a plot between the
    coefficient of variation and the percentile sequence will be made.}
}

\details{
  The function returned by \code{comp.SAM} calculates SAM statistics for
  each row of the microarray data matrix, with bindings for \code{L},
  \code{prob}, \code{B}, \code{stat.only}, \code{verbose},
  \code{deltas}, \code{s.step}, \code{alpha.step} and \code{plot.it}. If
  \code{quantile=NULL}, the fudge factor \eqn{s_0} is calculated as the
  percentile of the gene-wise standard deviations that minimizes the
  coefficient of variation of the statistics; otherwise \eqn{s_0} is set
  as the specified percentile of standard deviations. If
  \code{stat.only=T}, only SAM statistics are returned; otherwise,
  permutation will be carried out to calculate the FDRs for a set of
  \code{deltas} specified and a FDR table will be returned in addition
  to the SAM statistics.
  
}
\value{
  \code{SAM} returns a function (F) with bindings for a series of arguments.
  When \code{stat.only=T}, the function F when evaluated returns a
  numeric vector of SAM statistics;
  When \code{stat.only=F}, the function F when evaluated returns
  a list of the following components:
  \item{geneOrder}{Order of genes in terms of differential expression;}
  \item{sam}{Sorted SAM statistics;}
  \item{fdr.table}{A matrix with columns: delta, no.significance,
    no.positive, no.negatvie, FDR(50\%), FDR(90\%).} 
}

\references{
   Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays
   applied to the ionizing radiation response, \emph{PNAS}, 98,
   5116-5121.
 }

\author{
  Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.}

\seealso{\code{\link{comp.t}}}
\examples{
X <- matrix(rnorm(1000,0,0.5), nc=10)
L <- rep(0:1,c(5,5))

# genes 1-10 are differentially expressed
X[1:10,6:10]<-X[1:10,6:10]+1

# two sample test, statistics only
sam.fun <- comp.SAM(L)
sam.X <- sam.fun(X)

# two sample test, FDR
sam.fun <- comp.SAM(L, stat.only=FALSE, delta=c(0.1, 0.2, 0.5))
sam.X <- sam.fun(X)

}
\keyword{univar}




