\name{qqnorm-methods}
\docType{methods}
\alias{qqnorm.DEDS}
\alias{qqnorm-methods}
\alias{qqnorm}

\title{Normal Q-Q Plot for DEDS Objects}
\description{
  The function \code{qqnorm.DEDS} produces normal Quantile-Quantile
  plots of statistics for \code{\link{DEDS-class}} objects. The points
  corresponding to genes with DEDS q- or adjusted p-values less than
  a user defined threshold are highlighted. 
}
\usage{
\method{qqnorm}{DEDS}(y, subset=c(1:nrow(y$stats)),
xlab = "Quantiles of standard normal", thresh = 0.05, col = palette(), pch, ...)
}

\arguments{
  \item{y}{An object of \code{\link{DEDS}}, produced by
    \code{\link{deds.stat.linkC}} or \code{\link{deds.stat}}.}
  \item{subset}{A numeric vector indicating the subset of points to be plotted.}
  \item{xlab}{A title for the x axis}
  \item{thresh}{A numeric variable specifying the threshold of
    significance in differential expression (DE) for q- or 
    p-values of the DEDS object.}
  \item{col}{A specification for the colors to be used for plotting. It
    should have a length bigger than two. The first is used for points
    with q- or adjusted p-values smaller than the specified threshold
    (group I) and the second for points with q- or adjusted p-values
    bigger than the threshold (group II).}
  \item{pch}{A specification for the type of points to be used for
    plotting.  It should have a length bigger than two. The first
    parameter is used for group I genes, and the second for group II
    genes.}
  \item{\dots}{Extra parameters for plotting.}
}
\details{
  The function \code{qqnorm.DEDS} implements a S3 method of
  \code{\link{qqnorm}} for \code{\link{DEDS}}. The \code{\link{DEDS}}
  class is a simple list-based class to store DEDS results and
  \code{qqnorm.DEDS} is used for a DEDS object that is created  by
  functions \code{\link{deds.stat}}, \code{\link{deds.stat.linkC}}. The
  list contains a "stat" component, which stores statistics from
  various statistical tests. The function \code{qqnorm.DEDS} extracts the
  "stat" component and produces a normal QQ plot for each type of
  statistics. \code{qqnorm.DEDS} as a default highlights points
  (corresponding to genes) with DEDS adjusted p- or q-values less than
  a user defined threshold.

  For DEDS objects that are created by the function
  \code{\link{deds.pval}}, the "stat" matrix consists of unadjusted
  p-values from different statistical models. For graphical display of
  these p values, the user can use \code{\link{hist.DEDS}} and
  \code{\link{pairs.DEDS}}.
}

\author{Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.}

\seealso{\code{\link{deds.stat}}, \code{\link{deds.pval}},
  \code{\link{deds.stat.linkC}}, \code{\link{hist.DEDS}},
  \code{\link{qqnorm.DEDS}}}

\examples{
X <- matrix(rnorm(1000,0,0.5), nc=10)
L <- rep(0:1,c(5,5))

# genes 1-10 are differentially expressed
X[1:10,6:10]<-X[1:10,6:10]+1
# DEDS summarizing t, fc and sam
d <- deds.stat.linkC(X, L, B=200)

# qqnorm for t, fc and sam
qqnorm(d)
# change points color
qqnorm(d, col=c(2,3))
# change points type
qqnorm(d, pch=c(1,2))
}
\keyword{hplot}



