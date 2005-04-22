\name{hist-methods}
\docType{methods}
\alias{hist.DEDS}
\alias{hist-methods}
\alias{hist}

\title{Histogram for DEDS Objects}
\description{
  The function \code{hist.DEDS} produces histograms of unadjusted
  p-values for \code{\link{DEDS-class}} objects. 
}
\usage{
\method{hist}{DEDS}(x, subset=c(1:nrow(x$stats)), ...)
}

\arguments{
  \item{x}{An object of \code{\link{DEDS}}, produced by
    \code{\link{deds.pval}}.}
  \item{subset}{A numeric vector indicating the subset of points to be plotted.}
  \item{\dots}{Further graphical parameters, for example, \code{"col"},
    \code{"border"}, \code{"main"}, \code{"nclass"}.} 
}
\details{
  The function \code{hist.DEDS} implements a S3 method of
  \code{\link{hist}} for \code{\link{DEDS}}. The \code{\link{DEDS}}
  class is a simple list-based class to store DEDS results and
  \code{hist.DEDS} is used for a DEDS object that is created  by
  the function \code{\link{deds.pval}}. The list contains a "stat"
  component, which stores unadjusted p-values from various statistical
  models. The function \code{hist.DEDS} extracts the "stat" component
  and produces a histgram of the unadjusted p-values for each model.

   For DEDS objects that are created by functions
   \code{\link{deds.stat}} and \code{\link{deds.stat.linkC}}, the "stat"
   matrix consists of different types of statistics. For graphical
   dispaly of these statistics, the user can use \code{\link{qqnorm.DEDS}} and
  \code{\link{pairs.DEDS}}.
}

\author{Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.}

\seealso{\code{\link{deds.stat}}, \code{\link{deds.pval}},
  \code{\link{deds.stat.linkC}}, \code{\link{pairs.DEDS}},
  \code{\link{qqnorm.DEDS}}}

\examples{


}
\keyword{hplot}




