\name{deds.genExtra}
\alias{deds.genExtra}

\title{Generating Extra Parameters for Test Statistics Functions for
  Differential Expression}
\description{
  \code{deds.genExtra} is used to pass in extra arguments for
  \code{\link{comp.stat}} and \code{\link{deds.stat.linkC}},
  which computes various test statistics for differential expression
  in microarray data. 
}
\usage{
deds.genExtra(classlabel, tests)
}

\arguments{
  \item{classlabel}{A vector of integers corresponding to observation (column)
    class labels. For \eqn{k} classes, the labels must be integers
    between 0 and \eqn{k-1}.}
  \item{tests}{A vector of character string specifying the statistics to be
    used to test differential expression. The character string could be
    any of the followings: \cr
    \tabular{ll}{
      \code{test="t"}: \tab one- or two-sample t-statistics; \cr
      \code{test="f"}: \tab F-statistics;\cr
      \code{test="fc"}: \tab fold changes among classes;\cr
      \code{test="sam"}: \tab SAM-statistics; \cr
      \code{test="modt"}: \tab moderated t-statistics; \cr
      \code{test="modf"}: \tab moderated F-statistics; \cr
      \code{test="B"}: \tab B-statistics.}}
}
\details{
  Given the names of the test statistics, \code{deds.genExtra} generates
  extra parameters needed to be passed in the functions
  \code{comp.stat} and \code{deds.stat.linkC} for the assessment of
  differential expression. Both functions are interfaces to C functions.
  \code{deds.genExtra} generates default parameters as follows: \cr
  If \code{test="t" or "f", "fc", "modt", "modf"}, the extra parameter
  needed is the number of classes; \cr
  If \code{test="sam"}, the extra parameter needed is the percentile of
  within-gene standard deviations that the fudge factor \eqn{s_0} will
  be set at and the default is 0.5; \cr
  If \code{test="B"}, the extra parameter needed is the percentage of
  alternative hypotheses (differential expression) and the default is
  set at 0.01.
}

\value{
  A numeric vector, the length of which is determined by the length of
  the names of the test statistics for the argument \code{test}.
}

\author{Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.}

\seealso{\code{\link{comp.stat}}, \code{\link{deds.stat.linkC}}}
\examples{
## two sample test
L <- rep(0:1, c(5,5))
extras <- deds.genExtra(L, c("t","sam", "B"))
## extras will be c(2, 0.5, 0.01)
}
\keyword{manip}





