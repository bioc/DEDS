\name{comp.stat}
\alias{comp.stat}

\title{Computing Test Statistics for Differential Expression}
\description{
  This function computes test statistics, e.g., t-statitics, F-statistics, 
  SAM, fold changes, moderated t or F statistics, B statistics, 
  for each row of a microarray data matrix.
}

\usage{
comp.stat(X, L, test = c("t", "fc", "sam", "f", "modt", "modf", "B"), extra = NULL)
}
\arguments{
  \item{X}{A  matrix, with \eqn{m} rows corresponding to variables
    (hypotheses) and \eqn{n} columns to observations. In the case of gene 
    expression data, rows correspond to genes and columns to mRNA
    samples. The data can be read using \code{\link{read.table}}.}
  \item{L}{A vector of integers corresponding to observation (column)
    class labels. For \eqn{k} classes, the labels must be integers
    between 0 and \eqn{k-1}.}
  \item{test}{A character string specifying the statistic to be
    used to test the null hypothesis of no association between the
    variables and the class labels.\cr
    \tabular{ll}{
      \code{test="t"}: \tab t-statistics; \cr
      \code{test="f"}: \tab F-statistics;\cr
      \code{test="fc"}:\tab fold changes;\cr
      \code{test="sam"}: \tab SAM-statistics; \cr
      \code{test="modt"}: \tab moderated t-statistics;\cr
      \code{test="modf"}: \tab moderated F-statistics; \cr
      \code{test="B"}: \tab B-statistics.}}
  \item{extra}{Extra parameter needed for the test specified; see
    \code{\link{deds.genExtra}}.}
}

\details{
  The function \code{comp.stat} interfaces to a C function and
  computes various statistics for differential expression in the C
  environment and therefore faster than functions in R. However,
  functions in R that are implemented in the DEDS packages may have
  more flexibility in terms of specifications of arguments. Below is a
  table the details \code{comp.stat} and its equivalent R functions
  in the DEDS package. Note that all the R functions listed in the 2nd
  column of the table below return a function with bindings for a series
  of arguments which accept the microarray data matrix as its single
  argument and compute accordingly statistics. \cr
  
  \tabular{lll}{
    Interface to C \tab R functions \tab Statistics \cr
    deds.stat(X, L, test="t") \tab tTest(L=NULL, mu=0, var.equal=FALSE)
    \tab t statistics \cr
    deds.stat(X, L, test="fc") \tab FC(L=NULL, is.log=TRUE, FUN=mean)
    \tab fold change \cr
    deds.stat(X, L, test="sam") \tab Sam(L=NULL, prob=0.5, B=200,
    stat.only=TRUE, verbose=FALSE, deltas, s.step=0.01, alpha.step=0.01,
    plot.it=FALSE) \tab SAM statistics \cr
    deds.stat(X, L, test="f") \tab fTest(L=NULL) \tab F statistics \cr
    deds.stat(X, L, test="modt") \tab tmodTest(L=NULL) \tab moderated t
    statistics \cr
    deds.stat(X, L, test="modf") \tab fmodTest(L=NULL) \tab moderated F
    statistics \cr
    deds.stat(X, L, test="B") \tab BTest(L=NULL, proportion=0.01) \tab B
    statistics }
    
}
\value{
  A vector of test statistics for each row of the matrix.
}
\references{
  For references on B-statistics and moderated t and F statistics:
  
  Lönnstedt, I. and Speed, T. P. (2002). Replicated microarray
  data. \emph{Statistica Sinica} \bold{12}, 31-46.
  
  Smyth, G. K. (2003). Linear models and empirical Bayes methods for
  assessing differential expression in microarray
  experiments. http://www.statsci.org/smyth/pubs/ebayes.pdf
}

\author{Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.}

\seealso{\code{\link{deds.genExtra}},
  for B statistics: \code{\link[limma]{lm.series}} and
  \code{\link[limma]{ebayes}}
}
\examples{
X <- matrix(rnorm(1000,0,0.5), nc=10)
L <- rep(0:1,c(5,5))

# genes 1-10 are differentially expressed
X[1:10,6:10]<-X[1:10,6:10]+1

# t statistics
tstat <- comp.stat(X, L, test="t")

# SAM, fudge factor set as the median of pooled genewise standard deviations
samstat <- comp.stat(X, L, test="sam")
# SAM, fudge factor set as the 90\% of pooled genewise standard deviations
samstat <- comp.stat(X, L, test="sam", extra=c(0.9))

# moderated t
modtstat <- comp.stat(X, L, test="modt")

# B, proportion of differentially expressed genes is set at default, 1\%
Bstat <- comp.stat(X, L, test="B")
# B, proportion of differentially expressed genes is set at  10\%
Bstat <- comp.stat(X, L, test="B", extra=c(0.1))
}

\keyword{univar}




