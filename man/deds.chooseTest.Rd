\name{deds.chooseTest}
\alias{deds.chooseTest}

\title{Selection of the Most Common Statistics for Differential Expression}
\description{
  This function selects a set of functions of common statistics for differential
  expression in microarray data analysis, given specific observation
  class labels. As a default, t-statistics, fold change and SAM are selected.
}
\usage{
deds.chooseTest(L = NULL, tests = c("t", "sam", "fc"))
}

\arguments{
  \item{L}{A vector of integers corresponding to observation (column)
    class labels. For \eqn{k} classes, the labels must be integers
    between 0 and \eqn{k-1}.}
  \item{tests}{A character vector specifying the statistics to be
    used to test the null hypothesis of no association between the
    variables and the class labels. For DEDS, there should be more than
    one statistic chosen from the following:\cr
    \tabular{ll}{
      "t": \tab t-statistics; \cr
      "f": \tab F-statistics;\cr
      "fc":\tab fold changes;\cr
      "sam": \tab SAM-statistics; \cr
      "modt": \tab moderated t-statistics;\cr
      "modf": \tab moderated F-statistics; \cr
      "B": \tab B-statistics.}}
}
\details{
  \code{deds.chooseTest} can be used together with the function
  \code{deds.stat}. The user specifies the types of statistics needed for
  subsequent DEDS analysis by the argument \code{tests} and the function
  returns accordingly a list the statistics function, which could be
  used for input \code{testfun} in the function \code{deds.stat}.
}
\value{
  A list of statistics functions specified by the user which could be
  used for input in the function \code{deds.stat}.
}

\author{
  Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.
}

\seealso{\code{\link{comp.t}}, \code{\link{comp.FC}}, \code{\link{comp.SAM}}}
\examples{
X <- matrix(rnorm(1000,0,0.5), nc=10)
L <- rep(0:1,c(5,5))

# genes 1-10 are differentially expressed
X[1:10,6:10]<-X[1:10,6:10]+1

# as a default, chooses t, fc and sam
funcs <- deds.chooseTest(L)
deds.X <- deds.stat(X, L, testfun=funcs)

# chooses F statistic, SAM statistic, and moderated F statistic
L <- rep(0:2, c(3,3,4))
funcs <- deds.chooseTest(L, tests=c("f", "sam", "modf"))

}
\keyword{univar}





