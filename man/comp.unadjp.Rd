\name{comp.unadjp}
\alias{comp.unadjp}

\title{Computing permutation based unadjusted p values for each row of a matrix}
\description{
  This function computes permuation based unadjusted p values for a
  selected test statistic, e.g., one- or two-sample
  t-statitics, F-statistics, SAM, Fold change, for each row of a matrix.
}

\usage{
comp.unadjp(X, L, B = 1000, test = c("t", "fc", "sam", "f"), tail = c("abs", "lower", "higher"), extra = NULL)
}

\arguments{
  \item{X}{A  matrix, with \eqn{m} rows corresponding to variables
    (hypotheses) and\eqn{n} columns corresponding to observations.
    In the case of gene expression data, rows correspond to genes and
    columns to mRNA samples. The data can be read using \code{\link{read.table}}.
    }
  \item{L}{A vector of integers corresponding to observation (column)
    class labels. For \eqn{k} classes, the labels must be integers
    between 0 and \eqn{k-1}.
  }
  \item{B}{The number of permutations. For a complete enumeration,
    \code{B} should be 0 (zero) or any number not less than the total
    number of permutations.
  }
  \item{test}{A character string specifying the statistic to be
    used to test the null hypothesis of no association between the
    variables and the class labels.\cr
    If \code{test="t"}, for one-class, the tests are based on one-sample
    t-statistics; for two-class, the tests are based on two-sample t-statistics
    (unequal variances).  \cr
    If \code{test="f"}, the tests are based on F-statistics.\cr
    If \code{test="fc"}, the tests are based on fold changes among classes.\cr
    If \code{test="sam"}, the tests are based on SAM-statistics.
    }
  \item{tail}{A character string specifying the type of rejection region.\cr
    If \code{side="abs"}, two-tailed tests, the null hypothesis is rejected for large absolute values of the test statistic.\cr
    If \code{side="higher"}, one-tailed tests, the null hypothesis is rejected for large values of the test statistic.\cr
    If \code{side="lower"}, one-tailed tests,  the null hypothesis is rejected for small values of the test statistic.
   }
  \item{extra}{Extra parameter need for the test specified; see
    \code{\link{deds.genExtra}}.}
}
\details{
  The function \code{comp.unadjp} computes unadjusted \eqn{p} values using
  a permutation scheme.
}
\value{
  A vector of unadjusted \eqn{p} values for each row of the matrix.
  }
  
\author{Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.}

\seealso{\code{\link{deds.genExtra}}, \code{\link{comp.stat}}}
 
\examples{

X <- matrix(rnorm(1000,0,0.5), nc=10)
L <- rep(0:1,c(5,5))

# genes 1-10 are differentially expressed
X[1:10,6:10]<-X[1:10,6:10]+1

# t statistics
unadjp.t <- comp.unadjp(X, L, test="t")
}

\keyword{univar}





