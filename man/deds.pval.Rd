\name{deds.pval}
\alias{deds.pval}

\title{Differential Expression via Distance Summary of p Values from
  Multiple Models}
\description{
  \code{deds.pval} integrates different \eqn{p} values of differential
  expression (DE) to rank and select a set of DE genes.
}
\usage{
deds.pval(X, E = rep(0, ncol(X)), adj = c("fdr", "adjp"), B = 200)
}

\arguments{
  \item{X}{A matrix, with \eqn{m} rows corresponding to variables
    (hypotheses) and \eqn{n} columns corresponding to \eqn{p} values from
    different statistical models.}
  \item{E}{A numeric vector indicating the location of the most extreme
    \eqn{p} values in the direction of differential expression.}
  \item{adj}{A character string specifying the type of multiple testing
    adjustment. \cr
    If \code{adj="fdr"}, False Discovery Rate is controled and q values
    are returned. \cr
    If \code{adj="adjp"}, ajusted \eqn{p} values that controls family wise
      type I error rate is returned.}
  \item{B}{The number of permutations. For a complete enumeration,
      \code{B} should be 0 (zero) or any number not less than the total
      number of permutations.}
}
\details{
  \code{deds.pval} summarizes \eqn{p} values from multiple statistical models
  for the evidence of DE. The DEDS methodology treats each gene as
  a point corresponding to a gene's vector of DE measures. An "extreme
  origin" is defined as the point that indicate DE, typically a vector
  of zero \eqn{p} values. The  distance from all points to the extreme is
  computed and the ranking of a gene for DE is determined by the
  closeness of the gene to the extreme. To determine a cutoff for
  declaration of DE, null referent distributions are generated using an
  approach similar to the gap statistic (see Reference below). DEDS can also summarize
  different statistics, see \code{\link{deds.stat}} and
  \code{\link{deds.stat.linkC}}.
}
\value{
  An object of class \code{\link{DEDS}}. See \code{\link{DEDS-class}}.
}

\references{
  Tibshirani, R., Walther G., and Hastie T. (2000). Estimating the
  number of clusters in a dataset via the gap statistic. Department of
  Statistics,  Stanford University,
  http://www-stat.stanford.edu/~tibs/ftp/gap.ps

  Yang, Y.H., Xiao, Y. and Segal M.R.: Selecting differentially expressed
  genes from microarray experiment by sets of
  statistics. \emph{Bioinformatics} 2005 21:1084-1093.
}

\author{Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.
}

\seealso{\code{\link{deds.stat}}, \code{\link{deds.stat.linkC}}.}
\examples{


}
\keyword{htest}







