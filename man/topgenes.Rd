\name{topgenes}
\alias{topgenes}

\title{Table of Top Genes from DEDS}
\description{
  \code{topgenes} prints a table of top-ranked genes by
  \code{\link{DEDS}}.}

\usage{
topgenes(obj, number = 10, genelist = NULL, sort.by = c("deds", colnames(obj$stats[,-1])), tail = c("abs", "lower", "higher"))
}

\arguments{
  \item{obj}{An object of \code{\link{DEDS}}}
  \item{number}{A numeric variable specifying the number of top genes to
    be printed out.}
  \item{genelist}{A dataframe or a vector containing gene names.}
  \item{sort.by}{A character string specifying the name of the
    statistic to sort genes by. The default uses the DEDS result, the
    user can also choose from the names of the statistics (or unadjusted
    p values) that DEDS is used to summarize.}
  \item{tail}{A character string specifying the type of rejection
    region.\cr
    If \code{side="abs"}, two-tailed tests, genes are ranked by their
    absolute values. \cr
    If \code{side="higher"}, one-tailed tests, genes are ranked
    decreasingly. \cr
    If \code{side="lower"}, one-tailed tests, genes are ranked increasingly.}
}
\details{
  The function \code{topgenes} accepts a \code{\link{DEDS}} object as
  the first argument. The \code{\link{DEDS}} class is a simple
  list-based class to store DEDS results. The list contains a "stat"
  component, which stores statistics or unadjusted p-values from
  various statistical tests. The function \code{topgenes.DEDS} extracts the
  "stat" component and prints out the top genes according to the user
  defined criterion -- usually by DEDS or by a single statistical
  measure that DEDS summarizes. 
}

\value{
  A dataframe with rows for selected genes, and columns for the
  "stat" matrix and q- or adjusted p-values from DEDS.
}

\author{Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.}

\seealso{\code{\link{deds.stat.linkC}}, \code{\link{deds.pval}}}
\examples{
X <- matrix(rnorm(1000,0,0.5), nc=10)
L <- rep(0:1,c(5,5))

# genes 1-10 are differentially expressed
X[1:10,6:10]<-X[1:10,6:10]+1

}

\keyword{htest}



