\name{DEDS-class}
\docType{class}
\alias{DEDS-class}
\alias{DEDS}
\title{DEDS Result List - class}

\description{
  A simple list-based class to store DEDS results.\code{DEDS} objects
  are usually created  by functions \code{\link{deds.pval}},
  \code{\link{deds.stat}} and \code{\link{deds.stat.linkC}}.
}
 
\section{Slots/List Components}{
  \code{DEDS} objects can be created by \code{new("DEDS",deds)} where
  \code{deds} is a list. This class contains no slots, but objects
  should contain the following list components:
  \tabular{ll}{
    \code{E}: \tab A numeric vector of the most extreme point in the direction
      of differential expression. \cr
    \code{p}: \tab A numeric vector of q- or adjusted p-values.  \cr
    \code{geneOrder}:  \tab An integer vector giving the index of the top genes
      in terms of differential expression. \cr
    \code{stats}: \tab A matrix of p values or statistics. \cr
    \code{options}: \tab A character vector of options used in the test.
    }
}

\section{Methods}{
  This class inherits directly from class \code{list}, so any operation
  appropriate for lists will work on objects of this class. In addition,
  Other functions which operate on \code{DEDS} objects include
  \code{\link{pairs}} and \code{\link{hist}}.
}

\author{Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
    Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}.}

\seealso{
  \code{\link{deds.pval}}, \code{\link{deds.stat}},
  \code{\link{deds.stat.linkC}}
}

\keyword{classes}
\keyword{data}




