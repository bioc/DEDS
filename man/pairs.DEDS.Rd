\name{pairs-methods}
\docType{methods}
\alias{pairs.DEDS}
\alias{pairs}
\alias{pairs-methods}

\title{Pairs Plot for DEDS Objects}
\description{
  The function \code{pairs-DEDS} produces pairs plots of statistics or p
   values for \code{\link{DEDS-class}} objects.
}
\usage{
\method{pairs}{DEDS}(x, subset=c(1:nrow(x$stats)), labels =
colnames(x$stats[,-1]), logit = FALSE,
diagonal = c("qqnorm", "boxplot", "density", "histogram", "none"),
lower = c("cor", "none"), groups.by.deds = TRUE, thresh = 0.05, reg.line
= NULL, smooth = FALSE, line.by.group = FALSE, diag.by.group = TRUE, lower.by.group =
FALSE, col = palette(), pch = 1:n.groups, lwd = 1, legend.plot =
length(levels(groups)) > 1, ...)
}

\arguments{
  \item{x}{An object of \code{\link{DEDS}}.}
  \item{subset}{A numeric vector indicating the subset of points to be plotted.}
  \item{labels}{A character vector specifying the names of the variables.}
  \item{logit}{A logical variable, if \code{TRUE} the variables are
    logged, useful when plotting p values.}
  \item{diagonal}{A character string specifying the type of plot to be
    applied in the diagonal panels. \cr
    \tabular{ll}{
      \code{diagonal="qqnorm"}: \tab  \code{\link{qqnorm}} on the diagonal\cr
      \code{diagonal="boxplot"}: \tab \code{\link{boxplot}} on the diagonal\cr
      \code{diagonal="density"}: \tab \code{\link{density}} on the diagonal\cr
      \code{diagonal="histogram"}: \tab \code{\link{hist}} on the diagonal\cr
      \code{diagonal="none"}: \tab no special plot will be applied on
      the diagonal
    }
    }      
  \item{lower}{A character string specifying the function to be
    applied in the lower panels. \cr
    \tabular{ll}{
      \code{lower="cor"}: \tab  absolute correlation will be put on the
      lower panel\cr
      \code{none="cor"}; \tab no special function will be applied
      }
    }
  \item{groups.by.deds}{A logcial variable, if \code{TRUE}, points will
    be separated into groups according to their magnitude of q- or p-values
    by DEDS. }
  \item{thresh}{A numeric variable, if \code{thresh}<1, it specifies
    the threshold of significance in differential expression (DE) for q- or 
    p-values of the DEDS object; default is set at 0.05. If 
    \code{thresh}>1, it specifies the number of top DE genes to be
    highlighted. }
  \item{reg.line}{A function name specifying the type of regression line
    to be plotted in the scatter plots. If \code{reg.line=lm},
    linear regression line will be plotted; If \code{reg.line=NULL}, no
    regression line will be plotted in the scatter plot.}
  \item{smooth}{A logical variable specifying if smooth regression lines
  will be plotted in the scatter plots. If \code{smooth=TRUE}, a
  \code{\link{lowess}} line will be applied.}
  \item{line.by.group}{A logical variable specifying if the regression
    lines should be applied within groups.}
  \item{diag.by.group}{A logical variable specifying if the plot in the diagonal
    panels would be applied groupwise.}
  \item{lower.by.group}{A logical variable, if
    \code{lower.by.group=TRUE} and \code{lower="cor"}, correlation
    coefficients will be calculated and printed separated according to
    groups in the lower panels.}
  \item{col}{A specification for the colors to be used for plotting
    different groups, see \code{\link{par}}.}
  \item{pch}{A specification for the type of points to be used for
    plotting different groups, see \code{\link{par}}.}
  \item{lwd}{A specification for the width of lines to be used if lines
    are plotted; see \code{\link{par}}.}
  \item{legend.plot}{A logical variable specifying if the legend will be
    plotted.}
  \item{\ldots}{Extra parameters for plotting.}
}
\details{
  The function \code{pairs.DEDS} implements a S3 method of
  \code{\link{pairs}} for \code{\link{DEDS}}. The \code{\link{DEDS}}
  class is a simple list-based class to store DEDS results and it is
  usually created  by functions \code{\link{deds.pval}},
  \code{\link{deds.stat}}, \code{\link{deds.stat.linkC}}. The list
  contains a "stat" component, which stores statistics or p values from
  various statistical tests. The function \code{pairs.DEDS} extracts the
  "stat" component and produces a matrix of scatterplot.

  \code{pairs.DEDS} as a default highlights points (corresponding to
  genes) with adjusted p- or q-values less than a user defined
  threshold. The user can select among a series of options a plot for
  the diagonal panel; as a default, it produces a \code{\link{qqnorm}}
  for each column in the "stat" matrix. Both the diagonal and lower
  panels can be stratified by specifying the \code{diag.by.group} or
  \code{lower.by.group} arguments.
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
# DEDS
d <- deds.stat.linkC(X, L, B=200)

# pairs plot
pairs(d)
# plot regression line
pairs(d, reg.line=lm, lwd=2)
# histogram in the diagonal panel
pairs(d, diagonal="hist")
# boxplot on the diagonal panel and stratified
pairs(d, diagonal="boxplot", diag.by.group=TRUE)
}

\keyword{hplot}



