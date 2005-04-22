\name{ApoA1}
\alias{ApoA1}
\alias{ApoA1.L}

\title{Gene expression dataset from the ApoA1 Experiment}

\usage{
data(ApoA1)
}

\description{
  Gene expression data (6384 genes and 16 samples) from a
  study of a mouse model with very low HDL cholesterol levels
  described in Dudoit et al. (2002). Pre-processing
  was done as described in Dudoit et al. (2002). }

\value{
  \item{ApoA1}{matrix of gene expression levels measurements,
    rows correspond to genes (6384 genes) and columns to 16 samples.}
  \item{ApoA1.L}{numeric vector indicating the sample class, 8 (code 0)
    vs. 8 (code 1).}
}

\references{
  S. Dudoit, Y. H. Yang, T. P. Speed and M. J. Callow (2002) Statistical
  methods for identifying differentially expressed genes in replicated
  cDNA microarray experiments. \emph{Statistica Sinica}, Vol. 12, No. 1,
  pp. 111-139. }

\keyword{datasets}

