\name{affySpikeIn}
\alias{affySpikeIn}
\alias{affySpikeIn.L}
\alias{affySpikeIn.gnames}
\alias{spikedgene}

\title{Gene expression dataset from Affymetrix Spike-in Experiments}

\usage{
data(affySpikeIn)
}

\description{
The spike-in experiment represents a portion of the data used by Affymetrix
to develop their MAS 5.0 preprocessing algorithm.  Here we utilize the
RMA (Irizarry et. al) probe level summaries.  The data features 14 
human genes spiked-in at a series of 14 known concentrations 
(\eqn{0, 2^{-2}, 2^{-1}, \ldots, 2^{10}} pM) according to 
a Latin square design among 12612 null genes.  The data matrix
\code{affySpikeIn} represents the two array groups among the 14 array
groups that contain 12 replicates. Further details are available 
at \url{http://www.affymetrix.com/analysis/download_center2.affx}. }

\value{
  \item{affySpikeIn}{matrix of gene expression levels measurements,
    rows correspond to genes (12626 genes) and columns to 24 samples.}
  \item{affySpikeIn.L}{numeric vector indicating the sample class, 12 (code 0)
    vs. 12 (code 1).}
  \item{affySpikeIn.gnames}{character vecter containing the names of the
    12626 genes.}
  \item{spikedgene}{numeric vector given the location of the 14 spiked genes.}
}

\references{
R. A. Irizarry, B. M. Bolstad, F. Collin, L. Cope, B. Hobbs and
T. P. Speed (2003) Summeries of affymetrix genechip probe level data.
\emph{Nucleic Acide Research}, 31:e15. }

\keyword{datasets}

