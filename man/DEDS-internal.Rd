\name{DEDS-internal}
\alias{deds.checkB}
\alias{deds.checkothers}
\alias{deds.checkclasslabel}
\alias{deds.checkX}
\alias{test.checkX}
\alias{deds.calcFDR}
\alias{deds.calcAdjP}
\alias{deds.next.sample}
\alias{type2test}
\alias{sam.fdr}
\alias{sam.multiclass.func}
\alias{sam.twoclass.func}
\alias{sam.oneclass.func}
\alias{sam.s0}
\alias{euclidean}

\title{Internal DEDS functions}
\description{
  Internal DEDS functions}
\usage{
deds.checkB(classlabel, B, verbose = TRUE)
deds.checkothers(tail = "abs", distance = "weuclid", adj = "fdr")
deds.checkclasslabel(classlabel, tests)
deds.checkX(X, classlabel, tests, nsig = nrow(X), B)
test.checkX(X, classlabel, test)
deds.calcFDR(bD, D, R, K = length(D))
deds.calcAdjP(bD, D, R, K = length(D))
deds.next.sample(L)
type2test(tail = c("abs", "higher", "lower"), distance = c("weuclid",
"euclid"), adj = c("fdr", "adjp"), L)
sam.fdr(order.t, ordertB, deltas)
sam.multiclass.func(X, L, prob = 0.5, B = 200, stat.only = TRUE,
verbose = FALSE, s.step=0.01, alpha.step=0.01, plot.it=FALSE)
sam.twoclass.func(X, L, prob = 0.5, B = 200, stat.only = TRUE,
verbose = FALSE, s.step=0.01, alpha.step=0.01, plot.it=FALSE)
sam.oneclass.func(X, L, prob = 0.5, B = 200, stat.only = TRUE,
verbose = FALSE, s.step=0.01, alpha.step=0.01, plot.it=FALSE)
sam.s0(r, s, s.step = 0.01, alpha.step = 0.01, plot.it = FALSE)
euclidean(X, center, wval)
}

\details{
  These are not to be called directly by the user.
}
\keyword{internal}

