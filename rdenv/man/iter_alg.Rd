% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Estimation for IPE.R
\name{iter_alg}
\alias{iter_alg}
\title{The iteration algorithm. Obsolete now.}
\usage{
iter_alg(d, r, Gam1, Gam0, X, Y, ftol = 1e-05, beta.t, PRINT)
}
\arguments{
\item{ftol}{Should be less than 10.}
}
\value{
\item{beta.iter}{}
\item{eta1}
\item{eta2}
\item{C}
}
\description{
This algorithm is no good and is deprecated.
}
\details{
Note that \eqn{\beta}(input and output) is the non-transposed version: same to response envelopes.
}
