% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Preliminaries.R
\name{mpinv}
\alias{mpinv}
\title{MP Inverse}
\usage{
mpinv(X, rank = NULL, tol = sqrt(.Machine$double.eps))
}
\arguments{
\item{rank}{Specify the rank of the matrix. If null, then rank will be selected,
just as the original \code{ginv} function does.}
}
\description{
The Moore-Penrose inverse. Mainly copied from \code{\link[MASS]{ginv}}.
}
