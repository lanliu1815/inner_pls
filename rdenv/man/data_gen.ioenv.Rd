% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/For Simulations.R
\name{data_gen.ioenv}
\alias{data_gen.ioenv}
\title{For Inner Envelope Studying.}
\usage{
data_gen.ioenv(
  SAM = F,
  d,
  u2 = p - d + 1,
  u = r - d,
  p,
  r,
  n,
  x,
  dia,
  upper.eta = 10
)
}
\arguments{
\item{SAM}{Whether to use sample version.}

\item{u2, u}{u2 columns of Gamma0 is contained within beta.part. u is the dimension
of the envelope-structured data.}

\item{upper.eta}{The upper bound for the distribution of eta, which is a
uniform one, with the lower bound being 0.}
}
\description{
Generate two sets of inner envelope-structured data, the envelope structure of one is
 the full space and the other is not, and one set of envelope-structured data.
}
