% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/For Simulations.R
\name{data_gen.env}
\alias{data_gen.env}
\title{Data Generation for Envelope estimation.}
\usage{
data_gen.env(SAM = F, u, p, r, n, x, dia, upper.eta = 10)
}
\arguments{
\item{SAM}{Whether to generate data of sample version.}

\item{u}{u<r.}

\item{dia}{The first \code{u} elements of \code{dia} are taken as the
eigenvalues corresponding to the envelope subspace.}
}
\value{
\code{list(Sigma=Sigma, U=U, beta=beta, Gamma1=Gamma1,
Gamma0=Gamma0, eta=eta)} for population version.

\code{return(list(beta=beta, betaOLS=betaOLS, M=M, U=U, Y=Y, X=X, epsilon=eps,
 Gamma1=Gamma1, Gamma0=Gamma0, Sigma=Sigma))} for sample version.
}
\description{
\code{p<u} is OK.
}
\details{
\eqn{\eta} by independent uniform (0,10) variates.
}
