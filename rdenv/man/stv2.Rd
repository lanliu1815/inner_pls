% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Starting Values.R
\name{stv2}
\alias{stv2}
\title{Revised version of \code{\link{stv}}.}
\usage{
stv2(M, U, u, Gam.t, PRINT = F, c = NULL, p = NULL)
}
\arguments{
\item{c, p}{Scale added to \eqn{P_U} and the dimension of predictors, which is
only in need when specifying c.}
}
\value{
\code{list(init=init_min, obj=obj_min, init1=init1, init2=init2, init3=init3,
init4=init4, obj1=obj1, obj2=obj2, obj3=obj3, obj4=obj4,
idx1=idx1, idx2=idx2, idx3=idx3, idx4=idx4)}
}
\description{
This function is for studying and testing.
}
\details{
The 4 potential starting values are:
\itemize{
  \item Normal U by eigenvectors of \code{M+U},
  \item Standardized U by eigenvectors of \code{M+U},
  \item Normal U by eigenvectors of \code{M},
  \item Standardized U by eigenvectors of \code{M}.
}
We can choose to print the distance between the estimator and the true
parameter and the values of the objective function of each of the 4 sets
of choices.

Moreover, the 4 candidates, the values of the objective function on them and
 the indices will all be returned.

Finally, two other more candidates are added if \code{c} is not \code{NULL},
 for studying.
}
