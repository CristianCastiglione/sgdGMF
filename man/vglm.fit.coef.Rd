% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vglmfit.R
\name{vglm.fit.coef}
\alias{vglm.fit.coef}
\title{Estimate the coefficients of a vector generalized linear model}
\usage{
vglm.fit.coef(
  Y,
  X,
  family = gaussian(),
  weights = NULL,
  offset = NULL,
  parallel = FALSE,
  nthreads = 1,
  clust = NULL
)
}
\arguments{
\item{Y}{\eqn{n \times m} matrix of response variables}

\item{X}{\eqn{n \times p} matrix of covariates}

\item{family}{a \code{glm} family (see \code{\link{family}} for more details)}

\item{weights}{\eqn{n \times m} matrix of weighting values}

\item{offset}{\eqn{n \times m} matrix of offset values}

\item{parallel}{if \code{TRUE}, allows for parallel computing using the \code{foreach} package}

\item{nthreads}{number of cores to be used in parallel (only if \code{parallel=TRUE})}

\item{clust}{registered cluster to be used for distributing the computations (only if \code{parallel=TRUE})}
}
\description{
Estimate the coefficients of a vector generalized linear model via parallel
iterative re-weighted least squares. Computations can be performed in parallel
to speed up the execution.
}
\keyword{internal}
