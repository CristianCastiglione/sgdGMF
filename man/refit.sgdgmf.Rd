% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sgdGMF-class.R
\name{refit.sgdgmf}
\alias{refit.sgdgmf}
\title{Refine the final estimate of a GMF model}
\usage{
\method{refit}{sgdgmf}(
  object,
  normalize = TRUE,
  verbose = FALSE,
  parallel = FALSE,
  nthreads = 1
)
}
\arguments{
\item{object}{an object of class \code{sgdgmf}}

\item{normalize}{if \code{TRUE}, normalize \code{U} and \code{V} to uncorrelated Gaussian \code{U} and upper triangular \code{V} with positive diagonal}

\item{parallel}{if \code{TRUE}, use parallel computing using the \code{foreach} package}

\item{nthreads}{number of cores to be used in the \code{"glm"} method}
}
\description{
Refine the estimated latent scores of a GMF model via IRWLS
}
