% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init.R
\name{init.param.glm}
\alias{init.param.glm}
\title{GLM-SVD initialization}
\usage{
init.param.glm(
  Y,
  X = NULL,
  Z = NULL,
  ncomp = 2,
  family = gaussian(),
  type = c("deviance", "pearson", "working"),
  verbose = FALSE,
  parallel = FALSE,
  nthreads = 1
)
}
\description{
Initialize the parameters of a GMF model fitting a sequence of GLMs followed
by a residual SVD decomposition. See \code{\link{init.param}} for more details.
}
\keyword{internal}
