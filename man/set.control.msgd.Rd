% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/control.R
\name{set.control.msgd}
\alias{set.control.msgd}
\title{Check and set the control parameters for the memoized-SGD algorithm}
\usage{
set.control.msgd(
  normalize = TRUE,
  maxiter = 100,
  eps = 1e-08,
  nafill = 10,
  tol = 1e-05,
  size = 100,
  burn = 90,
  rate0 = 0.01,
  decay = 1,
  damping = 1e-04,
  rate1 = 0.05,
  rate2 = 0.01,
  parallel = FALSE,
  nthreads = 1,
  verbose = FALSE,
  frequency = 10,
  progress = FALSE
)
}
\description{
Check if the input control parameters are allowed and set them to default
values if they are not. Returns a list of well-defined control parameters.
}
