% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{partition}
\alias{partition}
\title{Split the data matrix in train and test sets}
\usage{
partition(y, p = 0.3)
}
\arguments{
\item{y}{input matrix to be split into train and test sets}

\item{p}{fraction of observations to be used for the test set}
}
\description{
Returns a list of two matrices \code{train} and \code{test}.
\code{train} corresponds to the input matrix with a fixed persentage of
entries masked by NA values. \code{test} is the complement of \code{train}
and contains the values of the input matrix in the cells where \code{train}
is NA, while all the other entries are filled by NA values.
}
\keyword{internal}
