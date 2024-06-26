% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{normalize.uv}
\alias{normalize.uv}
\title{Normalize the matrices U and V}
\usage{
normalize.uv(U, V, method = c("qr", "svd"))
}
\description{
Rotate U and V using either QR or SVD decompositions.
The QR methods rotate U and V in such a way to obtain an orthogonal U
and a lower triangular V.  The SVD method rotate U and V in such a way
to obtain an orthogonal U and a scaled orthogonal V.
}
\keyword{internal}
