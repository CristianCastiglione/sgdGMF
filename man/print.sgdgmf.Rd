% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sgdGMF-class.R
\name{print.sgdgmf}
\alias{print.sgdgmf}
\title{Print the fundamental characteristics of a GMF}
\usage{
\method{print}{sgdgmf}(x, ...)
}
\arguments{
\item{x}{an object of class \code{sgdgmf}}

\item{...}{further arguments passed to or from other methods}
}
\description{
Print some summary information of a GMF model.
}
\examples{
library(sgdGMF)

# Generate data from a Poisson model
data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())

# Fit a GMF model with 3 latent factors
gmf = sgdgmf.fit(data$Y, ncomp = 3, family = poisson())

# Print the GMF object
print(gmf)

}
