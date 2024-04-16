
#' @title Estimate the coefficients of a multivariate linear model
#'
#' @description
#' Estimate the coefficients of a multivariate linear model via ordinary least squares.
#'
#' @param Y \eqn{n \times m} matrix of response variables
#' @param X \eqn{n \times p} matrix of covariates
#' @param offset \eqn{n \times m} matrix of offset values
#'
#' @keywords internal
ols.fit.coef = function (
    Y, X, offset = NULL
) {
  # Set the offset matrix
  if (is.null(offset)) offset = 0

  # Parameter estimation
  XtX = crossprod(X)
  XtY = crossprod(X, Y - offset)
  coefs = t(solve(XtX, XtY))

  # Return the parameter estimates
  return (coefs)
}

#' @title Estimate the coefficients of a vector generalized linear model
#'
#' @description
#' Estimate the coefficients of a vector generalized linear model via parallel
#' iterative re-weighted least squares. Computations can be performed in parallel
#' to speed up the execution.
#'
#' @param Y \eqn{n \times m} matrix of response variables
#' @param X \eqn{n \times p} matrix of covariates
#' @param family a \code{glm} family (see \code{\link{family}} for more details)
#' @param offset \eqn{n \times m} matrix of offset values
#' @param parallel if \code{TRUE}, allows for parallel computing using the \code{foreach} package
#' @param nthreads number of cores to be used in parallel (only if \code{parallel=TRUE})
#' @param clust registered cluster to be used for distributing the computations (only if \code{parallel=TRUE})
#'
#' @keywords internal
vglm.fit.coef = function (
    Y, X, family = gaussian(), offset = NULL,
    parallel = FALSE, nthreads = 1, clust = NULL
) {
  # Set the model dimensions
  n = nrow(Y)
  m = ncol(Y)

  # Set the offset matrix
  if (is.null(offset))
    offset = matrix(0, nrow = n, ncol = m)

  # Register the clusters
  if (parallel) {
    nullclust = is.null(clust)
    if (nullclust) {
      ncores = parallel::detectCores() - 1
      ncores = max(1, min(nthreads, ncores))
      clust = parallel::makeCluster(ncores)
      doParallel::registerDoParallel(clust)
    }
  }

  if (!parallel) {
    # Sequential parameter estimation
    coefs = foreach(j = 1:m, .combine = "rbind") %do% {
      yj = as.vector(Y[,j])
      oj = as.vector(offset[,j])
      fit = stats::glm.fit(x = X, y = yj, family = family, offset = oj)
      t(fit$coefficients)
    }

    ## # As an alternative, we may use the following R code,
    ## # which does not depend on the foreach package
    ## coefs = matrix(NA, nrow = m, ncol = p)
    ## for (j in 1:m) {
    ##   yj = as.vector(Y[,j])
    ##   oj = as.vector(offset[,j])
    ##   fit = stats::glm.fit(x = X, y = yj, family = family, offset = oj)
    ##   coefs[j, ] = as.vector(fit$coefficients)
    ## }
  } else {
    # Parallel parameter estimation
    coefs = foreach(j = 1:m, .combine = "rbind") %dopar% {
      yj = as.vector(Y[,j])
      oj = as.vector(offset[,j])
      fit = stats::glm.fit(x = X, y = yj, family = family, offset = oj)
      t(fit$coefficients)
    }
  }

  # Close the connection to the clusters
  if (parallel) {
    if (nullclust) {
      parallel::stopCluster(clust)
    }
  }

  # Return the parameter estimates
  return (coefs)
}
