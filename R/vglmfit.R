
#' @title Estimate the coefficients of a multivariate linear model
#' @description ...
#' @keywords internal
ols.fit.coef = function (
    Y, X, family, offset = NULL
) {
  # Set the offset matrix
  if (is.null(offset)) {
    n = nrow(Y); m = ncol(Y)
    offset = matrix(0, nrow = n, ncol = m)
  }

  # Parameter estimation
  XtX = crossprod(X)
  XtY = crossprod(X, Y - offset)
  coefs = t(solve(XtX, XtY))

  # Return the parameter estimates
  return (coefs)
}

#' @title Estimate the coefficients of a vector generalized linear model
#' @description ...
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
