
#' @title Initialize the parameters of a generalized matrix factorization model
#'
#' @description
#' Provide four initialization methods to set the initial values of the parameters of
#' a generalized matrix factorization (GMF) model identified by a \code{\link{glm}} family
#' and a linear predictor of the form \eqn{g(\mu) = \eta = X B^\top + \Gamma Z^\top + U V^\top},
#' with bijective link function \eqn{g(\cdot)}.
#' See \code{\link{sgdgmf}} for more details on the model specification.
#'
#' @param Y matrix of responses (\eqn{n \times m})
#' @param X matrix of row fixed effects (\eqn{n \times p})
#' @param Z matrix of column fixed effects (\eqn{q \times m})
#' @param ncomp rank of the latent matrix factorization (default 2)
#' @param family a family as in the \code{\link{glm}} interface (default \code{gaussian()})
#' @param method optimization method:  \code{"svd"} (default), \code{"glm"}, \code{"random"}, \code{"values"}
#' @param niter number of iterations to refine the initial estimate (defult 0)
#' @param values a list of custom initial values for \code{B}, \code{A}, \code{U} and \code{V}, or a subset of them
#' @param verbose print the status of the initialization process
#'
#' @return
#' A list containing the initial values of \code{B}, \code{A}, \code{U}, \code{V} and \code{phi}.
#'
#' @details
#' \code{method="svd"}: the initialization is performed fitting a sequence of linear
#' regressions followed by a residual SVD decomposition.
#' To account for the non-Gaussian distribution of the data, regression and
#' decomposition are applied on the transformed response matrix \eqn{Y_h = g(h(Y))},
#' where \eqn{h(\cdot)} is a function which prevent \eqn{Y_h} to take infinite values.
#' For instance, in the Binomial case \eqn{h(y) = 2 (1-\epsilon) y + \epsilon},
#' while in the Poisson case \eqn{h(y) = y + \epsilon}, where \eqn{\epsilon} is a small
#' positive constant, typically \code{0.1} or \code{0.01}.
#'
#' \code{method="glm"}: the initialization is performed by fitting a sequence of
#' generalized linear models followed by a residual SVD decomposition.
#' In particular, we use independent GLM fit \eqn{y_j \sim X \beta_j} to set \eqn{\beta_j}.
#' Similarly, we fit the model \eqn{y_i \sim Z \gamma_i + o_i} with offset \eqn{o_i = B x_i}
#' to set \eqn{\gamma_j}. Then, \eqn{U} and \eqn{V} are obtained via SVD on the final
#' working residuals.
#'
#' \code{method="random"}: the initialization is performed using independent Gaussian
#' random values for all the parameters in the model
#'
#' \code{method="values"}: the initialization is performed using the user-specified
#' values provided as an input. The parameters not provided by the user are set
#' using the same strategy described in \code{method="svd"}.
#'
#'
#' @importFrom svd propack.svd
#'
#' @examples
#' ...
#'
#' @export
init.param = function (
    Y,
    X = NULL,
    Z = NULL,
    ncomp = 2,
    family = gaussian(),
    method = c("svd", "glm", "random", "values"),
    type = c("deviance", "pearson", "working"),
    niter = 0,
    values = list(),
    verbose = FALSE,
    parallel = FALSE,
    nthreads = 1
) {

  method = match.arg(method)

  # Initialize U, V and beta using the selected method
  init = NULL
  if (method == "svd") {
    init = init.param.svd(Y, X, Z, ncomp, family, niter, verbose)
  } else if (method == "glm") {
    if (parallel) {
      init = init.param.glm3(Y, X, Z, ncomp, family, type, verbose, nthreads)
    } else {
      init = init.param.glm2(Y, X, Z, ncomp, family, type, verbose)
    }
  } else if (method == "random") {
    init = init.param.random(Y, X, Z, ncomp)
  } else if (method == "values") {
    init = init.param.custom(Y, X, Z, ncomp, family, values, verbose)
  } else {
    stop("Not allowed initialization method.")
  }

  return (init)
}

#' @title Random initialization
#'
#' @description
#' Initialize the parameters of a GMF model sampling them from an independent
#' Gaussian distribution (see \code{\link{init.param}} for more details)
#'
#' @keywords internal
init.param.random = function (
    Y,
    X = NULL,
    Z = NULL,
    ncomp = 2,
    family = poisson()
) {

  # Derive data dimensions
  n = nrow(Y)
  m = ncol(Y)
  d = ncomp

  # Derive covariate dimensions
  p = 0
  q = 0
  if (!is.null(X)) p = ncol(X) # n x p matrix
  if (!is.null(Z)) q = ncol(Z) # m x q matrix

  # parameter dimensions
  dimU = c(n, d)
  dimV = c(m, d)
  dimB = c(m, p)
  dimA = c(n, q)

  # parameter generation
  sd = 1e-01
  U = array(rnorm(prod(dimU)) / prod(dimU) * sd, dimU)
  V = array(rnorm(prod(dimV)) / prod(dimV) * sd, dimV)
  B = array(rnorm(prod(dimB)) / prod(dimB) * sd, dimB)
  A = array(rnorm(prod(dimA)) / prod(dimA) * sd, dimA)
  phi = rep(1, length = m)

  # output
  list(U = U, V = V, A = A, B = B, phi = phi)
}

#' @title OLS-SVD initialization
#'
#' @description
#' Initialize the parameters of a GMF model fitting a sequence of multivariate
#' linear regression followed by a residual SVD decomposition. It allows to
#' recursively refine the initial estimate by repeating the process a pre-specified
#' number of times. See \code{\link{init.param}} for more details.
#'
#' @keywords internal
init.param.svd = function (
    Y,
    X = NULL,
    Z = NULL,
    ncomp = 2,
    family = poisson(),
    niter = 0,
    verbose = FALSE
) {

  n = nrow(Y)
  m = ncol(Y)
  d = ncomp

  # Select the data transformation to use for
  # the initialization of the working data
  f = set.jitter(family)

  # Compute the transformed data
  if (verbose) cat(" Initialization: working data \n")
  isna = is.na(Y)
  y = matrix(NA, nrow = n, ncol = m)
  # y[!isna] = f(Y[!isna])
  # y[isna] = mean(y[!isna])

  y[] = apply(Y, 2, function (x) {
    x[is.na(x)] = mean(x, na.rm = TRUE)
    return (f(x))
  })

  # Initialize the parameters and sufficient statistics
  # to NULL and zero, respectively
  A = B = U = V = NULL
  xtx = ztz = NULL
  xty = zty = NULL
  xb = az = uv = 0

  # Compute the initial column-specific regression parameters (if any)
  if (!is.null(X)) {
    if (verbose) cat(" Initialization: column-specific covariates \n")
    m = ncol(Y)
    p = ncol(X)
    B = matrix(NA, nrow = m, ncol = p)
    xtx = crossprod(X)
    xty = crossprod(X, y)
    B[] = t(solve(xtx, xty))
    xb = tcrossprod(X, B)
  }

  # Compute the initial row-specific regression parameter (if any)
  if (!is.null(Z)) {
    if (verbose) cat(" Initialization: row-specific covariates \n")
    n = nrow(Y)
    q = ncol(Z)
    A = matrix(NA, nrow = n, ncol = q)
    ztz = crossprod(Z)
    zty = crossprod(Z, t(y - xb))
    A[] = t(solve(ztz, zty))
    az = tcrossprod(A, Z)
  }

  # Compute the initial latent factors via incomplete SVD
  if (verbose) cat(" Initialization: latent scores and loadings \n")
  # sv = svd::propack.svd(y - xb - az, neig = ncomp)
  s = RSpectra::svds(y - xb - az, k = d, nu = d, nv = d)
  if (d == 1) {
    U = s$u %*% sqrt(s$d)
    V = s$v %*% sqrt(s$d)
  } else {
    U = s$u %*% diag(sqrt(s$d))
    V = s$v %*% diag(sqrt(s$d))
  }
  uv = tcrossprod(U, V)

  # Refinement loop (it might be useful if there are many missing values)
  if (niter > 0) {
    if (verbose) cat(" Refinement: |")
    for (iter in 1:niter) {

      if (verbose) cat("=")

      # Refine the initial matrix completion
      y[isna] = xb[isna] + az[isna] + uv[isna]

      # Refine the initial column-specific regression parameters (if any)
      if (!is.null(X)) {
        xty = crossprod(X, y - az - uv)
        B[] = t(solve(xtx, xty))
        xb = tcrossprod(X, B)
      }

      # Refine the initial row-specific regression parameter (if any)
      if (!is.null(Z)) {
        zty = crossprod(Z, t(y - xb - uv))
        A[] = t(solve(ztz, zty))
        az = tcrossprod(A, Z)
      }

      # Refine the initial latent factors via incomplete SV
      # s = svd::propack.svd(y - xb - az, neig = d)
      s = RSpectra::svds(y - xb - az, k = d, nu = d, nv = d)
      if (d == 1) {
        U = s$u * sqrt(s$d)
        V = s$v * sqrt(s$d)
      } else {
        U = s$u %*% diag(sqrt(s$d))
        V = s$v %*% diag(sqrt(s$d))
      }
      uv = tcrossprod(U, V)
    }
    if (verbose) cat("| \n")
  }

  # Compute the initial vector of dispersion parameters
  mu = family$linkinv(xb + az + uv)
  var = family$variance(mu)
  phi = colMeans((Y - mu)^2 / var, na.rm = TRUE)

  # Covariate matrices initialization when there are no covariates
  if (is.null(X)) B = matrix(0, nrow = m, ncol = 0)
  if (is.null(Z)) A = matrix(0, nrow = n, ncol = 0)

  # Return the obtained initial values
  list(U = U, V = V, A = A, B = B, phi = phi)
}


#' @title GLM-SVD initialization (DEPRECATED)
#'
#' @description
#' Initialize the parameters of a GMF model fitting a sequence of GLMs followed
#' by a residual SVD decomposition. See \code{\link{init.param}} for more details.
#'
#' @keywords internal
init.param.glm = function (
    Y,
    X = NULL,
    Z = NULL,
    ncomp = 2,
    family = poisson(),
    type = c("deviance", "pearson", "working"),
    verbose = FALSE
) {
  # Model dimensions
  n = nrow(Y)
  m = ncol(Y)
  d = ncomp

  # Set the residual type
  type = match.arg(type)

  # Compute the transformed data
  if (verbose) cat(" Initialization: working data \n")
  isna = is.na(Y)
  Y[] = apply(Y, 2, function (x) {
    x[is.na(x)] = mean(x, na.rm = TRUE)
    return (x)
  })

  # Initialize the mean and variance matrices
  eta = matrix(NA, nrow = n, ncol = m)
  mu = matrix(NA, nrow = n, ncol = m)
  var = matrix(NA, nrow = n, ncol = m)
  res = matrix(NA, nrow = n, ncol = m)

  # column-specific covariate vector initialization
  if (verbose) cat(" Initialization: column-specific covariates \n")
  B = c()
  for (j in 1:m) {
    yj = as.vector(Y[,j])
    fit = stats::glm.fit(x = X, y = yj, family = family)
    B = rbind(B, fit$coefficients)
  }

  # Update the linear predictor
  eta[] = tcrossprod(X, B)

  # row-specific covariate vector initialization
  if (verbose) cat(" Initialization: row-specific covariates \n")
  A = c()
  for (i in 1:n) {
    yi = as.vector(Y[i,])
    oi = as.vector(eta[i,])
    fit = stats::glm.fit(x = Z, y = yi, family = family, offset = oi)
    A = rbind(A, fit$coefficients)
  }

  # Update the linear predictor
  eta[] = eta + tcrossprod(A, Z)
  mu[] = family$linkinv(eta)

  # Compute the residuals
  if (type == "deviance") res[] = sign(Y - mu) * sqrt(family$dev.resids(Y, mu, 1))
  if (type == "pearson") res[] = (Y - mu) / sqrt(family$variance(mu))
  if (type == "working") res[] = (Y - mu) / family$mu.eta(eta)

  # Initialize the latent factors via residual SVD
  if (verbose) cat(" Initialization: latent scores \n")
  U = RSpectra::svds(res, k = ncomp)$u

  # Initialize the loading matrix via GLM regression
  if (verbose) cat(" Initialization: latent loadings \n")
  V = c()
  for (j in 1:m) {
    yj = as.vector(Y[,j])
    oj = as.vector(eta[,j])
    fit = stats::glm.fit(x = U, y = yj, family = family, offset = oj)
    V = rbind(V, fit$coefficients)
  }

  # Compute the initial vector of dispersion parameters
  eta[] = eta + tcrossprod(U, V)
  mu[] = family$linkinv(eta)
  var[] = family$variance(mu)
  phi = colMeans((Y - mu)^2 / var, na.rm = TRUE)

  # covariate matrices initialization when there are no covariates
  if (is.null(X)) B = matrix(0, nrow = ncol(Y), ncol = 0)
  if (is.null(Z)) A = matrix(0, nrow = nrow(Y), ncol = 0)

  # output
  list(U = U, V = V, A = A, B = B, phi = phi)
}


#' @title GLM-SVD initialization
#'
#' @description
#' Initialize the parameters of a GMF model fitting a sequence of GLMs followed
#' by a residual SVD decomposition. See \code{\link{init.param}} for more details.
#'
#' @keywords internal
init.param.glm2 = function (
    Y,
    X = NULL,
    Z = NULL,
    ncomp = 2,
    family = poisson(),
    type = c("deviance", "pearson", "working"),
    verbose = FALSE
) {
  # Require the needed packages
  # suppressPackageStartupMessages({
  #   require(parallel)
  #   require(doParallel)
  #   require(foreach)
  # })

  # Model dimensions
  n = nrow(Y)
  m = ncol(Y)
  d = ncomp

  # Set the residual type
  type = match.arg(type)

  # Compute the transformed data
  if (verbose) cat(" Initialization: working data \n")
  Y[] = apply(Y, 2, function (x) {
    x[is.na(x)] = mean(x, na.rm = TRUE)
    return (x)
  })

  # Set the covariate matrices
  if (is.null(X)) X = matrix(1, nrow = n, ncol = 1)
  if (is.null(Z)) Z = matrix(1, nrow = n, ncol = 1)

  # Initialize the mean and variance matrices
  eta = matrix(NA, nrow = n, ncol = m)
  mu  = matrix(NA, nrow = n, ncol = m)
  var = matrix(NA, nrow = n, ncol = m)
  res = matrix(NA, nrow = n, ncol = m)

  # Column-specific covariate vector initialization
  if (verbose) cat(" Initialization: column-specific covariates \n")
  B = foreach(j = 1:m, .combine = "rbind") %do% {
    yj = as.vector(Y[,j])
    fit = stats::glm.fit(x = X, y = yj, family = family)
    t(fit$coefficients)
  }

  # Update the linear predictor
  eta[] = tcrossprod(X, B)

  # Row-specific covariate vector initialization
  if (verbose) cat(" Initialization: row-specific covariates \n")
  A = foreach(i = 1:n, .combine = "rbind") %do% {
    yi = as.vector(Y[i,])
    oi = as.vector(eta[i,])
    fit = stats::glm.fit(x = Z, y = yi, family = family, offset = oi)
    fit$coefficients
  }

  # Update the linear predictor and the conditional mean matrix
  eta[] = eta + tcrossprod(A, Z)
  mu[] = family$linkinv(eta)

  # Compute the residuals
  if (type == "deviance") res[] = sign(Y - mu) * sqrt(family$dev.resids(Y, mu, 1))
  if (type == "pearson") res[] = (Y - mu) / sqrt(family$variance(mu))
  if (type == "working") res[] = (Y - mu) / family$mu.eta(eta)

  # Initialize the latent factors via residual SVD
  if (verbose) cat(" Initialization: latent scores \n")
  U = RSpectra::svds(res, k = ncomp)$u

  # Initialize the loading matrix via GLM regression
  if (verbose) cat(" Initialization: latent loadings \n")
  V = foreach(j = 1:m, .combine = "rbind") %do% {
    yj = as.vector(Y[,j])
    oj = as.vector(eta[,j])
    fit = stats::glm.fit(x = U, y = yj, family = family, offset = oj)
    t(fit$coefficients)
  }

  # Compute the initial vector of dispersion parameters
  eta[] = eta + tcrossprod(U, V)
  mu[] = family$linkinv(eta)
  var[] = family$variance(mu)
  phi = colMeans((Y - mu)^2 / var, na.rm = TRUE)

  # Initialization of the regression effects when there are no covariates
  if (is.null(X)) B = matrix(0, nrow = ncol(Y), ncol = 0)
  if (is.null(Z)) A = matrix(0, nrow = nrow(Y), ncol = 0)

  # Output
  list(U = U, V = V, A = A, B = B, phi = phi)
}


#' @title GLM-SVD initialization
#'
#' @description
#' Initialize the parameters of a GMF model fitting a sequence of GLMs followed
#' by a residual SVD decomposition. See \code{\link{init.param}} for more details.
#'
#' @keywords internal
init.param.glm3 = function (
    Y,
    X = NULL,
    Z = NULL,
    ncomp = 2,
    family = poisson(),
    type = c("deviance", "pearson", "working"),
    verbose = FALSE,
    nthreads = 1
) {
  # Require the needed packages
  # suppressPackageStartupMessages({
  #   require(parallel)
  #   require(doParallel)
  #   require(foreach)
  # })

  # Model dimensions
  n = nrow(Y)
  m = ncol(Y)
  d = ncomp

  # Set the residual type
  type = match.arg(type)

  # Compute the transformed data
  if (verbose) cat(" Initialization: working data \n")
  Y[] = apply(Y, 2, function (x) {
    x[is.na(x)] = mean(x, na.rm = TRUE)
    return (x)
  })

  # Set the covariate matrices
  if (is.null(X)) X = matrix(1, nrow = n, ncol = 1)
  if (is.null(Z)) Z = matrix(1, nrow = n, ncol = 1)

  # Initialize the mean and variance matrices
  eta = matrix(NA, nrow = n, ncol = m)
  mu  = matrix(NA, nrow = n, ncol = m)
  var = matrix(NA, nrow = n, ncol = m)
  res = matrix(NA, nrow = n, ncol = m)

  # Register the clusters
  ncores = parallel::detectCores() - 1
  ncores = max(1, min(nthreads, ncores))
  clust = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(clust)

  # Column-specific covariate vector initialization
  if (verbose) cat(" Initialization: column-specific covariates \n")
  B = foreach(j = 1:m, .combine = "rbind") %dopar% {
    yj = as.vector(Y[,j])
    fit = stats::glm.fit(x = X, y = yj, family = family)
    t(fit$coefficients)
  }

  # Update the linear predictor
  eta[] = tcrossprod(X, B)

  # Row-specific covariate vector initialization
  if (verbose) cat(" Initialization: row-specific covariates \n")
  A = foreach(i = 1:n, .combine = "rbind") %dopar% {
    yi = as.vector(Y[i,])
    oi = as.vector(eta[i,])
    fit = stats::glm.fit(x = Z, y = yi, family = family, offset = oi)
    fit$coefficients
  }

  # Update the linear predictor and the conditional mean matrix
  eta[] = eta + tcrossprod(A, Z)
  mu[] = family$linkinv(eta)

  # Compute the residuals
  if (type == "deviance") res[] = sign(Y - mu) * sqrt(family$dev.resids(Y, mu, 1))
  if (type == "pearson") res[] = (Y - mu) / sqrt(family$variance(mu))
  if (type == "working") res[] = (Y - mu) / family$mu.eta(eta)

  # Initialize the latent factors via residual SVD
  if (verbose) cat(" Initialization: latent scores \n")
  U = RSpectra::svds(res, k = ncomp)$u

  # Initialize the loading matrix via GLM regression
  if (verbose) cat(" Initialization: latent loadings \n")
  V = foreach(j = 1:m, .combine = "rbind") %dopar% {
    yj = as.vector(Y[,j])
    oj = as.vector(eta[,j])
    fit = stats::glm.fit(x = U, y = yj, family = family, offset = oj)
    t(fit$coefficients)
  }

  # Close the connection to the clusters
  parallel::stopCluster(clust)

  # Compute the initial vector of dispersion parameters
  eta[] = eta + tcrossprod(U, V)
  mu[] = family$linkinv(eta)
  var[] = family$variance(mu)
  phi = colMeans((Y - mu)^2 / var, na.rm = TRUE)

  # Initialization of the regression effects when there are no covariates
  if (is.null(X)) B = matrix(0, nrow = ncol(Y), ncol = 0)
  if (is.null(Z)) A = matrix(0, nrow = nrow(Y), ncol = 0)

  # Output
  list(U = U, V = V, A = A, B = B, phi = phi)
}


#' @title SVD initialization
#'
#' @description
#' Initialize the parameters of a GMF model using custom values provided by the user
#' and estimating the unspecified parameters using the same procedure described in
#' \code{\link{init.param.svd}}. See \code{\link{init.param}} for more details.
#'
#' @keywords internal
init.param.custom = function (
    Y,
    X = NULL,
    Z = NULL,
    ncomp = 2,
    family = poisson(),
    values = list(),
    verbose = FALSE
) {

  # Data dimensions
  n = nrow(Y)
  m = ncol(Y)
  d = ncomp

  # Select the data transformation to use for
  # the initialization of the working data
  f = set.jitter(family)

  # Compute the transformed data and fill the NA values with the column means
  if (verbose) cat(" Initialization: working data \n")
  isna = is.na(Y)
  y = matrix(NA, nrow = n, ncol = m)
  y[] = apply(Y, 2, function (x) {
    x[is.na(x)] = mean(x, na.rm = TRUE)
    return (f(x))
  })

  # Initialize the parameters and sufficient statistics
  A = B = U = V = NULL
  xtx = ztz = NULL
  xty = zty = NULL
  xb = az = uv = 0

  # Safety checks and parameter assignement
  if (is.list(values)) {
    if (check.dim(values$U, n, d)) U = values$U
    if (check.dim(values$V, m, d)) V = values$V
    if (check.dim(values$A, n, ncol(Z))) A = values$A
    if (check.dim(values$B, m, ncol(X))) B = values$B
    if (!is.null(U) & !is.null(V)) uv = tcrossprod(U, V)
    if (!is.null(A) & !is.null(Z)) az = tcrossprod(A, Z)
    if (!is.null(X) & !is.null(B)) xb = tcrossprod(X, B)

  }

  # Compute the initial column-specific regression parameters (if any)
  if (!is.null(X)) {
    if (is.null(B)) {
      if (verbose) cat(" Initialization: column-specific covariates \n")
      m = ncol(Y)
      p = ncol(X)
      B = matrix(NA, nrow = m, ncol = p)
      xtx = crossprod(X)
      xty = crossprod(X, y - az - uv)
      B[] = t(solve(xtx, xty))
    }
    xb = tcrossprod(X, B)
  }

  # Compute the initial row-specific regression parameter (if any)
  if (!is.null(Z)) {
    if (is.null(A)) {
      if (verbose) cat(" Initialization: row-specific covariates \n")
      n = nrow(Y)
      q = ncol(Z)
      A = matrix(NA, nrow = m, ncol = q)
      ztz = crossprod(Z)
      zty = crossprod(Z, t(y - xb - uv))
      A[] = t(solve(ztz, zty))
    }
    az = tcrossprod(A, Z)
  }

  # If both U and V are provided, orthogonalize them via incomplete SVD
  if (!is.null(U) & !is.null(V)) {
    uv = tcrossprod(U, V)
    s = RSpectra::svds(uv, k = d)
    U = s$u %*% diag(sqrt(s$d))
    V = s$v %*% diag(sqrt(s$d))
  }

  # If U is unspecified, compute it via penalized least squares
  if (!is.null(U) & is.null(V)) {
    if (verbose) cat(" Initialization: latent scores \n")
    # Compute the unnormalized U via penalized least squares
    vtv = crossprod(V) + diag(d)
    vty = crossprod(V, t(y - xb - za))
    U = t(solve(vtv, vty))
    # Orthogonalize U and V via incomplete SVD
    uv = tcrossprod(U, V)
    s = RSpectra::svds(uv, k = d)
    U = s$u %*% diag(sqrt(s$d))
    V = s$v %*% diag(sqrt(s$d))
  }

  # If V is unspecified, compute it via least squares
  if (is.null(U) & !is.null(V)) {
    if (verbose) cat(" Initialization: latent loadings \n")
    # Compute the unnormalized V via penalized least squares
    utu = crossprod(U)
    uty = crossprod(U, t(y - xb - za))
    V = t(solve(utu, uty))
    # Orthogonalize U and V via incomplete SVD
    uv = tcrossprod(U, V)
    s = RSpectra::svds(uv, k = d)
    U = s$u %*% diag(sqrt(s$d))
    V = s$v %*% diag(sqrt(s$d))
  }

  # If both U and V are unspecified, compute them via incomplete SVD
  # calculated over the regression residuals
  if (is.null(U) & is.null(V)) {
    if (verbose) cat(" Initialization: latent scores and loadings \n")
    s = RSpectra::svds(y - xb - az, k = d)
    U = s$u %*% diag(sqrt(s$d))
    V = s$v %*% diag(sqrt(s$d))
    uv = tcrossprod(U, V)
  }

  # Compute the initial vector of dispersion parameters
  mu = family$linkinv(xb + az + uv)
  var = family$variance(mu)
  phi = colMeans((Y - mu)^2 / var, na.rm = TRUE)

  # Regression matrices initialization when there are no covariates
  if (is.null(X)) B = matrix(0, nrow = m, ncol = 0)
  if (is.null(Z)) A = matrix(0, nrow = n, ncol = 0)

  # Return the obtained initial values
  list(U = U, V = V, A = A, B = B, phi = phi)
}

