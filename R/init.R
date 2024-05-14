
#' @title Initialize the parameters of a generalized matrix factorization model
#'
#' @description
#' Provide four initialization methods to set the initial values of the parameters of
#' a generalized matrix factorization (GMF) model identified by a \code{\link{glm}} family
#' and a linear predictor of the form \eqn{g(\mu) = \eta = X B^\top + A Z^\top + U V^\top},
#' with bijective link function \eqn{g(\cdot)}.
#' See \code{\link{sgdgmf}} for more details on the model specification.
#'
#' @param Y matrix of responses (\eqn{n \times m})
#' @param X matrix of row-specific fixed effects (\eqn{n \times p})
#' @param Z matrix of column-specific fixed effects (\eqn{q \times m})
#' @param ncomp rank of the latent matrix factorization
#' @param family a family as in the \code{\link{glm}} interface
#' @param method optimization method to be used for the initial fit
#' @param type type of residuals to be used for initializing \code{U} via incomplete SVD decomposition
#' @param niter number of iterations to refine the initial estimate (only if \code{method="ols"} or \code{"svd"})
#' @param values a list of custom initial values for \code{B}, \code{A}, \code{U} and \code{V}
#' @param verbose print the status of the initialization process
#' @param parallel if \code{TRUE}, allows for parallel computing using the \code{foreach} package
#' @param nthreads number of cores to be used in parallel (only if \code{parallel=TRUE})
#' @param savedata if \code{TRUE}, stores a copy of the input data
#'
#' @return
#' A list containing the initial values of \code{B}, \code{A}, \code{U}, \code{V} and \code{phi}.
#'
#' @details
#' If \code{method = "svd"}, the initialization is performed fitting a sequence of linear
#' regressions followed by a residual SVD decomposition.
#' To account for non-Gaussian distribution of the data, regression and
#' decomposition are applied on the transformed response matrix \eqn{Y_h = (g \circ h)(Y)},
#' where \eqn{h(\cdot)} is a function which prevent \eqn{Y_h} to take infinite values.
#' For instance, in the Binomial case \eqn{h(y) = 2 (1-\epsilon) y + \epsilon},
#' while in the Poisson case \eqn{h(y) = y + \epsilon}, where \eqn{\epsilon} is a small
#' positive constant, typically \code{0.1} or \code{0.01}.
#'
#' If \code{method = "glm"}, the initialization is performed by fitting a sequence of
#' generalized linear models followed by a residual SVD decomposition.
#' In particular, we use independent GLM fit with \eqn{y_j \sim X \beta_j} to set \eqn{\beta_j}.
#' Similarly, we fit the model \eqn{y_i \sim Z \alpha_i + o_i} with offset \eqn{o_i = B x_i}
#' to set \eqn{\alpha_j}. Then, \eqn{U} and \eqn{V} are obtained via SVD on the final
#' working residuals.
#'
#'If  \code{method = "random"}, the initialization is performed using independent Gaussian
#' random values for all the parameters in the model
#'
#' If \code{method = "values"}, the initialization is performed using user-specified
#' values provided as an input.
#'
#' @examples
#' ...
#'
#' @export init.gmf.param
init.gmf.param = function (
    Y,
    X = NULL,
    Z = NULL,
    ncomp = 2,
    family = gaussian(),
    method = c("svd", "ols", "glm", "random", "values"),
    type = c("deviance", "pearson", "working", "link"),
    niter = 0,
    values = list(),
    verbose = FALSE,
    parallel = FALSE,
    nthreads = 1,
    savedata = TRUE
) {
  # Set the initialization method
  method = match.arg(method)

  # Initialize the parameters using the selected method
  init = switch(method,
    "svd" = init.param.ols(Y, X, Z, ncomp, family, type, verbose),
    "ols" = init.param.ols(Y, X, Z, ncomp, family, type, verbose),
    "glm" = init.param.glm(Y, X, Z, ncomp, family, type, verbose, parallel, nthreads),
    "random" = init.param.random(Y, X, Z, ncomp),
    "values" = init.param.custom(Y, X, Z, ncomp, family, values, verbose))

  # Save all the initialization options
  init$method = method
  init$family = family
  init$ncomp = ncomp
  init$type = type
  init$verbose = verbose
  init$parallel = parallel
  init$nthreads = nthreads
  init$savedata = savedata

  if (savedata) {
    n = nrow(Y)
    m = ncol(Y)
    p = ifelse(is.null(X), 1, ncol(X))
    q = ifelse(is.null(Z), 1, ncol(Z))

    init$Y = matrix(NA, nrow = n, ncol = m)
    init$X = matrix(NA, nrow = n, ncol = p)
    init$Z = matrix(NA, nrow = m, ncol = q)

    init$Y[] = Y
    init$X[] = if (!is.null(X)) X else matrix(1, nrow = n, ncol = p)
    init$Z[] = if (!is.null(Z)) Z else matrix(1, nrow = m, ncol = q)
  }

  # Set the initialization class
  class(init) = "initgmf"

  # Return the initial estimates
  return (init)
}

# @title OLS-SVD initialization
#
# @description
# Initialize the parameters of a GMF model fitting a sequence of multivariate
# linear regression followed by a residual SVD decomposition. It allows to
# recursively refine the initial estimate by repeating the process a pre-specified
# number of times. See \code{\link{init.gmf.param}} for more details.

#' @rdname init.gmf.param
#' @keywords internal
init.param.ols = function (
    Y,
    X = NULL,
    Z = NULL,
    ncomp = 2,
    family = gaussian(),
    type = c("deviance", "pearson", "working", "link"),
    verbose = FALSE
) {
  # Set the residual type
  type = match.arg(type)

  # Set the covariate matrices
  if (is.null(X)) X = matrix(1, nrow = nrow(Y), ncol = 1)
  if (is.null(Z)) Z = matrix(1, nrow = ncol(Y), ncol = 1)

  # Set the minimum and maximum of the data
  minY = min(Y)
  maxY = max(Y)

  # Set the model dimensions
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  q = ncol(Z)
  d = ncomp

  # Set the data transformation to use for
  # the initialization of the working data
  family = set.family(family)

  # Compute the transformed data
  if (verbose) cat(" Initialization: working data \n")
  if (anyNA(Y)) {
    isna = is.na(Y)
    Y[] = apply(Y, 2, function (x) {
      x[is.na(x)] = mean(x, na.rm = TRUE)
      return (x)
    })
  }
  gY = family$transform(Y)

  # Initialize the parameters
  B = matrix(NA, nrow = m, ncol = p)
  A = matrix(NA, nrow = n, ncol = q)
  U = matrix(NA, nrow = n, ncol = d)
  V = matrix(NA, nrow = m, ncol = d)

  # Compute the initial column-specific regression parameters
  if (verbose) cat(" Initialization: column-specific covariates \n")
  B[] = ols.fit.coef(gY, X, offset = NULL)
  eta = tcrossprod(X, B)

  # Compute the initial row-specific regression parameter
  if (verbose) cat(" Initialization: row-specific covariates \n")
  A[] = ols.fit.coef(t(gY), Z, offset = t(eta))
  eta = eta + tcrossprod(A, Z)

  # Compute the GLM residuals to be decompose via PCA
  mu = family$linkinv(eta)
  mu[mu > maxY] = maxY
  mu[mu < minY] = minY

  res = switch(type,
    "deviance" = sign(Y - mu) * sqrt(abs(family$dev.resids(Y, mu, 1))),
    "pearson" = (Y - mu) / sqrt(family$variance(mu)),
    "working" = (Y - mu) * family$mu.eta(eta) / abs(family$variance(mu)),
    "link" = (gY - eta))

  if (anyNA(res)) res[is.nan(res) | is.na(res)] = 0

  # Compute the initial latent factors via incomplete SVD
  if (verbose) cat(" Initialization: latent scores \n")
  U[] = RSpectra::svds(res, ncomp)$u

  # Compute the initial loading matrix via OLS
  if (verbose) cat(" Initialization: latent loadings \n")
  V[] = ols.fit.coef(gY, U, offset = eta)
  eta = eta + tcrossprod(U, V)

  # Compute the initial vector of dispersion parameters
  mu = family$linkinv(eta)
  var = family$variance(mu)
  phi = colMeans((Y - mu)^2 / var, na.rm = TRUE)

  # Return the obtained initial values
  list(U = U, V = V, A = A, B = B, phi = phi)
}


# @title GLM-SVD initialization
#
# @description
# Initialize the parameters of a GMF model fitting a sequence of GLMs followed
# by a residual SVD decomposition. See \code{\link{init.gmf.param}} for more details.

#' @rdname init.gmf.param
#' @keywords internal
init.param.glm = function (
    Y,
    X = NULL,
    Z = NULL,
    ncomp = 2,
    family = gaussian(),
    type = c("deviance", "pearson", "working", "link"),
    verbose = FALSE,
    parallel = FALSE,
    nthreads = 1
) {

  # Model dimensions
  n = nrow(Y)
  m = ncol(Y)
  d = ncomp

  # Set the residual type
  type = match.arg(type)

  # Set the data transformation to use for
  # the initialization of the working data
  family = set.family(family)

  # Compute the transformed data
  if (verbose) cat(" Initialization: working data \n")
  if (anyNA(Y)) {
    isna = is.na(Y)
    Y[] = apply(Y, 2, function (x) {
      x[is.na(x)] = mean(x, na.rm = TRUE)
      return (x)
    })
  }

  # Set the covariate matrices
  if (is.null(X)) X = matrix(1, nrow = nrow(Y), ncol = 1)
  if (is.null(Z)) Z = matrix(1, nrow = ncol(Y), ncol = 1)

  # Set the minimum and maximum of the data
  minY = min(Y)
  maxY = max(Y)

  # Initialize the mean and variance matrices
  eta = matrix(NA, nrow = n, ncol = m)
  mu  = matrix(NA, nrow = n, ncol = m)
  var = matrix(NA, nrow = n, ncol = m)
  res = matrix(NA, nrow = n, ncol = m)

  # Register and open the connection to the clusters
  clust = NULL
  if (parallel) {
    ncores = parallel::detectCores() - 1
    ncores = max(1, min(nthreads, ncores))
    clust = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(clust)
  }

  # Column-specific covariate vector initialization
  if (verbose) cat(" Initialization: column-specific covariates \n")
  B = vglm.fit.coef(Y, X, family = family, offset = NULL,
                    parallel = parallel, nthreads = nthreads, clust = clust)

  # Update the linear predictor
  eta[] = tcrossprod(X, B)

  # Row-specific covariate vector initialization
  if (verbose) cat(" Initialization: row-specific covariates \n")
  A = vglm.fit.coef(t(Y), Z, family = family, offset = t(eta),
                    parallel = parallel, nthreads = nthreads, clust = clust)

  # Update the linear predictor and the conditional mean matrix
  eta[] = eta + tcrossprod(A, Z)
  mu[] = family$linkinv(eta)

  mu[mu > maxY] = maxY
  mu[mu < minY] = minY

  # Compute the residuals
  res[] = switch(type,
    "deviance" = sign(Y - mu) * sqrt(abs(family$dev.resids(Y, mu, 1))),
    "pearson" = (Y - mu) / sqrt(abs(family$variance(mu))),
    "working" = (Y - mu) * family$mu.eta(eta) / abs(family$variance(mu)),
    "link" = (family$transform(Y) - eta))

  if (anyNA(res)) res[is.na(res) | is.nan(res)] = 0

  # Initialize the latent factors via residual SVD
  if (verbose) cat(" Initialization: latent scores \n")
  U = RSpectra::svds(res, k = ncomp)$u

  # Initialize the loading matrix via GLM regression
  if (verbose) cat(" Initialization: latent loadings \n")
  V = vglm.fit.coef(Y, U, family = family, offset = eta,
                    parallel = parallel, nthreads = nthreads, clust = clust)

  # Close the connection to the clusters
  if (parallel) parallel::stopCluster(clust)

  # Compute the initial vector of dispersion parameters
  eta[] = eta + tcrossprod(U, V)
  mu[] = family$linkinv(eta)
  var[] = family$variance(mu)
  phi = colMeans((Y - mu)^2 / var, na.rm = TRUE)

  # Output
  list(U = U, V = V, A = A, B = B, phi = phi)
}


#' @rdname init.gmf.param
#' @keywords internal
init.param.random = function (
    Y,
    X = NULL,
    Z = NULL,
    ncomp = 2,
    family = gaussian(),
    sigma = 1
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
  sd = 1e-01 * sigma
  U = array(rnorm(prod(dimU)) / prod(dimU) * sd, dimU)
  V = array(rnorm(prod(dimV)) / prod(dimV) * sd, dimV)
  B = array(rnorm(prod(dimB)) / prod(dimB) * sd, dimB)
  A = array(rnorm(prod(dimA)) / prod(dimA) * sd, dimA)
  phi = rep(1, length = m)

  # output
  list(U = U, V = V, A = A, B = B, phi = phi)
}

# @title SVD initialization
#
# @description
# Initialize the parameters of a GMF model using custom values provided by the user
# and estimating the unspecified parameters using the same procedure described in
# \code{\link{init.param.svd}}. See \code{\link{init.gmf.param}} for more details.

#' @rdname init.gmf.param
#' @keywords internal
init.param.custom = function (
    Y,
    X = NULL,
    Z = NULL,
    ncomp = 2,
    family = gaussian(),
    values = list(),
    verbose = FALSE
) {

  # Set the covariate matrices
  if (is.null(X)) X = matrix(1, nrow = n, ncol = 1)
  if (is.null(Z)) Z = matrix(1, nrow = n, ncol = 1)

  # Data dimensions
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  q = ncol(Z)
  d = ncomp

  # Set the error message
  error.message = function (mat, nr, nc)
    gettextf("Incompatible dimensions: dim(%s) != c(%d, %d).", mat, nr, nc)

  # Check if all the dimensions are compatible
  if (any(dim(values$B) != c(m, p))) stop(error.message("B", m, p), call. = FALSE)
  if (any(dim(values$A) != c(n, q))) stop(error.message("A", n, q), call. = FALSE)
  if (any(dim(values$U) != c(n, d))) stop(error.message("U", n, d), call. = FALSE)
  if (any(dim(values$V) != c(m, d))) stop(error.message("V", m, d), call. = FALSE)

  # Compute the initial vector of dispersion parameters
  eta = tcrossprod(cbind(X, values$A, values$U), cbind(values$B, Z, values$V))
  mu = family$linkinv(eta)
  var = family$variance(mu)
  values$phi = colMeans((Y - mu)^2 / var, na.rm = TRUE)

  # Return the obtained initial values
  return(values)
}

