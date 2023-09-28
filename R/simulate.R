
#' @title Simulate a matrix of non-Gaussian observations
#'
#' @description
#' ...
#'
#' @param n ...
#' @param m ...
#' @param ncomp ...
#' @param ngroup ...
#' @param nx ...
#' @param nz ...
#' @param intercept_row ...
#' @param intercept_col ...
#' @param weight ...
#' @param offset ...
#' @param family ...
#'
#' @return
#' ...
#'
#' @details
#' ...
#'
#' @references
#' ...
#'
#' @import
#' svd
#'
#' @examples
#' ...
#'
#' @export
gmf.simulate = function (
    n = 100,
    m = 20,
    ncomp = 2,
    ngroup = 3,
    nx = NULL,
    nz = NULL,
    intercept_row = FALSE,
    intercept_col = FALSE,
    weight = 0.3,
    offset = NULL,
    family = poisson()
    ){

  # Define the groups
  idx = 1:n
  size = floor(n / ngroup) + 1
  chunks = split(idx, ceiling(seq_along(idx) / size))
  ngroup = length(chunks)

  # Create a selection matrix allocating the i-th obs to the j-th group
  S = matrix(0, nrow = n, ncol = ngroup)
  for (j in 1:ngroup) {
    S[chunks[[j]],j] = 1
  }

  # Create the group vector
  groups = drop(rowSums(S %*% diag(1:ngroup)))

  # Create a correlation matrix and extract its Cholesky factor
  C = matrix(rnorm(ngroup * ngroup), nrow = ngroup, ncol = ngroup)
  C = crossprod(C)
  C = 0.8 * C / norm(C, "F")
  C = diag(diag(C))
  L = chol(C)

  # Simulate the latent variables
  V = matrix(rnorm(m * ncomp), nrow = m, ncol = ncomp) # latent loadings
  U = matrix(rnorm(n * ncomp), nrow = n, ncol = ncomp) # latent scores

  # Rotate the latent variables
  for (j in 1:ncomp) {
    ej = rnorm(ngroup, mean = 0, sd = 1)
    mj = drop(S %*% L %*% ej)
    U[,j] = (1 - weight) * mj + weight * U[,j]
  }

  # Simulate columns-specific covariates
  if (is.null(nx)) {
    X = matrix(0, nrow = n, ncol = 0)
    B = matrix(0, nrow = m, ncol = 0)
  } else {
    X = matrix(rnorm(n * nx), nrow = n, ncol = nx)
    B = matrix(runif(m * nx, -1, 1), nrow = m, ncol = nx)
  }

  if (intercept_col) {
    X = cbind(1, X)
    B = cbind(runif(m, -1, +1), B)
  }

  # Simulate row-specific covariates
  if (is.null(nz)) {
    Z = matrix(0, nrow = m, ncol = 0)
    A = matrix(0, nrow = n, ncol = 0)
  } else {
    Z = matrix(rnorm(m * nz), m, nz)
    A = matrix(runif(n * nz, -1, 1), n, nz)
  }

  if (intercept_row) {
    Z = cbind(1, Z)
    A = cbind(runif(n, -1, +1), A)
  }

  # Set the offset component
  offset = set.offset(offset, n, m)

  # Compute the linear predictor
  xb = tcrossprod(X, B)
  az = tcrossprod(A, Z)
  uv = tcrossprod(U, V)
  eta = xb + az + uv + offset

  # Define the random generator
  if (family$family == "poisson"){
    rr = function(n, l) rpois(n, lambda = l)
  }
  if (family$family == "binomial") {
    rr = function(n, p) rbinom(n, size = 1, prob = p)
  }
  if (family$family == "Gamma") {
    rr = function(n, p) rgamma(n, shape = 1, scale = p)
  }
  if (family$family == "gaussian"){
    rr = function(n, m) rnorm(n, mean = m, sd = 0.5)
  }
  if (substr(family$family, 1, 17) == "Negative Binomial"){
    theta = as.numeric(substr(family$family, 19, nchar(family$family) - 1))
    rr = function(n, p) { rnegbin(n, p, theta) }
  }

  # Simulate the response matrix
  Y = apply(eta, 1:2, function(x) rr(1, family$linkinv(x)))

  # Return the simulated data
  list(y = Y, # response data matrix
       u = U, # latent scores
       v = V, # latent loadings
       x = X, # column-specific covariates
       z = Z, # row-specific covariates
       beta.x = B, # column-specific coefficient matrix
       beta.z = A, # row-specific coefficient matrix
       offset = offset, # offset
       eta = eta, # linear predictor
       n = n, # n rows
       m = m, # n columns
       p = nx, # n column-covariate
       q = nz, # n row-covariates
       ncomp = ncomp, # n latent components
       ngroup = ngroup, # number of sub-groups
       groups = groups, # group indicator
       family = family # family
  )
}

set.offset = function (offset, n, m) {

  if (is.null(offset)) {
    # Check whether offset is NULL
    offset = 0
  } else if (is.matrix(offset)) {
    # Check whether offset is a matrix
    if (all(dim(offset) == c(m, n))) {
      offset = t(offset)
    }
    if (any(dim(offset) != c(n, m))) {
      warning("If offset is a matrix, it must have dimensions n x m.")
      offset = 0
    }
  } else if (is.vector(offset)) {
    # Check whether offset is a vector (or a scalar)
    if (!(length(offset) %in% c(1, n, m))) {
      warning("If offset is a vector, it must have dimensions n x 1 or m x 1.")
      offset = 0
    }
    if (length(offset) == 1) {
      offset = drop(offset)
    }
    if (length(offset) == n) {
      offset = tcrossprod(offset, rep(1, length = m))
    }
    if (length(offset) == m) {
      offset = tcrossprod(rep(1, length = n), offset)
    }
  } else {
    # Check whether offset is not a vector, a matrix, or a NULL value
    warning("The parameter offset must be a vector, a matrix or a NULL value.")
    offset = 0
  }

  # Return the fixed offset
  return (offset)
}
