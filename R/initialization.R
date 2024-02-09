
#' @title  Initialization of the generalized matrix factorization model
#' @description ...
#' @import svd
#' @keywords internal
init.param = function (
    Y,
    X = NULL,
    Z = NULL,
    ncomp = 2,
    family = gaussian(),
    method = "svd",
    niter = 0,
    values = list(),
    verbose = FALSE
) {

  # Initialize U, V and beta using the selected method
  init = NULL
  if (method == "glm") {
    init = init.param.glm(Y, X, Z, ncomp, family, verbose)
  } else if (method == "svd") {
    init = init.param.svd(Y, X, Z, ncomp, family, niter, verbose)
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
#' @description ...
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
#' @description ...
#' @import svd
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
  s = svd::propack.svd(y - xb - az, neig = d)
  U = s$u %*% diag(sqrt(s$d))
  V = s$v %*% diag(sqrt(s$d))
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
      s = svd::propack.svd(y - xb - az, neig = d)
      U = s$u %*% diag(sqrt(s$d))
      V = s$v %*% diag(sqrt(s$d))
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


#' @title GLM-SVD initialization
#' @description ...
#' @import svd
#' @keywords internal
init.param.glm = function (
    Y,
    X = NULL,
    Z = NULL,
    ncomp = 2,
    family = poisson(),
    verbose = FALSE
) {

  n = nrow(Y)
  m = ncol(Y)
  d = ncomp

  # column-specific covariate vector initialization
  if (verbose) cat(" Initialization: column-specific covariates \n")
  res = c()
  bx = c()
  for (j in 1:ncol(Y)) {
    yj = Y[,j]
    if (family$family == "binomial") {
      isnaj = is.na(yj)
      mj = mean(yj, na.rm = TRUE)
      yj[isnaj] = rbinom(n = sum(isnaj), size = 1, prob = mj)
      yj = as.factor(yj)
    }

    if (!is.null(X)) {
      m = glm(yj ~ X - 1, family = family)
      bx = rbind(bx, m$coefficients)
      res = cbind(res, m$residuals)
    } else {
      res = cbind(res, yj)
    }
  }

  # partial linear predictor
  eta = NULL
  if (!is.null(X)) {
    eta = tcrossprod(X, bx)
  }

  # row-specific covariate vector initialization
  if (verbose) cat(" Initialization: row-specific covariates \n")
  res = c()
  bz = c()
  for (i in 1:nrow(Y)) {
    yi = Y[i,]
    if (family$family == "binomial") {
      isnai = is.na(yi)
      mi = mean(yi, na.rm = TRUE)
      yi[isnai] = rbinom(n = sum(isnai), size = 1, prob = mi)
      yi = as.factor(yi)
    }

    if (!is.null(Z)) {
      m = glm(yi ~ Z - 1, family = family, offset = eta[i,])
      bz = rbind(bz, m$coefficients)
      res = rbind(res, m$residuals)
    } else {
      res = rbind(res, yi)
    }
  }

  # partial linear predictor
  eta = NULL
  if (!is.null(Z)) {
    eta = tcrossprod(bz, Z)
  }

  # residual matrix factorization for initializing the latent components U * Vt
  if (verbose) cat(" Initialization: latent scores and loadings \n")
  s = svd(res, nu = d, nv = d)
  u = s$u %*% diag(sqrt(s$d[1:d]))
  v = s$v %*% diag(sqrt(s$d[1:d]))

  # Compute the initial vector of dispersion parameters
  eta = tcrossprod(u, v)
  if (!is.null(X)) eta = eta + tcrossprod(X, bx)
  if (!is.null(Z)) eta = eta + tcrossprod(bz, Z)
  mu = family$linkinv(eta)
  var = family$variance(mu)
  phi = colMeans((Y - mu)^2 / var, na.rm = TRUE)

  # covariate matrices initialization when there are no covariates
  if (is.null(X)) bx = matrix(0, nrow = ncol(Y), ncol = 0)
  if (is.null(Z)) bz = matrix(0, nrow = nrow(Y), ncol = 0)

  # output
  list(U = u, V = v, A = bz, B = bx, phi = phi)
}


#' @title SVD initialization
#' @description ...
#' @import svd
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
    s = svd::propack.svd(uv, neig = d)
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
    s = svd::propack.svd(uv, neig = d)
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
    s = svd::propack.svd(uv, neig = d)
    U = s$u %*% diag(sqrt(s$d))
    V = s$v %*% diag(sqrt(s$d))
  }

  # If both U and V are unspecified, compute them via incomplete SVD
  # calculated over the regression residuals
  if (is.null(U) & is.null(V)) {
    if (verbose) cat(" Initialization: latent scores and loadings \n")
    s = svd::propack.svd(y - xb - az, neig = d)
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

