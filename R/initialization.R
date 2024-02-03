
#' @title  Initialization of the generalized matrix factorization model
#' @description ...
#' @import svd
#' @keywords internal
gmf.init = function (
    Y,
    X = NULL,
    Z = NULL,
    d = min(dim(Y)),
    family = gaussian(),
    method = "svd",
    niter = 0,
    values = list(),
    verbose = FALSE
) {

  # Initialize U, V and beta using the selected method
  init = NULL
  if (method == "glm-svd") {
    init = gmf.init.glm(Y, X, Z, d, family, verbose)
  } else if (method == "ls-svd") {
    init = gmf.init.svd(Y, X, Z, d, family, niter, verbose)
  } else if (method == "svd") {
    init = gmf.init.svd(Y, X, Z, d, family, niter, verbose)
  } else if (method == "random") {
    init = gmf.init.random(Y, X, Z, d)
  } else if (method == "values") {
    init = gmf.init.custom(Y, X, Z, d, family, values, verbose)
  } else {
    stop("Not allowed initialization method.")
  }

  return (init)
}

#' @title Random initialization
#' @description ...
#' @keywords internal
gmf.init.random = function (
    Y,
    X = NULL,
    Z = NULL,
    d = min(dim(Y)),
    family = poisson()
) {

  # Derive data dimensions
  n = nrow(Y)
  m = ncol(Y)

  # Derive covariate dimensions
  p = 0
  q = 0
  if (!is.null(X)) p = ncol(X) # n x p matrix
  if (!is.null(Z)) q = ncol(Z) # m x q matrix

  # parameter dimensions
  dim_u = c(n, d)
  dim_v = c(m, d)
  dim_bx = c(m, p)
  dim_bz = c(n, q)

  # parameter generation
  sd = 1e-01
  u = array(rnorm(prod(dim_u)) / prod(dim_u) * sd, dim_u)
  v = array(rnorm(prod(dim_v)) / prod(dim_v) * sd, dim_v)
  bx = array(rnorm(prod(dim_bx)) / prod(dim_bx) * sd, dim_bx)
  bz = array(rnorm(prod(dim_bz)) / prod(dim_bz) * sd, dim_bz)
  phi = rep(1, length = m)

  # output
  list(u = u, v = v, bz = bz, bx = bx, phi = phi)
}

#' @title OLS-SVD initialization
#' @description ...
#' @import svd
#' @keywords internal
gmf.init.svd = function (
    Y,
    X = NULL,
    Z = NULL,
    d = min(dim(Y)),
    family = poisson(),
    niter = 0,
    verbose = FALSE
) {

  # Select the data transformation to use for
  # the initialization of the working data
  f = set.jitter(family)

  # Compute the transformed data
  if (verbose) cat(" Initialization: working data \n")
  isna = is.na(Y)
  y = matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
  # y[!isna] = f(Y[!isna])
  # y[isna] = mean(y[!isna])

  y[] = apply(Y, 2, function (x) {
    x[is.na(x)] = mean(x, na.rm = TRUE)
    return (f(x))
  })

  # Initialize the parameters and sufficient statistics
  # to NULL and zero, respectively
  a = b = u = v = NULL
  xtx = ztz = NULL
  xty = zty = NULL
  xb = az = uv = 0

  # Compute the initial column-specific regression parameters (if any)
  if (!is.null(X)) {
    if (verbose) cat(" Initialization: column-specific covariates \n")
    xtx = crossprod(X)
    xty = crossprod(X, y)
    b = t(solve(xtx, xty))
    xb = tcrossprod(X, b)
  }

  # Compute the initial row-specific regression parameter (if any)
  if (!is.null(Z)) {
    if (verbose) cat(" Initialization: row-specific covariates \n")
    ztz = crossprod(Z)
    zty = crossprod(Z, t(y - xb))
    a = t(solve(ztz, zty))
    az = tcrossprod(a, Z)
  }

  # Compute the initial latent factors via incomplete SVD
  if (verbose) cat(" Initialization: latent scores and loadings \n")
  s = svd::propack.svd(y - xb - az, neig = d)
  u = s$u %*% diag(sqrt(s$d))
  v = s$v %*% diag(sqrt(s$d))
  uv = tcrossprod(u, v)

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
        b = t(solve(xtx, xty))
        xb = tcrossprod(X, b)
      }

      # Refine the initial row-specific regression parameter (if any)
      if (!is.null(Z)) {
        zty = crossprod(Z, t(y - xb - uv))
        a = t(solve(ztz, zty))
        az = tcrossprod(a, Z)
      }

      # Refine the initial latent factors via incomplete SV
      s = svd::propack.svd(y - xb - az, neig = d)
      u = s$u %*% diag(sqrt(s$d))
      v = s$v %*% diag(sqrt(s$d))
      uv = tcrossprod(u, v)
    }
    if (verbose) cat("| \n")
  }

  # Compute the initial vector of dispersion parameters
  mu = family$linkinv(xb + az + uv)
  var = family$variance(mu)
  phi = colMeans((Y - mu)^2 / var, na.rm = TRUE)

  # Covariate matrices initialization when there are no covariates
  if (is.null(X)) bv = matrix(0, nrow = ncol(Y), ncol = 0)
  if (is.null(Z)) bu = matrix(0, nrow = nrow(Y), ncol = 0)

  # Return the obtained initial values
  list(u = u, v = v, bz = a, bx = b, phi = phi)
}


#' @title GLM-SVD initialization
#' @description ...
#' @import svd
#' @keywords internal
gmf.init.glm = function (
    Y,
    X = NULL,
    Z = NULL,
    d = min(dim(Y)),
    family = poisson(),
    verbose = FALSE
) {

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
  list(u = u, v = v, bz = bz, bx = bx, phi = phi)
}


#' @title SVD initialization
#' @description ...
#' @import svd
#' @keywords internal
gmf.init.custom = function (
    Y,
    X = NULL,
    Z = NULL,
    d = min(dim(Y)),
    family = poisson(),
    values = list(),
    verbose = FALSE
) {

  # Safety checks
  check = function (object, n, m) {
    flag = FALSE
    # check if object exists
    if (!is.null(object)) {
      # check if object is a numeric matrix
      if (is.numeric(object) & is.matrix(object)) {
        # check is object is a matrix of appropriate dimensions
        if (nrow(object) == n & ncol(object) == m) {
          flag = TRUE
        }
      }
    }
    return (flag)
  }

  # Select the data transformation to use for
  # the initialization of the working data
  f = set.jitter(family)

  # Compute the transformed data and fill the NA values with the column means
  if (verbose) cat(" Initialization: working data \n")
  isna = is.na(Y)
  y = matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
  # y[!isna] = f(Y[!isna])
  # y[isna] = mean(y[!isna])

  y[] = apply(Y, 2, function (x) {
    x[is.na(x)] = mean(x, na.rm = TRUE)
    return (f(x))
  })

  # Initialize the parameters and sufficient statistics
  a = b = u = v = NULL
  xtx = ztz = NULL
  xty = zty = NULL
  xb = az = uv = 0

  # Safety checks and parameter assignement
  u = v = a = b = NULL
  if (is.list(values)) {
    if (check(values$u, n, d)) u = values$u
    if (check(values$v, m, d)) v = values$v
    if (!is.null(Z)) if (check(values$a, n, ncol(Z))) a = values$a
    if (!is.null(X)) if (check(values$b, m, ncol(X))) b = values$b
  }

  # Compute the initial column-specific regression parameters (if any)
  if (!is.null(X)) {
    if (is.null(b)) {
      if (verbose) cat(" Initialization: column-specific covariates \n")
      xtx = crossprod(X)
      xty = crossprod(X, y)
      b = t(solve(xtx, xty))
    }
    xb = tcrossprod(X, b)
  }

  # Compute the initial row-specific regression parameter (if any)
  if (!is.null(Z)) {
    if (is.null(a)) {
      if (verbose) cat(" Initialization: row-specific covariates \n")
      ztz = crossprod(Z)
      zty = crossprod(Z, t(y - xb))
      a = t(solve(ztz, zty))
    }
    az = tcrossprod(a, Z)
  }

  # If both U and V are provided, orthogonalize them via incomplete SVD
  if (!is.null(u) & !is.null(v)) {
    uv = tcrossprod(u, v)
    s = svd::propack.svd(uv, neig = d)
    u = s$u %*% diag(sqrt(s$d))
    v = s$v %*% diag(sqrt(s$d))
  }

  # If U is unspecified, compute it via penalized least squares
  if (!is.null(u) & is.null(v)) {
    if (verbose) cat(" Initialization: latent scores \n")
    # Compute the unnormalized U via penalized least squares
    vtv = crossprod(v) + diag(d)
    vty = crossprod(v, t(y - xb - za))
    u = t(solve(vtv, vty))
    # Orthogonalize U and V via incomplete SVD
    uv = tcrossprod(u, v)
    s = svd::propack.svd(uv, neig = d)
    u = s$u %*% diag(sqrt(s$d))
    v = s$v %*% diag(sqrt(s$d))
  }

  # If V is unspecified, compute it via least squares
  if (is.null(u) & !is.null(v)) {
    if (verbose) cat(" Initialization: latent loadings \n")
    # Compute the unnormalized V via penalized least squares
    utu = crossprod(u)
    uty = crossprod(u, t(y - xb - za))
    v = t(solve(utu, uty))
    # Orthogonalize U and V via incomplete SVD
    uv = tcrossprod(u, v)
    s = svd::propack.svd(uv, neig = d)
    u = s$u %*% diag(sqrt(s$d))
    v = s$v %*% diag(sqrt(s$d))
  }

  # If both U and V are unspecified, compute them via incomplete SVD
  # calculated over the regression residuals
  if (is.null(u) & is.null(v)) {
    if (verbose) cat(" Initialization: latent scores and loadings \n")
    s = svd::propack.svd(y - xb - az, neig = d)
    u = s$u %*% diag(sqrt(s$d))
    v = s$v %*% diag(sqrt(s$d))
    uv = tcrossprod(u, v)
  }

  # Compute the initial vector of dispersion parameters
  mu = family$linkinv(xb + az + uv)
  var = family$variance(mu)
  phi = colMeans((Y - mu)^2 / var, na.rm = TRUE)

  # Regression matrices initialization when there are no covariates
  if (is.null(X)) bv = matrix(0, nrow = ncol(Y), ncol = 0)
  if (is.null(Z)) bu = matrix(0, nrow = nrow(Y), ncol = 0)

  # Return the obtained initial values
  list(u = u, v = v, bz = a, bx = b, phi = phi)
}

