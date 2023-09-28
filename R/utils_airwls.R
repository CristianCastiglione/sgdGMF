
#' @title GLM step
#' @description ...
#' @keywords internal
glm.step = function (X, Y,
                     beta,
                     family,
                     offset = 0,
                     stepsize = 1,
                     penalized = 0) {

  eta = offset + X %*% beta
  mu = family$linkinv(eta)

  mu.eta = family$mu.eta(eta)
  var.mu = family$variance(mu)

  Winv = var.mu / mu.eta**2
  W = mu.eta**2 / var.mu

  Z = (eta - offset) + (Y - mu) * Winv

  # Use only rows that give reasonable values
  thresh = 1e20
  keep.rows = !is.nan(c(W)) &
    !is.infinite(c(W)) &
    !is.nan(c(Z)) &
    !is.infinite(c(Z)) &
    (c(W) > 1/thresh) &
    (c(W) < thresh) &
    (c(abs(Z)) > 1/thresh) &
    (c(abs(Z)) < thresh)

  if (sum(keep.rows) < ncol(X)+1)
  {
    stop("Too many rows with infinite values")
  }
  Xfull = X[keep.rows,,drop=FALSE]
  Yfull = Z[keep.rows,,drop=FALSE]
  Wfull = c(W)[keep.rows]

  if (penalized) {
    Yfull = rbind(Yfull, matrix(0, ncol(Xfull), 1))
    Xfull = rbind(Xfull, diag(1, ncol(Xfull)))
    Wfull = c(Wfull, rep(penalized, ncol(Xfull)))
  }

  fit = suppressWarnings(lsfit(Xfull, Yfull, Wfull, intercept = FALSE))$coefficients

  return(stepsize * fit + (1 - stepsize) * beta)
}

#' @title GLM basic update
#' @description ...
#' @keywords internal
glm.basic = function (X, Y,
                      beta = NULL,
                      family = gaussian(),
                      tol = 1e-5,
                      offset = 0,
                      stepsize = 1,
                      penalized = 0,
                      steps = 1) {

  if (is.null(beta))
    beta = matrix(0, nrow = ncol(X), ncol = 1)
  for (i in 1:steps) {
    betaold = beta
    beta = glm.step(X, Y,
                    beta,
                    family,
                    offset = offset,
                    stepsize = stepsize,
                    penalized = penalized)
    if (norm_vec(betaold - beta) / norm_vec(beta) < tol) {
      break
    }
  }

  return(beta)
}

#' @title GLM update computed slice-wise
#' @description ...
#' @keywords internal
slice.glm = function(X, Y,
                     slices,
                     coefs,
                     offset = NULL,
                     penalized = 0,
                     parallel = 1,
                     family = poisson(),
                     method = "step",
                     stepsize = 1){
  res = c()
  steps = 10
  if (method == "step") {
    steps = 1
  }

  res = parallel::mclapply(1:slices, function(i){
    if (is.null(offset)){
      offset.slice = 0
    } else {
      offset.slice = offset[, i, drop = FALSE]
    }
    ok = !is.na(Y[,i])
    glm.basic(X[ok,  , drop = FALSE],
              Y[ok, i, drop = FALSE],
              beta = coefs[, i],
              family = family,
              offset = offset.slice,
              stepsize = stepsize,
              penalized = penalized,
              steps = steps)
  }, mc.cores = parallel)
  arr = simplify2array(res)
  if (is.vector(arr))
    arr = rbind(arr, NULL)

  return(arr)
}

