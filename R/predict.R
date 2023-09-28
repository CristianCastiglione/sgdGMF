
#' @title Predict method for GMF fits
#'
#' @description ...
#'
#' @export
predict.sgdgmf = function (
    object,
    newdata = NULL,
    idx = NULL,
    type = c("link", "response", "terms"), ...) {

  # Check function for the index vector
  check.idx = function (idx, n) {
    flag = FALSE
    if (!is.null(idx)) {
      if (is.vector(idx)) {
        if (is.integer(idx)) {
          if (all(idx <= n) && all(idx >= 1)) {
            flag = TRUE
          }
        }
      }
    }
    return (flag)
  }

  # Check function for the design matrices
  check.mat = function (u, dims) {
    flag = FALSE
    if (!is.null(newdata$u)) {
      if (is.matrix(newdata$u)) {
        if (dim(newdata$u) == dims) {
          flag = TRUE
        }
      }
    }
    return (flag);
  }

  # Data and model dimensions
  n = object$dims$nrow
  m = object$dims$ncol
  d = object$dims$ncomp
  p = object$dims$nx
  q = object$dims$nz

  # safety check for the index vector
  if (check.idx(idx, n)) idx = 1:n

  # Completed factor and loading matrices
  U = cbind(object$data$X, object$coef$betaZ, object$coef$U)[idx, ]
  V = cbind(object$coef$betaX, object$data$Z, object$coef$V)

  # Overwrite the new-data matrices
  ni = length(idx)
  if (!is.null(newdata) & is.list(newdata)) {
    if (p>0) if (check.mat(newdata$X, c(ni, p))) U[,1:p] = newdata$X
    if (q>0) if (check.mat(newdata$Z, c(m, q))) V[,(p+1):(p+q)] = newdata$Z
    if (check.mat(newdata$U, c(ni, d))) U[,(p+q+1):(p+q+d)] = newdata$U
  }

  # Compute the linear predictor
  eta = NULL
  if (type == "link") {
    # Response scale predictions
    out = tcrossprod(U, V)
  } else if (type == "response") {
    # Linear predictor
    eta = tcrossprod(U, V)
    out = object$family$linkinv(eta)
  } else if (type == "terms") {
    # Response decomposition
    out = list(XB = tcrossprod(U, V[,1:p]),
               AZ = tcrossprod(U, V[,(p+1):(p+q)]),
               UV = tcrossprod(U, V[,(p+q+1):(p+q+d)]))
  }

  # Output
  return (out)
}

factors = function(object, newdata, idx) UseMethod("factors")
loadings = function(object) UseMethod("loadings")

#' @title Return the estimated/predicted factor matrix U
#'
#' @description ...
#'
#' @export
factors.sgdgmf = function (object, newdata = NULL, idx = NULL) {

  # Check function for the index vector
  check.idx = function (idx, n) {
    flag = FALSE
    if (!is.null(idx)) {
      if (is.vector(idx)) {
        if (is.integer(idx)) {
          if (all(idx <= n) && all(idx >= 1)) {
            flag = TRUE
          }
        }
      }
    }
    return (flag)
  }

  # Check function for the design matrices
  check.mat = function (u, dims) {
    flag = FALSE
    if (!is.null(newdata$u)) {
      if (is.matrix(newdata$u)) {
        if (dim(newdata$u) == dims) {
          flag = TRUE
        }
      }
    }
    return (flag);
  }

  # Data and model dimensions
  n = object$dims$nrow
  m = object$dims$ncol
  d = object$dims$ncomp
  p = object$dims$nx
  q = object$dims$nz

  # Check whether to return the estimated factors or
  # it is needed to compute an out-of-sample prediction
  out = NULL
  if (is.null(newdata)) {
    # Just return the estimated factors
    if (!check.idx(idx, n)) idx = 1:n
    out = object$cief$U[idx, ]
  } else {
    if (is.null(newdata$Y)) {
      stop("A response matrix is needed to estimate the latent factors.")
    } else {
      ny = nrow(newdata$Y)
      my = ncol(newdata$Y)
      X = Z = V = NULL
      if (p > 0) {
        flag = check.mat(newdata$X, c(ny, p))
        if (flag) {X = newdata$X} else {stop("X has the wrong dimensions.")}
      }
      if (q > 0) {
        flag = check.mat(newdata$Z, c(my, q))
        if (flag) {Z = newdata$Z} else {stop("Z has the wrong dimensions.")}
      }
      offset = tcrossprod(X, object$coef$betaX)
      V = cbind(Z, object$coef$V)
      U = matrix(0, nrow = ny, ncol = d)

      # HERE WE STILL HAVE TO IMPLEMENT THE CONDITIONAL ESTIMATION OF U GIVEN V
    }
  }
  return (out)
}

#' @title Return the estimated loading matrix V
#'
#' @description ...
#'
#' @export
loadings.sgdgmf = function (object) {
  return (object$coef$V)
}



