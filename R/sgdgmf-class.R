

#' @title S3 class: sgdgmf
#'
#' @description
#' A short description...
#'
#' @slot method estimation method to minimize the negative penalized log-likelihood
#' @slot family a \code{glm} family (see \code{\link{family}} for more details)
#' @slot ncomp rank of the latent matrix factorization
#' @slot control.init list of control parameters for the initialization (see \code{\link{set.control.init}} for more details)
#' @slot control.alg list of control parameters for the optimization (see \code{\link{set.control.alg}} for more details)
#' @slot control.cv list of control parameters for the cross-validation (see \code{\link{set.control.cv}} for more details)
#' @slot Y matrix of responses (\eqn{n \times m})
#' @slot X matrix of row fixed effects (\eqn{n \times p})
#' @slot Z matrix of column fixed effects (\eqn{m \times q})
#' @slot A matrix of row-specific regression effects (\eqn{n \times q})
#' @slot B matrix of column-specific regression effects (\eqn{m \times p})
#' @slot U matrix of latent scores (\eqn{n \times d})
#' @slot V matrix of factor loadings (\eqn{m \times d})
#' @slot eta matrix of linear predictor (\eqn{n \times m})
#' @slot mu matrix of fitted means (\eqn{n \times m})
#' @slot var matrix of fitted variances (\eqn{n \times m})
#' @slot phi scalar dispersion parameter
#' @slot penalty penalty  obtained at the end of the estimation
#' @slot deviance deviance obtained at the end of the estimation
#' @slot objective objective function at the end of the estimation
#' @slot aic Akaike information criterion of the estimated model
#' @slot bic Bayesian information criterion of the estimated model
#' @slot cbic Corrected Bayesian information criterion of the estimated model
#' @slot exe.time Final execution time in seconds
#' @slot trace data frame collecting the optimization history
#' @slot summary.cv data frame collecting the cross-validation history
#'
#' @import methods
#' @export
setClass("sgdgmf",
  slots = list(
    method = "character",
    family = "list",
    ncomp = "numeric",
    control.init = "list",
    control.alg = "list",
    control.cv = "list",
    Y = "matrix",
    X = "matrix",
    Z = "matrix",
    A = "matrix",
    B = "matrix",
    U = "matrix",
    V = "matrix",
    eta = "matrix",
    mu = "matrix",
    var = "matrix",
    phi = "numeric",
    penalty = "numeric",
    deviance = "numeric",
    objective = "numeric",
    aic = "numeric",
    bic = "numeric",
    cbic = "numeric",
    exe.time = "numeric",
    trace = "data.frame",
    summary.cv = "data.frame"
))

#' @title Extract the coefficient of a GMF model
#'
#' @description ...
#'
#' @param object an object of class \code{sgdgmf}
#' @param type the type of coefficients which should be returned
#'
#' @method coefficients sgdgmf
#' @export
coefficients.sgdgmf = function (
    object, type = c("all", "colreg", "rowreg", "scores", "loadings")
) {
  type = match.arg(type)
  switch(type,
    "all" = list(colcoef = object$B, rowcoef = object$A,
                 scores = object$U, loadings = object$V),
    "colreg" = object$B,
    "rowreg" = object$A,
    "scores" = object$U,
    "loadings" = object$V)
}

#' @title Extract the residuals of a GMF model
#'
#' @description ...
#'
#' @param object an object of class \code{sgdgmf}
#' @param type the type of residuals which should be returned
#'
#' @method residuals sgdgmf
#' @export
residuals.sgdgmf = function (
    object, type = c("deviance", "pearson", "working", "response")
) {

  type = match.arg(type)
  y = object$Y
  mu = object$mu
  family = object$family

  switch(type,
    "deviance" = sign(y - mu) * sqrt(family$dev.resids(y, mu, 1)),
    "pearson" = (y - mu) / sqrt(family$variance(mu)),
    "working" = (y - mu) / family$mu.eta(mu),
    "response" = y - mu)
}


#' @title Extract the fitted values of a GMF models
#'
#' @description
#' Computes the predictions from a fitted generalized matrix factorization model (GMF).
#'
#' @param object an object of class \code{sgdgmf}
#' @param type the type of fitted values which should be returned
#'
#' @method fitted sgdgmf
#' @export
fitted.sgdgmf = function (
    object, type = c("link", "response", "terms")
) {
  # Set the fitted value type
  type = match.arg(type)

  # Return the fitted values depending on the prediction type
  switch(type,
    "link" = object$eta,
    "response" = object$mu,
    "terms" = list(
      XB = tcrossprod(object$X, object$B),
      AZ = tcrossprod(object$A, object$Z),
      UV = tcrossprod(object$U, object$V)))
}

#' @title Predict method for GMF models
#'
#' @description
#' Computes the predictions from a fitted generalized matrix factorization model (GMF).
#'
#' @param object an object of class \code{sgdgmf}
#' @param newdata optionally, a list containing new values for \code{X} and \code{Z}
#' @param type the type of prediction which should be returned
#' @param parallel if \code{TRUE}, allows for parallel computing using the package \code{foreach}
#' @param nthreads number of cores to be used in parallel (only if \code{parallel=TRUE})
#'
#' @details
#' If \code{newY} and \code{newX} are omitted, the predictions are based on the data
#' used for the fit. In that case, the predictions corresponds to the fitted values.
#' If \code{newY} and \code{newX} are provided, a corresponding set of \code{A} and
#' \code{U} are estimated via maximum likelihood using the \code{glm.fit} function.
#' By doing so, \code{B} and \code{V} are kept fixed.
#'
#' @method predict sgdgmf
#' @export
predict.sgdgmf = function (
    object, newY = NULL, newX = NULL,
    type = c("link", "response", "terms", "coef"),
    parallel = FALSE, nthreads = 1
) {

  # Set the prediction type
  type = match.arg(type)

  # Set the model family
  family = object$family

  # Set the model dimensions
  n = nrow(object$Y)
  m = ncol(object$Y)
  p = ncol(object$X)
  q = ncol(object$Z)
  d = object$ncomp

  # Check whether in-sample or out-of-sample prediction is required
  if (is.null(newY) && is.null(newX)) {
    # Compute the fitted values depending on the prediction type
    pred = switch(type,
      "link" = object$eta,
      "response" = object$mu,
      "terms" = list(
        XB = tcrossprod(object$X, object$B),
        AZ = tcrossprod(object$A, object$Z),
        UV = tcrossprod(object$U, object$V)),
      "coef" = list(
        B = object$B, A = object$A,
        U = object$U, V = object$V))
  } else {
    # Check the input matrices
    if (!is.numeric(newY)) stop("Type error: 'newY' is not numeric.")
    if (!is.numeric(newX)) stop("Type error: 'newX' is not numeric.")
    if (!is.matrix(newY)) stop("Type error: 'newY' is not a matrix.")
    if (!is.matrix(newX)) stop("Type error: 'newX' is not a matrix.")

    # Check the dimensions of the input matrices
    ny = nrow(newY); my = ncol(newY)
    nx = nrow(newX); px = ncol(newX)

    if (my != m ) stop("Incompatible dimensions: 'newY' has wrong dimentions.")
    if (px != p ) stop("Incompatible dimensions: 'newX' has wrong dimentions.")
    if (nx != ny) stop("Incompatible dimensions: 'newX' has wrong dimentions.")

    # Check the parallelization settings
    if (!is.logical(parallel)) stop("'parallel' must be a logical values")
    if (!is.numeric(nthreads)) stop("'nthreads' mus be a positive integer value")
    if (floor(nthreads) < 1) stop("'nthreads' mus be a positive integer value")

    idA = 1:q
    idU = (q+1):(q+d)
    if (parallel) {
      # Register the clusters
      ncores = parallel::detectCores() - 1
      ncores = max(1, min(floor(nthreads), ncores))
      clust = parallel::makeCluster(ncores)
      doParallel::registerDoParallel(clust)

      # Estimate the latent scores independently and using a GLM fitting strategy
      # (the row-by-row estimation is performed in parallel)
      newV = cbind(object$Z, object$V)
      newO = tcrossprod(newX, object$B)
      newU = foreach(i = 1:ny, .combine = "rbind") %dopar% {
        yi = as.vector(newY[i,])
        oi = as.vector(newO[i,])
        fit = stats::glm.fit(x = newV, y = yi, family = family, offset = oi)
        fit$coefficients
      }

      # Close the connection to the clusters
      parallel::stopCluster(clust)
    } else {
      # Estimate the latent scores independently and a parallel GLM fitting strategy
      # (the row-by-row estimation is performed sequentially)
      newV = cbind(object$Z, object$V)
      newO = tcrossprod(newX, object$B)
      newU = foreach(i = 1:ny, .combine = "rbind") %do% {
        yi = as.vector(newY[i,])
        oi = as.vector(newO[i,])
        fit = stats::glm.fit(x = newV, y = yi, family = family, offset = oi)
        fit$coefficients
      }
    }

    # Split the fixed-effect parameters from the factor scores
    newA = newU[, idA, drop = FALSE]
    newU = newU[, idU, drop = FALSE]

    # Compute the out-of-sample values depending on the
    pred = switch (type,
      "link" = tcrossprod(
        cbind(newX, newA, newU),
        cbind(object$B, object$Z, object$V)),
      "response" = f$linkinv(tcrossprod(
        cbind(newX, newA, newU),
        cbind(object$B, object$Z, object$V))),
      "terms" = list(
        XB = tcrossprod(newX, object$B),
        AZ = tcrossprod(newA, object$Z),
        UV = tcrossprod(newU, object$V)),
      "coef" = list(
        B = object$B, A = newA,
        U = newU, V = object$V))
  }

  # Output
  return (pred)
}



