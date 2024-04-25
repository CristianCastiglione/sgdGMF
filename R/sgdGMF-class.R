

#' @title S3 class: sgdgmf
#'
#' @description
#' A short description...
#'
#' @slot method estimation method to minimize the negative penalized log-likelihood
#' @slot family a \code{glm} family (see \code{\link{family}} for more details)
#' @slot ncomp rank of the latent matrix factorization
#' @slot npar number of parameters in the model
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
## @import methods
#' @export
setClass("sgdgmf",
  slots = list(
    method = "character",
    family = "list",
    ncomp = "numeric",
    npar = "numeric",
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
    sic = "numeric",
    exe.time = "vector",
    trace = "data.frame",
    summary.cv = "data.frame"
))

#' @title Refine the final estimate of a GMF model
#'
#' @description Refine the estimated latent scores of a GMF model via IRWLS
#'
#' @param object an object of class \code{sgdgmf}
#' @param normalize if \code{TRUE}, normalize \code{U} and \code{V} to uncorrelated Gaussian \code{U} and upper triangular \code{V} with positive diagonal
#' @param parallel if \code{TRUE}, use parallel computing using the \code{foreach} package
#' @param nthreads number of cores to be used in the \code{"glm"} method
#'
#' @method refit sgdgmf
#' @export
refit.sgdgmf = function (
    object,
    normalize = TRUE,
    verbose = FALSE,
    parallel = FALSE,
    nthreads = 1
) {

  # Default error message
  message = function (var)
    warning(paste0("Refit control: '", var,"' was set to default value."),
            call. = FALSE, immediate. = TRUE, domain = NULL)

  # Safety checks
  if (!is.logical(normalize)) {message("normalize"); normalize = TRUE}
  if (!is.logical(verbose)) {message("verbose"); verbose = FALSE}
  if (!is.logical(parallel)) {message("parallel"); paralle = FALSE}
  if (!is.numeric(nthreads) | nthreads < 1) {message("nthreads"); nthreads = 1}

  # Get the parameter dimensions
  idxA = seq(from = 1, to = ncol(object$A))
  idxU = seq(from = 1, to = ncol(object$U)) + ncol(object$A)

  # Fill the missing values with the predictions
  Y = object$Y
  if (anyNA(Y)) {
    isna = is.na(Y)
    Y[isna] = object$mu[isna]
  }

  # Refit A and U via IRWLS
  coefs = vglm.fit.coef(
    Y = t(Y), X = cbind(object$Z, object$V),
    family = object$family, offset = tcrossprod(object$B, object$X),
    parallel = parallel, nthreads = as.integer(nthreads), clust = NULL)

  # Set the final estimates
  object$A = coefs[, idxA]
  object$U = coefs[, idxU]

  # Recompute the linear predictor
  object$eta = tcrossprod(cbind(object$X, object$A, object$U),
                          cbind(object$B, object$Z, object$V))

  # Recompute the conditional mean and variance matrices
  object$mu = object$family$linkinv(object$eta)
  object$var = object$family$variance(object$mu)

  # Normalize the latent factors
  if (normalize) {
    uv = normalize.uv(object$U, object$V, method = "qr")
    object$U = uv$U
    object$V = uv$V
  }

  # Return the refined object
  return (object)
}

#' @title Compute the deviance of a GMF model
#'
#' @description Compute the deviance of an estimated GMF object
#'
#' @param object an object of class \code{sgdgmf}
#'
#' @method deviance sgdgmf
#' @export
deviance.sgdgmf = function (object, normalize = FALSE) {
  dev = matrix.deviance(object$mu, object$Y, object$family)
  if (normalize) {
    mu0 = matrix(mean(Y, na.rm = TRUE), nrow(Y), ncol(Y))
    dev0 = matrix.deviance(mu0, Y, object$family)
    dev = dev / dev0
  }
  return (dev)
}

#' @title Compute the AIC of a GMF model
#'
#' @description Compute the Akaike information criterion (AIC) of an estimated GMF object
#'
#' @param object an object of class \code{sgdgmf}
#'
#' @method AIC sgdgmf
#' @export
AIC.sgdgmf = function (object) {
  dev = deviance(object, normalize = FALSE)
  df = object$npar
  aic = dev + 2 * df
  return (aic)
}

#' @title Compute the BIC of a GMF model
#'
#' @description Compute the Bayesian information criterion (BIC) of an estimated GMF object
#'
#' @param object an object of class \code{sgdgmf}
#'
#' @method BIC sgdgmf
#' @export
BIC.sgdgmf = function (object) {
  dev = deviance(object, normalize = FALSE)
  df = object$npar
  nm = prod(dim(object$Y)) - sum(is.na(object$Y))
  bic = dev + df * log(nm)
  return (bic)
}

#' @title Compute the SIC of a GMF model
#'
#' @description Compute the Bayesian information criterion (SIC) of an estimated GMF object
#'
#' @param object an object of class \code{sgdgmf}
#'
#' @method SIC sgdgmf
#' @export
SIC.sgdgmf = function (object) {
  dev = deviance(object, normalize = FALSE)
  df = object$npar
  nm = prod(dim(object$Y)) - sum(is.na(object$Y))
  sic = dev + df * log(nm) / nm
  return (sic)
}

#' @title Print the fundamental characteristics of a GMF
#'
#' @description Print some summary information of a GMF model.
#'
#' @param object an object of class \code{sgdgmf}
#'
#' @method print sgdgmf
#' @export
print.sgdgmf = function (object) {
  # Percentage of explained deviance
  dev.fit  = matrix.deviance(object$mu, object$Y, object$family)
  dev.null = matrix.deviance(mean(object$Y, na.rm = TRUE), object$Y, object$family)
  dev.exp  = 100 * (1 - dev.fit / dev.null)

  # Elapsed execution time
  time.init = object$exe.time[1]
  time.opt = object$exe.time[2]
  time.tot = object$exe.time[3]

  # Print the output
  cat(gettextf("\n Number of samples: %d", nrow(object$Y)))
  cat(gettextf("\n Number of features: %d", ncol(object$Y)))
  cat(gettextf("\n Data sparsity: %.2f %%", 100 * mean(is.na(object$Y))))
  cat(gettextf("\n Column covariates: %d", ncol(object$X)))
  cat(gettextf("\n Row covariates: %d", ncol(object$Z)))
  cat(gettextf("\n Latent space rank: %d", object$ncomp))
  cat(gettextf("\n Number of parameters: %d", object$npar))
  cat(gettextf("\n Model family: %s", object$family$family))
  cat(gettextf("\n Model link: %s", object$family$link))
  cat(gettextf("\n Estimation method: %s", object$method))
  cat(gettextf("\n Explained deviance: %.2f %%", dev.exp))
  cat(gettextf("\n Initialization exe. time: %.2f s (%.2f m)", time.init, time.init/60))
  cat(gettextf("\n Optimization exe. time: %.2f s (%.2f m)", time.opt, time.opt/60))
  cat(gettextf("\n Total execution time: %.2f s (%.2f m)", time.tot, time.tot/60))
  cat("\n")
}

#' @title Extract the coefficient of a GMF model
#'
#' @description
#' Return the estimated coefficients of a GMF model, i.e., the row- and column-specific
#' regression effects, the latent scores and loadings.
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
    "colreg" = object$B,
    "rowreg" = object$A,
    "scores" = object$U,
    "loadings" = object$V,
    "all" = list(colcoef = object$B, rowcoef = object$A,
                 scores = object$U, loadings = object$V))
}

#' @rdname coefficients.sgdgmf
#' @method coef sgdgmf
#' @export
coef.sgdgmf = function (
    object, type = c("all", "colreg", "rowreg", "scores", "loadings")
) {
  coefficients.sgdgmf(object, type = type)
}

#' @title Extract the residuals of a GMF model
#'
#' @description
#' Extract the residuals of a GMF model and, if required, compute the eigenvalues
#' of the residuals covariance/correlation matrix.
#' Moreover, if required, return the partial residual of the model obtained by
#' excluding the matrix decomposition from the linear predictor.
#'
#' @param object an object of class \code{sgdgmf}
#' @param type the type of residuals which should be returned
#' @param partial if \code{TRUE}, computes the residuals excluding the matrix factorization from the linear predictor
#' @param normalize if \code{TRUE}, standardize the residuals column-by-column
#' @param fillna if \code{TRUE}, fills \code{NA} values column-by-column
#' @param spectrum if \code{TRUE}, returns the eigenvalues of the residual covariance matrix
#' @param ncomp number of eigenvalues to be calculated (only if \code{spectrum=TRUE})
#'
#' @details
#' Let \eqn{g(\mu) = \eta = X B^\top + \Gamma Z^\top + U V^\top} be the linear predictor of a
#' GMF model. Let \eqn{R = (r_{ij})} be the correspondent residual matrix.
#' The following residuals can be considered:
#' \itemize{
#' \item deviance: \eqn{r_{ij}^{_D} = \textrm{sign}(y_{ij} - \mu_{ij}) \sqrt{D(y_{ij}, \mu_{ij})}};
#' \item Pearson: \eqn{r_{ij}^{_P} = (y_{ij} - \mu_{ij}) / \sqrt{\nu(\mu_{ij})}};
#' \item working: \eqn{r_{ij}^{_W} = (y_{ij} - \mu_{ij}) / \{g'(\mu_{ij}) \,\nu(\mu_{ij})\}};
#' \item response: \eqn{r_{ij}^{_R} = y_{ij} - \mu_{ij}};
#' \item link: \eqn{r_{ij}^{_G} = g(y_{ij}) - \eta_{ij}}.
#' }
#' If \code{partial=TRUE}, \eqn{mu} is computed excluding the latent matrix decomposition
#' from the linear predictor, so as to obtain the partial residuals.
#'
#' Let \eqn{\Sigma} be the empirical variance-covariance matrix of \eqn{R}, being
#' \eqn{\sigma_{ij} = \textrm{Cov}(r_{:i}, r_{:j})}. Then, the latent spectrum of
#' the model is the collection of eigenvalues of \eqn{\Sigma}.
#'
#' Notice that, in case of Gaussian data, the latent spectrum corresponds to the principal
#' component analysis on the regression residuals, whose eigenvalues can be used to
#' infer the amount of variance explained by each principal component. Similarly,
#' we can use the (partial) latent spectrum in non-Gaussian data settings to infer
#' the correct number of principal components to include into the GMF model or to
#' detect some residual dependence structures not already explained by the model.
#'
#' @method residuals sgdgmf
#' @export
residuals.sgdgmf = function (
    object, type = c("deviance", "pearson", "working", "response", "link"),
    partial = FALSE, normalize = FALSE, fillna = FALSE, spectrum = FALSE, ncomp = 50
) {
  # Set the residual type
  type = match.arg(type)

  # Compute the predicted values
  y = object$Y
  family = object$family
  if (partial) {
    U = cbind(object$X, object$A)
    V = cbind(object$B, object$Z)
    eta = tcrossprod(U, V)
    mu = family$linkinv(eta)
  } else {
    eta = object$eta
    mu = object$mu
  }

  # Compute the residuals
  res = switch(type,
    "deviance" = sign(y - mu) * sqrt(abs(family$dev.resids(y, mu, 1))),
    "pearson" = (y - mu) / sqrt(abs(family$variance(mu))),
    "working" = (y - mu) * family$mu.eta(mu) / abs(family$variance(mu)),
    "response" = (y - mu),
    "link" = (family$transform(y) - eta))

  # Fill the missing values using Gaussian random values
  if (anyNA(res) & (fillna | spectrum)) {
    res = apply(res, 2, function (x) {
      if (anyNA(x)) {
        na = which(is.na(x) | is.nan(x))
        r = length(na)
        m = mean(x, na.rm = TRUE)
        s = sd(x, na.rm = TRUE)
        x[na] = rnorm(r, mean = m, sd = s)
      }
      return (x)
    })
  }

  # Standardize the residuals column-by-column
  if (normalize) {
    res = scale(res, center = TRUE, scale = TRUE)
  }

  # Decompose the residuals using incomplete SVD
  if (spectrum) {
    rcov = cov(res)
    ncomp = max(1, min(ncomp, ncol(res)))
    pca = RSpectra::eigs_sym(rcov, ncomp)

    # Estimate the explained and residual variance
    var.eig = pca$values
    var.tot = sum(diag(rcov))
    var.exp = sum(var.eig)
    var.res = var.tot - var.exp
  }

  # Return the residuals and the corresponding spectrum
  if (!spectrum) {
    return (res)
  } else {
    return (
      list(residuals = res,
           lambdas = var.eig,
           explained.var = var.exp,
           reminder.var = var.res,
           total.var = var.tot))
  }
}

#' @rdname residuals.sgdgmf
#' @method resid sgdgmf
#' @export
resid.sgdgmf = function (
    object, type = c("deviance", "pearson", "working", "response", "link"),
    partial = FALSE, normalize = FALSE, fillna = FALSE, spectrum = FALSE, ncomp = 50
) {
  residuals.sgdgmf(object, type = type, partial = partial, normalize = normalize,
                   fillna = fillna, spectrum = spectrum, ncomp = ncomp)
}


#' @title Extract the fitted values of a GMF models
#'
#' @description Computes the fitted values of a GMF model.
#'
#' @param object an object of class \code{sgdgmf}
#' @param type the type of fitted values which should be returned
#' @param partial if \code{TRUE}, returns the partial fitted values
#'
#' @method fitted sgdgmf
#' @export
fitted.sgdgmf = function (
    object, type = c("link", "response", "terms"), partial = FALSE
) {
  # Set the fitted value type
  type = match.arg(type)

  # Return the fitted values depending on the prediction type
  if (!partial) {
    switch(type,
      "link" = object$eta,
      "response" = object$mu,
      "terms" = list(
        XB = tcrossprod(object$X, object$B),
        AZ = tcrossprod(object$A, object$Z),
        UV = tcrossprod(object$U, object$V)))
  } else {
    XB = tcrossprod(object$X, object$B)
    AZ = tcrossprod(object$A, object$Z)
    switch(type,
      "link" = XB + AZ,
      "response" = object$family$linkinv(XB + AZ),
      "terms" = list(XB = XB, AZ = AZ, UV = NULL))
  }
}

#' @title Predict method for GMF models
#'
#' @description
#' Computes the predictions of a GMF model. Out-of-sample predictions for a new
#' set of responses and covariates are computed via MLE, by keeping fixed the values
#' of the estimated \code{B} and \code{V} and maximizing the likelihood with respect
#' to \code{A} and \code{U}.
#'
#' @param object an object of class \code{sgdgmf}
#' @param newY optionally, a matrix of new response variable
#' @param newX optionally, a matrix of new covariate values
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
      "response" = family$linkinv(tcrossprod(
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


#' @title Simulate method for GMF models
#'
#' @description
#' Simulate new data from a fitted generalized matrix factorization models
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
#' @method simulate sgdgmf
#' @export
simulate.sgdgmf = function (
    object, newY = NULL, newX = NULL,
    type = c("link", "response", "terms", "coef"),
    parallel = FALSE, nthreads = 1
) {
  type = match.arg(type)

  # ...
  # ...
  # ...
}

#' @title Spectrum method for GMF models
#'
#' @description
#' Compute the latent spectrum of a GMF model evaluated on the GLM residual scale.
#'
#' @param object an object of class \code{sgdgmf}
#' @param ncomp number of eigenvalues to compute
#' @param type the type of residual which should be returned
#' @param normalize if \code{TRUE}, standardize the residuals column-by-column
#'
#' @details
#' Let \eqn{g(\mu) = \eta = X B^\top + \Gamma Z^\top} be the linear predictor of a
#' GMF model where the latent factorization \eqn{U V^\top} is not included into the
#' model specification. Let \eqn{R} be the correspondent partial residual matrix.
#' Either deviance, \eqn{r_{ij} = \textrm{sign}(y_{ij} - \mu_{ij}) \sqrt{D(y_{ij}, \mu_{ij})}},
#' or Pearson, \eqn{r_{ij} = (y_{ij} - \mu_{ij}) / \sqrt{\nu(\mu_{ij})}}, residuals
#' can be considered. Finally, we define \eqn{\Sigma} as the empirical variance-covariance
#' matrix of \eqn{R}, being \eqn{\sigma_{ij} = \textrm{Cov}(r_{:i}, r_{:j})}. Then, we define
#' the latent spectrum of the model as the collection of eigenvalues of \eqn{\Sigma}.
#' Notice that, in case of Gaussian data, the latent spectrum corresponds to the principal
#' component analysis on the regression residuals, whose eigenvalues can be used to
#' infer how much signal can be explained by each principal component. Similarly,
#' we can use the latent spectrum in non-Gaussian data settings to infer the correct
#' number of principal components to include into the GMF model.
#'
#' @keywords internal
eigenval.sgdgmf = function (
    object, ncomp = object$ncomp,
    type = c("deviance", "pearson", "working", "link"),
    normalize = FALSE
) {
  # Set the type and stat parameters
  type = match.arg(type)
  stat = match.arg(stat)

  # Compute the model residuals
  family = object$family
  eta = tcrossprod(cbind(object$X, object$A), cbind(object$B, object$Z))
  mu = family$linkinv(eta)
  res = switch(type,
    "deviance" = sign(object$Y - mu) * sqrt(abs(family$dev.resid(object$Y, mu, 1))),
    "pearson" = (object$Y - mu) / sqrt(abs(family$variance(mu))),
    "working" = (object$Y - mu) * family$mu.eta(eta) / abs(family$variance(mu)),
    "link" = (family$transform(objec$Y) - eta))

  # Fill the missing values using Gaussian random values
  if (anyNA(res)) {
    res = apply(res, 2, function (x) {
      na = which(is.na(x) | is.nan(x))
      r = length(na)
      m = mean(x, na.rm = TRUE)
      s = sd(x, na.rm = TRUE)
      x[na] = rnorm(r, mean = m, sd = s)
      return (x)
    })
  }

  # Standardize the residuals column-by-column
  if (normalize) {
    res = scale(res, center = TRUE, scale = TRUE)
  }

  # Decompose the residuals using incomplete SVD
  pca = RSpectra::eigs_sym(cov(res), ncomp)

  # Estimate the explained and residual variance
  var.eig = pca$values
  var.tot = sum(diag(S))
  var.exp = sum(var.eig)
  var.res = var.tot - var.exp

  # Return the spectrum
  list(lambdas = var.eig, explained = var.exp,
       reminder = var.res, total = var.tot)
}

#' @title Plot diagnostics for a GMF model
#'
#' @description
#' Plots (one of) six diagnostics to graphically analyze the marginal and conditional
#' distribution of the residuals of a GMF model. Currently, the following plots are
#' available: residuals against observation indices, residuals agains fitted values,
#' absolute square-root residuals against fitted values, histogram of the residuals,
#' residual QQ-plot, residual ECDF-plot.
#'
#' @param object an object of class \code{sgdgmf}
#' @param type the type of plot which should be returned
#' @param resid the type of residuals which should be used
#' @param subsample if \code{TRUE}, computes the residuals over o small fraction of the data
#' @param sample.size the dimension of the sub-sample which should be used
#' @param partial if \code{TRUE}, computes the partial residuals
#' @param normalize if \code{TRUE}, standardizes the residuals column-by-column
#' @param fillna if \code{TRUE}, fills the \code{NA} values with \code{0}
#'
#' @method plot sgdgmf
#' @export
plot.sgdgmf = function (
    object,
    type = c("res-idx", "res-fit", "std-fit", "hist", "qq", "ecdf"),
    resid = c("deviance", "pearson", "working", "response", "link"),
    subsample = FALSE, sample.size = 500, partial = FALSE,
    normalize = FALSE, fillna = FALSE
) {
  type = match.arg(type)
  resid = match.arg(resid)

  # Get the fitted values
  fit = switch(resid,
    "deviance" = fitted(object, type = "response", partial = partial),
    "pearson" = fitted(object, type = "response", partial = partial),
    "working" = fitted(object, type = "response", partial = partial),
    "response" = fitted(object, type = "response", partial = partial),
    "link" = fitted(object, type = "link", partial = partial))

  # Get the residual values
  res = residuals(object, type = resid, partial = partial,
                  normalize = normalize, fillna = fillna, spectrum = FALSE)

  # Eventually, extract only a fraction of the data
  if (subsample) {
    n = nrow(object$Y)
    m = ncol(object$Y)
    if (sample.size < n*m) {
      idx = cbind(
        row = sample.int(n = n, size = sample.size, replace = TRUE),
        col = sample.int(n = m, size = sample.size, replace = TRUE))
      fit = fit[idx]
      res = res[idx]
    }
  }

  # Create the selected plot
  plt = switch(type,
    "res-idx" = {
      df = data.frame(residuals = c(res), index = c(1:prod(dim(res))))
      ggplot(data = df, map = aes(x = index, y = residuals)) +
      geom_point(alpha = 0.5) + geom_hline(yintercept = 0, col = 2, lty = 2) +
        labs(x = "Index", y = "Residuals", title = "Residuals vs Fitted values")
    },
    "res-fit" = {
      df = data.frame(residuals = c(res), fitted = c(fit))
      ggplot(data = df, map = aes(x = fitted, y = residuals)) +
        geom_point(alpha = 0.5) + geom_hline(yintercept = 0, col = 2, lty = 2) +
        labs(x = "Fitted values", y = "Residuals", title = "Residuals vs Fitted values")
    },
    "std-fit" = {
      res = sqrt(abs(res))
      df = data.frame(residuals = c(res), fitted = c(fit))
      ggplot(data = df, map = aes(x = fitted, y = residuals)) +
        geom_point(alpha = 0.5) + geom_hline(yintercept = mean(res), col = 2, lty = 2) +
        labs(x = "Fitted values", y = "|Residuals|", title = "Residuals vs Fitted values")
    },
    "hist" = {
      df = data.frame(residuals = c(res))
      ggplot(data = df, map = aes(x = residuals, y = after_stat(density))) +
        geom_histogram(bins = 30) + geom_vline(xintercept = 0, col = 2, lty = 2) +
        labs(x = "Residuals", y = "Frequency", title = "Histogram of the residuals")
    },
    "qq" = {
      df = list2DF(qqnorm(scale(c(res)), plot.it = FALSE))
      ggplot(data = df, map = aes(x = x, y = y)) +
        geom_abline(intercept = 0, slope = 1, color = 2, lty = 2) + geom_point(alpha = 0.5) +
        labs(x = "Theoretical quantiles", y = "Empirical quantiles", title = "Residual QQ-plot")
    },
    "ecdf" = {
      zn = scale(c(res))
      zz = seq(from = min(zn), to = max(zn), length = 100)
      df1 = data.frame(x = zn, y = ecdf(zn)(zn))
      df2 = data.frame(x = zz, y = pnorm(zz))
      ggplot() +
        geom_line(data = df2, map = aes(x = x, y = y), color = 2) +
        geom_point(data = df1, map = aes(x = x, y = y), alpha = 0.5) +
        labs(x = "Standardized residuals", y = "Empirical CDF", title = "Residual ECDF plot")
    })

  # Return the ggplot object
  return (plt)
}

#' @title Screeplot for the residuals of a GMF model
#'
#' @description
#' Plots the variances of the principal components of the residuals against the
#' number of principal component.
#'
#' @param object an object of class \code{sgdgmf}
#' @param ncomp number of components to be plotted
#' @param type the type of residuals which should be used
#' @param partial if \code{TRUE}, plots the eigenvalues of the partial residuals
#' @param normalize if \code{TRUE}, plots the eigenvalues of the standardized residuals
#' @param cumulative if \code{TRUE}, plots the cumulative sum of the eigenvalues
#' @param proportion if \code{TRUE}, plots the fractions of explained variance
#'
#' @method screeplot sgdgmf
#' @export
screeplot.sgdgmf = function (
    object, ncomp = 20,
    type = c("deviance", "pearson", "working", "response", "link"),
    partial = FALSE, normalize = FALSE,
    cumulative = FALSE, proportion = FALSE
) {

  # Compute the spctrum of the residuals
  ncomp = max(1, min(ncomp, nrow(object$U), nrow(object$V)))
  res = residuals(object, type = type, partial = partial,
                  normalize = normalize, fillna = TRUE,
                  spectrum = TRUE, ncomp = ncomp)

  # Eventually, normalize the eigenvalues
  lambdas = res$lambdas
  if (cumulative) lambdas = cumsum(lambdas)
  if (proportion) lambdas = lambdas / res$total.var

  # Draw the screeplot
  df = data.frame(components = 1:ncomp, lambdas = lambdas)
  plt = ggplot(data = df, map = aes(x = components, y = lambdas)) + geom_col() +
    labs(x = "Components", y = "Eigenvalues", title = "Residual screeplot")

  # Return the ggplot object
  return (plt)
}

#' @title Biplot of a GMF model
#'
#' @description
#' Plot the observations on a two-dimensional projection determined by the
#' estimated score matrix
#'
#' @param object an object of class \code{sgdgmf}
#' @param choices a length 2 vector specifying the components to plot
#' @param normalize if \code{TRUE}, orthogonalizes the scores using SVD
#' @param labels a vector of labels which should be plotted
#' @param palette the color-palette which should be used
#'
#' @method biplot sgdgmf
#' @export
biplot.sgdgmf = function (
    object, choices = 1:2, normalize = FALSE, labels = NULL, palette = NULL
) {
  # Get the data dimensions
  n = nrow(object$U)
  m = nrow(object$V)

  # Set the coordinates to draw
  i = max(1, min(choices[1], object$ncomp))
  j = max(1, min(choices[2], object$ncomp))

  # Set the point labels
  if (is.null(labels) | !is.list(labels)) {
    if (!is.list(labels)) labels = list()
    if (is.null(labels$scores)) labels$scores = c(1:n)
    if (is.null(labels$loadings)) labels$loadings = c(1:m)
  }

  # Normalize the columns of U and V
  if (normalize) {
    pca = RSpectra::svds(tcrossprod(object$U, object$V), object$ncomp)
    u = scale(pca$u[,c(i,j)])
    v = scale(pca$v[,c(i,j)])
  } else {
    u = scale(object$U[,c(i,j)])
    v = scale(object$V[,c(i,j)])
  }

  # Create the score and loading data-frames
  scores = data.frame(idx = labels$scores, pc1 = c(u[,1]), pc2 = c(u[,2]))
  loadings = data.frame(idx = labels$loadings, pc1 = c(v[,1]), pc2 = c(v[,2]))

  # Set the color palettes
  if (is.null(palette) | !is.list(palette)) {
    palette = list()
    palette$scores = viridisLite::viridis(100, option = "viridis")
    palette$loadings = viridisLite::viridis(100, option = "inferno")
  }
  if (is.list(palette)) {
    if (is.null(palette$scores)) palette$scores = viridisLite::viridis(100, option = "viridis")
    if (is.null(palette$loadings)) palette$loadings = viridisLite::viridis(100, option = "inferno")
  }

  # Create the score biplot
  plt.scores =
    ggplot(data = scores, map = aes(x = pc1, y = pc2, color = idx, label = idx)) +
    geom_hline(yintercept = 0, lty = 2, color = "grey40") +
    geom_vline(xintercept = 0, lty = 2, color = "grey40") +
    geom_point() + geom_text(color = 1, size = 2.5, nudge_x = -0.1, nudge_y = +0.1) +
    scale_colour_gradientn(colours = palette$scores) + theme(legend.position = "bottom") +
    labs(x = paste("PC", i), y = paste("PC", j), color = "Index", title = "Scores")

  # Create the loading biplot
  plt.loadings =
    ggplot(data = loadings, map = aes(x = pc1, y = pc2, color = idx, label = idx)) +
    geom_hline(yintercept = 0, lty = 2, color = "grey40") +
    geom_vline(xintercept = 0, lty = 2, color = "grey40") +
    geom_point() + geom_text(color = 1, size = 2.5, nudge_x = -0.1, nudge_y = +0.1) +
    scale_colour_gradientn(colours = palette$loadings) + theme(legend.position = "bottom") +
    labs(x = paste("PC", i), y = paste("PC", j), color = "Index", title = "Loadings")

  # Return the ggplot objects
  list(scores = plt.scores, loadings = plt.loadings)
}

#' @title Heatmap of a GMF model
#'
#' @description
#' Plots a heatmap of either the data, the fitted values, or the residual values
#' of a GMF model allowing for different types of transformations and normalizations.
#' Moreover, it also permits to plot the latent score and loading matrices.
#'
#' @param object an object of class \code{sgdgmf}
#' @param type the type of data/predictions/residuals which should be returned
#' @param resid if \code{TRUE}, plots the residual values
#' @param symmetric if \code{TRUE}, symmetrizes the color limits
#' @param transpose if \code{TRUE}, transposes the matrix before plotting it
#' @param limits the color limits which should be used
#' @param palette the color-palette which should be used
#'
#' @method image sgdgmf
#' @export
image.sgdgmf = function (
    object,
    type = c("data", "response", "link", "scores", "loadings", "deviance", "pearson", "working"),
    resid = FALSE, symmetric = FALSE, transpose = FALSE, limits = NULL, palette = NULL
) {
  type = match.arg(type)

  # Safety checks
  if (resid) {
    if (type == "data") stop("type='data' is not allowed with resid=TRUE", call. = FALSE)
    if (type == "scores") stop("type='scores' is not allowed with resid=TRUE", call. = FALSE)
    if (type == "loadings") stop("type='loadings' is not allowed with resid=TRUE", call. = FALSE)

    mat = residuals(object, type = type)
  } else {
    if (type == "deviance") stop("type='deviance' is not allowed with resid=FALSE", call. = FALSE)
    if (type == "pearson") stop("type='pearson' is not allowed with resid=FALSE", call. = FALSE)
    if (type == "working") stop("type='working' is not allowed with resid=FALSE", call. = FALSE)

    mat = switch(type,
                 "data" = object$Y,
                 "response" = object$mu,
                 "link" = object$eta,
                 "scores" = object$U,
                 "loadings" = object$V)
  }

  if (transpose) mat = t(mat)

  df = reshape2::melt(mat, varnames = c("sample", "variable"))

  if (is.null(limits)) {
    limits = range(df$value, na.rm = TRUE)
    if (symmetric) limits = c(-1,+1) * max(abs(limits))
  }

  if (is.null(palette)) {
    palette = viridisLite::viridis(100, option = "viridis")
    if (resid) palette = grDevices::hcl.colors(100, palette = "RdBu")
    if (symmetric) palette = grDevices::hcl.colors(100, palette = "RdBu")
  }

  plt = ggplot(data = df, map = aes(x = variable, y = sample, fill = value)) +
    geom_raster() + scale_fill_gradientn(colours = palette, limits = limits) +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    theme(legend.position = "bottom", panel.grid = element_blank()) +
    labs(x = "Variables", y = "Samples", fill = "Intensity")

  return (plt)
}


