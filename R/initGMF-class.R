
#' @title S3 class: initgmf
#'
#' @description
#' A short description...
#'
#' @slot method initialization method to approximatelly minimize the negative penalized log-likelihood
#' @slot family a \code{glm} family (see \code{\link{family}} for more details)
#' @slot ncomp rank of the latent matrix factorization
#' @slot type residual type to be used for initializing the latent scores
#' @slot verbose if \code{TRUE}, print the optimization status (default \code{TRUE})
#' @slot parallel if \code{TRUE}, allows for parallel computing using the package \code{foreach} (only if \code{method="glm"})
#' @slot nthreads number of cores to be used in parallel (only if \code{parallel=TRUE})
#' @slot Y matrix of responses (\eqn{n \times m})
#' @slot X matrix of row fixed effects (\eqn{n \times p})
#' @slot Z matrix of column fixed effects (\eqn{m \times q})
#' @slot A matrix of row-specific regression effects (\eqn{n \times q})
#' @slot B matrix of column-specific regression effects (\eqn{m \times p})
#' @slot U matrix of latent scores (\eqn{n \times d})
#' @slot V matrix of factor loadings (\eqn{m \times d})
#' @slot phi scalar dispersion parameter
#'
#' @export
setClass("initgmf",
  slots = list(
    method = "character",
    family = "list",
    ncomp = "numeric",
    type = "character",
    verbose = "logical",
    parallel = "logical",
    nthreads = "numeric",
    savedata = "logical",
    Y = "matrix",
    X = "matrix",
    Z = "matrix",
    A = "matrix",
    B = "matrix",
    U = "matrix",
    V = "matrix",
    phi = "vector"
))

#' @method deviance initgmf
#' @export
deviance.initgmf = function (object, normalize = FALSE) {
  if (!object$savedata) {
    stop("'object' does not contain the data matrices Y, X and Z.", call. = FALSE)
  }
  U = cbind(object$X, object$A, object$U)
  V = cbind(object$B, object$Z, object$V)
  mu = object$family$linkinv(tcrossprod(U, V))
  dev = sum(object$family$dev.resids(object$Y, mu, 1), na.rm = TRUE)
  if (normalize) {
    n = nrow(object$Y)
    m = ncol(object$Y)
    mu0 = matrix(mean(object$Y, na.rm = TRUE), nrow = n, ncol = m)
    # mu0 = tcrossprod(rep(1, length = n), colMeans(object$Y, na.rm = TRUE))
    dev0 = sum(object$family$dev.resids(object$Y, mu0, 1), na.rm = TRUE)
    dev = dev / dev0
  }
  return (dev)
}

#' @method AIC initgmf
#' @export
AIC.initgmf = function (object) {
  if (!object$savedata) {
    stop("'object' does not contain the data matrices Y, X and Z.", call. = FALSE)
  }
  dev = deviance(object, normalize = FALSE)
  df = prod(dim(object$B)) + prod(dim(object$A)) + prod(dim(object$U)) + prod(dim(object$V))
  aic = dev + 2 * df
  return (aic)
}

#' @method BIC initgmf
#' @export
BIC.initgmf = function (object) {
  if (!object$savedata) {
    stop("'object' does not contain the data matrices Y, X and Z.", call. = FALSE)
  }
  dev = deviance(object, normalize = FALSE)
  df = prod(dim(object$B)) + prod(dim(object$A)) + prod(dim(object$U)) + prod(dim(object$V))
  nm = prod(dim(object$Y)) - sum(is.na(object$Y))
  bic = dev + df * log(nm)
  return (bic)
}

#' @method SIC initgmf
#' @export
SIC.initgmf = function (object) {
  if (!object$savedata) {
    stop("'object' does not contain the data matrices Y, X and Z.", call. = FALSE)
  }
  dev = deviance(object, normalize = FALSE)
  df = prod(dim(object$B)) + prod(dim(object$A)) + prod(dim(object$U)) + prod(dim(object$V))
  nm = prod(dim(object$Y)) - sum(is.na(object$Y))
  sic = dev + df * log(nm) / nm
  return (sic)
}

#' @title Extract the coefficient of an initialized GMF model
#'
#' @description
#' Return the initialized coefficients of a GMF model, i.e., the row- and column-specific
#' regression effects, the latent scores and loadings.
#'
#' @param object an object of class \code{initgmf}
#' @param type the type of coefficients which should be returned
#'
#' @seealso \code{\link{coefficients.sgdgmf}} and \code{\link{coef.sgdgmf}}.
#'
#' @method coefficients initgmf
#' @export
coefficients.initgmf = function (
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

#' @rdname coefficients.initgmf
#' @method coef initgmf
#' @export
coef.initgmf = function (
    object, type = c("all", "colreg", "rowreg", "scores", "loadings")
) {
  # Call the coefficients method for initgmf objects
  coefficients.initgmf(object, type = type)
}

#' @title Extract the residuals of an initialized GMF model
#'
#' @description
#' Extract the residuals of an initialized GMF model and, if required, compute
#' the eigenvalues of the residuals covariance/correlation matrix.
#' Moreover, if required, return the partial residual of the model obtained by
#' excluding the matrix decomposition from the linear predictor.
#'
#' @param object an object of class \code{initgmf}
#' @param type the type of residuals which should be returned
#' @param partial if \code{TRUE}, computes the residuals excluding the matrix factorization from the linear predictor
#' @param normalize if \code{TRUE}, standardize the residuals column-by-column
#' @param fillna if \code{TRUE}, fills \code{NA} values column-by-column
#' @param spectrum if \code{TRUE}, returns the eigenvalues of the residual covariance matrix
#' @param ncomp number of eigenvalues to be calculated (only if \code{spectrum=TRUE})
#'
#' @seealso \code{\link{residuals.sgdgmf}} and \code{\link{resid.sgdgmf}} for more details on the residual computation.
#'
#' @method residuals initgmf
#' @export
residuals.initgmf = function (
    object, type = c("deviance", "pearson", "working", "response", "link"),
    partial = FALSE, normalize = FALSE, fillna = FALSE, spectrum = FALSE, ncomp = 50
) {
  # Set the residual type
  type = match.arg(type)

  # Check if the data are available
  if (!object$savedata) {
    stop("'object' does not contain the data matrices Y, X and Z.", call. = FALSE)
  }

  # Compute the predicted values
  family = object$family
  if (partial) {
    U = cbind(object$X, object$A)
    V = cbind(object$B, object$Z)
  } else {
    U = cbind(object$X, object$A, object$U)
    V = cbind(object$B, object$Z, object$V)
  }
  eta = tcrossprod(U, V)
  mu = family$linkinv(eta)

  # Compute the residuals
  res = switch(type,
    "deviance" = sign(Y - mu) * sqrt(abs(family$dev.resids(Y, mu, 1))),
    "pearson" = (Y - mu) / sqrt(abs(family$variance(mu))),
    "working" = (Y - mu) * family$mu.eta(mu) / abs(family$variance(mu)),
    "response" = (Y - mu),
    "link" = (family$transform(Y) - eta))

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
    return (list(
      residuals = res,
      lambdas = var.eig, explained.var = var.exp,
      reminder.var = var.res, total.var = var.tot))
  }
}

#' @rdname residuals.initgmf
#' @method resid initgmf
#' @export
resid.initgmf = function (
    object, type = c("deviance", "pearson", "working", "response", "link"),
    partial = FALSE, normalize = FALSE, fillna = FALSE, spectrum = FALSE, ncomp = 50
) {
  # Call the residuals method for initgmf objects
  residuals.initgmf(object, type = type, partial = partial, normalize = normalize,
                    fillna = fillna, spectrum = spectrum, ncomp = ncomp)
}

#' @method fitted initgmf
#' @export
fitted.initgmf = function (
    object, type = c("link", "response", "terms"), partial = FALSE
) {
  # Set the fitted value type
  type = match.arg(type)

  # Check if the data are available
  if (!object$savedata) {
    stop("'object' does not contain the data matrices Y, X and Z.", call. = FALSE)
  }

  # Return the fitted values depending on the prediction type
  XB = tcrossprod(object$X, object$B)
  AZ = tcrossprod(object$A, object$Z)
  if (!partial) {
    UV = tcrossprod(object$U, object$V)
    switch(type,
      "link" = XB + AZ + UV,
      "response" = object$family$linkinv(XB + AZ + UV),
      "terms" = list(XB = XB, AZ = AZ, UV = UV))
  } else {
    switch(type,
      "link" = XB + AZ,
      "response" = object$family$linkinv(XB + AZ),
      "terms" = list(XB = XB, AZ = AZ, UV = NULL))
  }
}

#' @method plot initgmf
#' @export
plot.initgmf = function (
    object,
    type = c("res-idx", "res-fit", "std-fit", "hist", "qq", "ecdf"),
    resid = c("deviance", "pearson", "working", "response", "link"),
    subsample = FALSE, sample.size = 500, partial = FALSE,
    normalize = FALSE, fillna = FALSE
) {
  # Check if the data are available
  if (!object$savedata) {
    stop("'object' does not contain the data matrices Y, X and Z.", call. = FALSE)
  }

  # Call the plot method for sgdgmf objects
  plot.sgdgmf(object = object, type = type, resid = resid,
              subsample = subsample, sample.size = sample.size,
              partial = partial, normalize = normalize, fillna = fillna)
}

#' @method screeplot initgmf
#' @export
screeplot.initgmf = function (
    object, ncomp = 20,
    type = c("deviance", "pearson", "working", "response", "link"),
    partial = FALSE, normalize = FALSE,
    cumulative = FALSE, proportion = FALSE
) {

  # Check if the data are available
  if (!object$savedata) {
    stop("'object' does not contain the data matrices Y, X and Z.", call. = FALSE)
  }

  # Call the screeplot method for sgdgmf objects
  screeplot.sgdgmf(object = object, ncomp = ncomp, type = type,
                   partial = partial, normalize = normalize,
                   cumulative = cumulative, proportion = proportion)
}

#' @method biplot initgmf
#' @export
biplot.initgmf = function (
    object, choices = 1:2, normalize = FALSE, labels = NULL, palette = NULL
) {
  # Call the biplot method for sgdgmf objects
  biplot.sgdgmf(object = object, choices = choices,
                normalize = normalize, labels = labels, palette = palette)
}

#' @method image initgmf
#' @export
image.initgmf = function (
    object,
    type = c("data", "response", "link", "scores", "loadings", "deviance", "pearson", "working"),
    resid = FALSE, symmetric = FALSE, transpose = FALSE, limits = NULL, palette = NULL
) {
  type = match.arg(type)

  # Check if the data are available
  if (!object$savedata) {
    stop("'object' does not contain the data matrices Y, X and Z.", call. = FALSE)
  }

  # Store the predictions in the object
  object$eta = fitted.initgmf(object, type = "link")
  object$mu = object$family$linkinv(object$eta)

  # Call the image method for the sgdgmf object
  image.sgdgmf(object, type = type, resid = resid, symmetric = symmetric,
               transpose = transpose, limits = limits, palette = palette)
}

