
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

#' @title Compute deviance, AIC and BIC of an initialized GMF model
#'
#' @description Compute deviance, AIC and BIC of an initialized GMF object
#'
#' @param object an object of class \code{initgmf}
#' @param normalize if \code{TRUE}, normalize the result using the null-deviance
#' @param k the penalty parameter to be used for AIC; the default is \code{k = 2}
#'
#' @seealso \code{\link{deviance.sgdgmf}}, \code{\link{AIC.sgdgmf}} and \code{\link{AIC.sgdgmf}}.
#'
#' @examples
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model with 3 latent factors
#' gmf = sgdgmf.fit(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the GMF deviance, AIC and BIC
#' deviance(gmf)
#' AIC(gmf)
#' BIC(gmf)
#'
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
    dev0 = sum(object$family$dev.resids(object$Y, mu0, 1), na.rm = TRUE)
    dev = dev / dev0
  }
  return (dev)
}

#' @rdname deviance.initgmf
#' @method AIC initgmf
#' @export
AIC.initgmf = function (object, k = 2) {
  if (!object$savedata) {
    stop("'object' does not contain the data matrices Y, X and Z.", call. = FALSE)
  }
  dev = deviance(object, normalize = FALSE)
  df = prod(dim(object$B)) + prod(dim(object$A)) + prod(dim(object$U)) + prod(dim(object$V))
  aic = dev + k * df
  return (aic)
}

#' @rdname deviance.initgmf
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
#' @examples
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model with 3 latent factors
#' init = init.gmf.param(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the estimated coefficients of a GMF model
#' coefficients(init) # returns all the coefficients
#' coefficients(init, type = "scores") # returns only the scores, say U
#' coefficients(init, type = "loadings") # returns only the loadings, say V
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
#' @examples
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model with 3 latent factors
#' init = init.gmf.param(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the deviance residuals of a GMF model
#' residuals(init) # returns the overall deviance residuals
#' residuals(init, partial = TRUE) # returns the partial residuals
#' residuals(init, spectrum = TRUE) # returns the eigenvalues of the residual var-cov matrix
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

#' @title Extract the fitted values of an initialized GMF model
#'
#' @description Computes the fitted values of an initialized GMF model.
#'
#' @param object an object of class \code{initgmf}
#' @param type the type of fitted values which should be returned
#' @param partial if \code{TRUE}, returns the partial fitted values
#'
#' @seealso \code{\link{fitted.sgdgmf}}.
#'
#' @examples
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model with 3 latent factors
#' init = init.gmf.param(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the fitted values of a GMF model
#' fitted(init) # returns the overall fitted values in link scale
#' fitted(init, type = "response") # returns the overall fitted values in response scale
#' fitted(init, partial = TRUE) # returns the partial fitted values in link scale
#'
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

#' @title Plot diagnostics for an initialized GMF model
#'
#' @description
#' Plots (one of) six diagnostics to graphically analyze the marginal and conditional
#' distribution of the residuals of a GMF model. Currently, the following plots are
#' available: residuals against observation indices, residuals agains fitted values,
#' absolute square-root residuals against fitted values, histogram of the residuals,
#' residual QQ-plot, residual ECDF-plot.
#'
#' @param object an object of class \code{initgmf}
#' @param type the type of plot which should be returned
#' @param resid the type of residuals which should be used
#' @param subsample if \code{TRUE}, computes the residuals over o small fraction of the data
#' @param sample.size the dimension of the sub-sample which should be used
#' @param partial if \code{TRUE}, computes the partial residuals
#' @param normalize if \code{TRUE}, standardizes the residuals column-by-column
#' @param fillna if \code{TRUE}, fills the \code{NA} values with \code{0}
#'
#' @seealso \code{\link{plot.sgdgmf}}.
#'
#' @examples
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model
#' init = init.gmf.param(data$Y, ncomp = 3, family = poisson())
#'
#' # Plot the residual-based GMF diagnostics
#' plot(init, type = "res-fit") # Residuals vs fitted values
#' plot(init, type = "std-fit") # Abs-sqrt-transformed residuals vs fitted values
#' plot(init, type = "qq") # Residual QQ-plot
#' plot(init, type = "hist") # Residual histogram
#'
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

#' @title Screeplot for the residuals of an initialized GMF model
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
#' @seealso \code{\link{screeplot.sgdgmf}}.
#'
#' @examples
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model
#' init = init.gmf.param(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the partial residual spectrum of a GMF model
#' screeplot(init) # screeplot of the var-cov matrix of the deviance residuals
#' screeplot(init, partial = TRUE) # screeplot of the partial residuals
#' screeplot(init, cumulative = TRUE) # cumulative screeplot
#' screeplot(init, proportion = TRUE) # proportion of explained residual variance
#'
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


#' @title Biplot of an initialized GMF model
#'
#' @description
#' Plot the observations on a two-dimensional projection determined by the
#' estimated score matrix
#'
#' @param object an object of class \code{initgmf}
#' @param choices a length 2 vector specifying the components to plot
#' @param normalize if \code{TRUE}, orthogonalizes the scores using SVD
#' @param labels a vector of labels which should be plotted
#' @param palette the color-palette which should be used
#'
#' @seealso \code{\link{biplot.sgdgmf}}.
#'
#' @examples
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model
#' init = init.gmf.param(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the biplot of a GMF model
#' biplot(init) # 1st vs 2nd principal components
#' biplot(init, choices = 2:3) #2nd vs 3rd principal components
#'
#' @method biplot initgmf
#' @export
biplot.initgmf = function (
    object, choices = 1:2, normalize = FALSE, labels = NULL, palette = NULL
) {
  # Call the biplot method for sgdgmf objects
  biplot.sgdgmf(object = object, choices = choices,
                normalize = normalize, labels = labels, palette = palette)
}

#' @title Heatmap of an initialized GMF model
#'
#' @description
#' Plots a heatmap of either the data, the fitted values, or the residual values
#' of a GMF model allowing for different types of transformations and normalizations.
#' Moreover, it also permits to plot the latent score and loading matrices.
#'
#' @param object an object of class \code{initgmf}
#' @param type the type of data/predictions/residuals which should be returned
#' @param resid if \code{TRUE}, plots the residual values
#' @param symmetric if \code{TRUE}, symmetrizes the color limits
#' @param transpose if \code{TRUE}, transposes the matrix before plotting it
#' @param limits the color limits which should be used
#' @param palette the color-palette which should be used
#'
#' @examples
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model
#' init = init.gmf.param(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the heatmap of a GMF model
#' image(init, type = "data") # original data
#' image(init, type = "response") # fitted values in response scale
#' image(init, type = "scores") # estimated score matrix
#' image(init, type = "loadings") # estimated loading matrix
#' image(init, type = "deviance", resid = TRUE) # deviance residual matrix
#'
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

