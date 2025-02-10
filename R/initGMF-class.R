
#' @title Compute deviance, AIC and BIC of an initialized GMF model
#'
#' @description Compute deviance, AIC and BIC of an initialized GMF object
#'
#' @param object an object of class \code{initgmf}
#' @param ... further arguments passed to or from other methods
#' @param normalize if \code{TRUE}, normalize the result using the null-deviance
#' @param k the penalty parameter to be used for AIC; the default is \code{k = 2}
#'
#' @returns The value of the deviance extracted from a \code{initgmf} object.
#'
#' @seealso \code{\link{deviance.sgdgmf}}, \code{\link{AIC.sgdgmf}} and \code{\link{AIC.sgdgmf}}.
#'
#' @examples
#' # Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model with 3 latent factors
#' init = sgdgmf.init(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the GMF deviance, AIC and BIC
#' deviance(init)
#' AIC(init)
#' BIC(init)
#'
#' @method deviance initgmf
#' @export
deviance.initgmf = function (object, ..., normalize = FALSE) {
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
AIC.initgmf = function (object, ..., k = 2) {
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
BIC.initgmf = function (object, ...) {
  if (!object$savedata) {
    stop("'object' does not contain the data matrices Y, X and Z.", call. = FALSE)
  }
  dev = deviance(object, normalize = FALSE)
  df = prod(dim(object$B)) + prod(dim(object$A)) + prod(dim(object$U)) + prod(dim(object$V))
  nm = prod(dim(object$Y)) - sum(is.na(object$Y))
  bic = dev + df * log(nm)
  return (bic)
}

#' @title Print the fundamental characteristics of an initialized GMF
#'
#' @description Print some summary information of an initialized GMF model.
#'
#' @param x an object of class \code{initgmf}
#' @param ... further arguments passed to or from other methods
#'
#' @returns No return value, called only for printing.
#'
#' @examples
#' # Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model with 3 latent factors
#' init = sgdgmf.init(data$Y, ncomp = 3, family = poisson())
#'
#' # Print the GMF object
#' print(init)
#'
#' @method print initgmf
#' @export
print.initgmf = function (x, ...) {

  object = x

  # Percentage of explained deviance
  dev = 100 * (1 - deviance(object, normalize = TRUE))

  # Elapsed execution time
  # time.init = object$exe.time[1]
  # time.opt = object$exe.time[2]
  # time.tot = object$exe.time[3]

  # Print the output
  cat(gettextf("\n Number of samples: %d", nrow(object$U)))
  cat(gettextf("\n Number of features: %d", nrow(object$V)))
  cat(gettextf("\n Data sparsity: %.2f %%", 100 * mean(is.na(object$Y))))
  cat(gettextf("\n Column covariates: %d", ncol(object$X)))
  cat(gettextf("\n Row covariates: %d", ncol(object$Z)))
  cat(gettextf("\n Latent space rank: %d", object$ncomp))
  cat(gettextf("\n Number of parameters: %d", object$npar))
  cat(gettextf("\n Model family: %s", object$family$family))
  cat(gettextf("\n Model link: %s", object$family$link))
  cat(gettextf("\n Estimation method: %s", object$method))
  cat(gettextf("\n Explained deviance: %.2f %%", dev))
  # cat(gettextf("\n Initialization exe. time: %.2f s (%.2f m)", time.init, time.init/60))
  # cat(gettextf("\n Optimization exe. time: %.2f s (%.2f m)", time.opt, time.opt/60))
  # cat(gettextf("\n Total execution time: %.2f s (%.2f m)", time.tot, time.tot/60))
  cat("\n")
}

#' @title Extract the coefficient of an initialized GMF model
#'
#' @description
#' Return the initialized coefficients of a GMF model, i.e., the row- and column-specific
#' regression effects, the latent scores and loadings.
#'
#' @param object an object of class \code{initgmf}
#' @param ... further arguments passed to or from other methods
#' @param type the type of coefficients which should be returned
#'
#' @return
#' If \code{type="all"}, a list of coefficients containing the fields \code{B}, \code{A}, \code{U} and \code{V}.
#' Otherwise, a matrix of coefficients, corresponding to the selected \code{type}.
#'
#' @seealso \code{\link{coefficients.sgdgmf}} and \code{\link{coef.sgdgmf}}.
#'
#' @examples
#' # Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model with 3 latent factors
#' init = sgdgmf.init(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the estimated coefficients of a GMF model
#' str(coefficients(init)) # returns all the coefficients
#' str(coefficients(init, type = "scores")) # returns only the scores, say U
#' str(coefficients(init, type = "loadings")) # returns only the loadings, say V
#'
#' @method coefficients initgmf
#' @export
coefficients.initgmf = function (
    object, ...,
    type = c("all", "colreg", "rowreg", "scores", "loadings")
) {
  type = match.arg(type)
  switch(type,
    "colreg" = object$B,
    "rowreg" = object$A,
    "scores" = object$U,
    "loadings" = object$V,
    "all" = list(B = object$B, A = object$A,
                 U = object$U, V = object$V))
}

#' @rdname coefficients.initgmf
#' @method coef initgmf
#' @export
coef.initgmf = function (
    object, ...,
    type = c("all", "colreg", "rowreg", "scores", "loadings")
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
#' @param ... further arguments passed to or from other methods
#' @param type the type of residuals which should be returned
#' @param partial if \code{TRUE}, computes the residuals excluding the matrix factorization from the linear predictor
#' @param normalize if \code{TRUE}, standardize the residuals column-by-column
#' @param fillna if \code{TRUE}, fills \code{NA} values column-by-column
#' @param spectrum if \code{TRUE}, returns the eigenvalues of the residual covariance matrix
#' @param ncomp number of eigenvalues to be calculated (only if \code{spectrum=TRUE})
#'
#' @return
#' If \code{spectrum=FALSE}, a matrix containing the selected residuals.
#' If \code{spectrum=TRUE}, a list containing the residuals (\code{res}), the first \code{ncomp}
#' eigenvalues of the residual covariance matrix, say (\code{lambdas}), the variance explained by the first
#' \code{ncomp} principal component of the residuals (\code{explained.var}), the variance not
#' explained by the first \code{ncomp} principal component of the residuals (\code{residual.var}),
#' the total variance of the residuals (\code{total.var}).
#'
#' @seealso \code{\link{residuals.sgdgmf}} and \code{\link{resid.sgdgmf}} for more details on the residual computation.
#'
#' @examples
#' # Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model with 3 latent factors
#' init = sgdgmf.init(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the deviance residuals of a GMF model
#' str(residuals(init)) # returns the overall deviance residuals
#' str(residuals(init, partial = TRUE)) # returns the partial residuals
#' str(residuals(init, spectrum = TRUE)) # returns the eigenvalues of the residual var-cov matrix
#'
#' @method residuals initgmf
#' @export
residuals.initgmf = function (
    object, ...,
    type = c("deviance", "pearson", "working", "response", "link"),
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
  Y = object$Y
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
        s = stats::sd(x, na.rm = TRUE)
        x[na] = stats::rnorm(r, mean = m, sd = s)
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
    rcov = stats::cov(res)
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
      residual.var = var.res, total.var = var.tot))
  }
}

#' @rdname residuals.initgmf
#' @method resid initgmf
#' @export
resid.initgmf = function (
    object, ...,
    type = c("deviance", "pearson", "working", "response", "link"),
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
#' @param ... further arguments passed to or from other methods
#' @param type the type of fitted values which should be returned
#' @param partial if \code{TRUE}, returns the partial fitted values
#'
#' @return
#' If \code{type="terms"}, a list of fitted values containing the fields \code{XB},
#' \code{AZ} and \code{UV}. Otherwise, a matrix of fitted values in the link or
#' response scale, depending on the selected \code{type}.
#'
#' @seealso \code{\link{fitted.sgdgmf}}.
#'
#' @examples
#' # Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model with 3 latent factors
#' init = sgdgmf.init(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the fitted values of a GMF model
#' str(fitted(init)) # returns the overall fitted values in link scale
#' str(fitted(init, type = "response")) # returns the overall fitted values in response scale
#' str(fitted(init, partial = TRUE)) # returns the partial fitted values in link scale
#'
#' @method fitted initgmf
#' @export
fitted.initgmf = function (
    object, ...,
    type = c("link", "response", "terms"), partial = FALSE
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
#' @param x an object of class \code{initgmf}
#' @param ... further arguments passed to or from other methods
#' @param type the type of plot which should be returned
#' @param resid the type of residuals which should be used
#' @param subsample if \code{TRUE}, computes the residuals over o small fraction of the data
#' @param sample.size the dimension of the sub-sample which should be used
#' @param partial if \code{TRUE}, computes the partial residuals
#' @param normalize if \code{TRUE}, standardizes the residuals column-by-column
#' @param fillna if \code{TRUE}, fills the \code{NA} values with \code{0}
#'
#' @returns A ggplot object showing the selected diagnostic plot.
#'
#' @seealso \code{\link{plot.sgdgmf}}.
#'
#' @examples
#' \donttest{# Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model
#' init = sgdgmf.init(data$Y, ncomp = 3, family = poisson())
#'
#' # Plot the residual-based GMF diagnostics
#' plot(init, type = "res-fit") # Residuals vs fitted values
#' plot(init, type = "std-fit") # Abs-sqrt-transformed residuals vs fitted values
#' plot(init, type = "qq") # Residual QQ-plot
#' plot(init, type = "hist") # Residual histogram
#' }
#' @method plot initgmf
#' @export
plot.initgmf = function (
    x, ...,
    type = c("res-idx", "res-fit", "std-fit", "hist", "qq", "ecdf"),
    resid = c("deviance", "pearson", "working", "response", "link"),
    subsample = FALSE, sample.size = 500, partial = FALSE,
    normalize = FALSE, fillna = FALSE
) {
  # Check if the data are available
  if (!x$savedata) {
    stop("'object' does not contain the data matrices Y, X and Z.", call. = FALSE)
  }

  # Call the plot method for sgdgmf objects
  plot.sgdgmf(x, type = type, resid = resid,
              subsample = subsample, sample.size = sample.size,
              partial = partial, normalize = normalize, fillna = fillna)
}

#' @title Screeplot for the residuals of an initialized GMF model
#'
#' @description
#' Plots the variances of the principal components of the residuals against the
#' number of principal component.
#'
#' @param x an object of class \code{sgdgmf}
#' @param ... further arguments passed to or from other methods
#' @param ncomp number of components to be plotted
#' @param type the type of residuals which should be used
#' @param partial if \code{TRUE}, plots the eigenvalues of the partial residuals
#' @param normalize if \code{TRUE}, plots the eigenvalues of the standardized residuals
#' @param cumulative if \code{TRUE}, plots the cumulative sum of the eigenvalues
#' @param proportion if \code{TRUE}, plots the fractions of explained variance
#'
#' @returns A ggplot object showing the residual screeplot of the model.
#'
#' @seealso \code{\link{screeplot.sgdgmf}}.
#'
#' @examples
#' \donttest{# Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model
#' init = sgdgmf.init(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the partial residual spectrum of a GMF model
#' screeplot(init) # screeplot of the var-cov matrix of the deviance residuals
#' screeplot(init, partial = TRUE) # screeplot of the partial residuals
#' screeplot(init, cumulative = TRUE) # cumulative screeplot
#' screeplot(init, proportion = TRUE) # proportion of explained residual variance
#' }
#' @method screeplot initgmf
#' @export
screeplot.initgmf = function (
    x, ...,
    ncomp = 20,
    type = c("deviance", "pearson", "working", "response", "link"),
    partial = FALSE, normalize = FALSE,
    cumulative = FALSE, proportion = FALSE
) {

  # Check if the data are available
  if (!x$savedata) {
    stop("'object' does not contain the data matrices Y, X and Z.", call. = FALSE)
  }

  # Call the screeplot method for sgdgmf objects
  screeplot.sgdgmf(x, ncomp = ncomp, type = type,
                   partial = partial, normalize = normalize,
                   cumulative = cumulative, proportion = proportion)
}


#' @title Biplot of an initialized GMF model
#'
#' @description
#' Plot the observations on a two-dimensional projection determined by the
#' estimated score matrix
#'
#' @param x an object of class \code{initgmf}
#' @param ... further arguments passed to or from other methods
#' @param choices a length 2 vector specifying the components to plot
#' @param arrange if \code{TRUE}, return a single plot with two panels
#' @param byrow if \code{TRUE}, the panels are arranged row-wise (if \code{arrange=TRUE})
#' @param normalize if \code{TRUE}, orthogonalizes the scores using SVD
#' @param labels a vector of labels which should be plotted
#' @param palette the color-palette which should be used
#'
#' @return
#' If \code{arrange=TRUE}, a single ggplot object with the selected biplots,
#' otherwise, a list of two ggplot objects showing the row and column latent variables.
#'
#' @seealso \code{\link{biplot.sgdgmf}}.
#'
#' @examples
#' \donttest{# Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model
#' init = sgdgmf.init(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the biplot of a GMF model
#' biplot(init) # 1st vs 2nd principal components
#' biplot(init, choices = 2:3) #2nd vs 3rd principal components
#' }
#' @method biplot initgmf
#' @export
biplot.initgmf = function (
    x, ...,
    choices = 1:2, arrange = TRUE, byrow = FALSE,
    normalize = FALSE, labels = NULL, palette = NULL
) {
  # Call the biplot method for sgdgmf objects
  biplot.sgdgmf(x, choices = choices, arrange = arrange,
                byrow = byrow, normalize = normalize, labels = labels,
                palette = palette)
}

#' @title Heatmap of an initialized GMF model
#'
#' @description
#' Plots a heatmap of either the data, the fitted values, or the residual values
#' of a GMF model allowing for different types of transformations and normalizations.
#' Moreover, it also permits to plot the latent score and loading matrices.
#'
#' @param x an object of class \code{initgmf}
#' @param ... further arguments passed to or from other methods
#' @param type the type of data/predictions/residuals which should be returned
#' @param resid if \code{TRUE}, plots the residual values
#' @param symmetric if \code{TRUE}, symmetrizes the color limits
#' @param transpose if \code{TRUE}, transposes the matrix before plotting it
#' @param limits the color limits which should be used
#' @param palette the color-palette which should be used
#'
#' @returns A ggplot object showing the selected heatmap.
#'
#' @examples
#' \donttest{# Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model
#' init = sgdgmf.init(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the heatmap of a GMF model
#' image(init, type = "data") # original data
#' image(init, type = "response") # fitted values in response scale
#' image(init, type = "scores") # estimated score matrix
#' image(init, type = "loadings") # estimated loading matrix
#' image(init, type = "deviance", resid = TRUE) # deviance residual matrix
#' }
#' @method image initgmf
#' @export
image.initgmf = function (
    x, ...,
    type = c("data", "response", "link", "scores", "loadings", "deviance", "pearson", "working"),
    resid = FALSE, symmetric = FALSE, transpose = FALSE, limits = NULL, palette = NULL
) {
  type = match.arg(type)

  # Check if the data are available
  if (!x$savedata) {
    stop("'object' does not contain the data matrices Y, X and Z.", call. = FALSE)
  }

  # Store the predictions in the object
  x$eta = fitted.initgmf(x, type = "link")
  x$mu = x$family$linkinv(x$eta)

  # Call the image method for the sgdgmf object
  image.sgdgmf(x, type = type, resid = resid, symmetric = symmetric,
               transpose = transpose, limits = limits, palette = palette)
}

