
#' @title Refine the final estimate of a GMF model
#'
#' @description Refine the estimated latent scores of a GMF model via IRWLS
#'
#' @param object an object of class \code{sgdgmf}
#' @param ... further arguments passed to or from other methods
#' @param normalize if \code{TRUE}, normalize \code{U} and \code{V} to uncorrelated Gaussian \code{U} and upper triangular \code{V} with positive diagonal
#' @param verbose if \code{TRUE}, print the optimization status
#' @param parallel if \code{TRUE}, use parallel computing using the \code{foreach} package
#' @param nthreads number of cores to be used in the \code{"glm"} method
#'
#' @returns An \code{sgdgmf} object containing the re-fitted model.
#'
#' @seealso \code{\link{sgdgmf.fit}}
#'
#' @examples
#' \donttest{# Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model using SGD
#' gmf_old = sgdgmf.fit(data$Y, ncomp = 3, family = poisson(), method = "sgd")
#'
#' # Refine the score matrix estimate
#' gmf_new = refit(gmf_old)
#'
#' # Get the fitted values in the link and response scales
#' mu_hat_old = fitted(gmf_old, type = "response")
#' mu_hat_new = fitted(gmf_new, type = "response")
#'
#' # Compare the results
#' oldpar = par(); par(mfrow = c(2,2), mar = c(1,1,3,1))
#' image(data$Y, axes = FALSE, main = expression(Y))
#' image(data$mu, axes = FALSE, main = expression(mu))
#' image(mu_hat_old, axes = FALSE, main = expression(hat(mu)[old]))
#' image(mu_hat_new, axes = FALSE, main = expression(hat(mu)[new]))
#' par(oldpar)
#' }
#' @method refit sgdgmf
#' @export
refit.sgdgmf = function (
    object, ...,
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
  if (!is.numeric(nthreads)) {message("nthreads"); nthreads = 1}

  # Set the number of threads
  ncores = parallel::detectCores() - 1
  nthreads = floor(max(1, min(nthreads, ncores)))

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
    family = object$family, weights = t(object$weights),
    offset = t(object$offset) + tcrossprod(object$B, object$X),
    parallel = parallel, nthreads = nthreads, clust = NULL)

  # Set the final estimates
  object$A = coefs[, idxA]
  object$U = coefs[, idxU]

  # Recompute the linear predictor
  object$eta = object$offset +
    tcrossprod(cbind(object$X, object$A, object$U),
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

#' @title Compute deviance, AIC and BIC of a GMF model
#'
#' @description Compute deviance, AIC and BIC of a GMF object
#'
#' @param object an object of class \code{sgdgmf}
#' @param ... further arguments passed to or from other methods
#' @param normalize if \code{TRUE}, normalize the result using the null-deviance
#' @param k the penalty parameter to be used for AIC; the default is \code{k = 2}
#'
#' @returns The value of the deviance extracted from a \code{sgdgmf} object.
#'
#' @examples
#' # Load the sgdGMF package
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
#' @method deviance sgdgmf
#' @export
deviance.sgdgmf = function (object, ..., normalize = FALSE) {
  dev = sum(object$family$dev.resids(object$Y, object$mu, 1), na.rm = TRUE)
  if (normalize) {
    n = nrow(object$Y)
    m = ncol(object$Y)
    mu0 = matrix(mean(object$Y, na.rm = TRUE), nrow = n, ncol = m)
    dev0 = sum(object$family$dev.resids(object$Y, mu0, 1), na.rm = TRUE)
    dev = dev / dev0
  }
  return (dev)
}

#' @rdname deviance.sgdgmf
#' @method AIC sgdgmf
#' @export
AIC.sgdgmf = function (object, ..., k = 2) {
  dev = deviance(object, normalize = FALSE)
  df = object$npar
  aic = dev + k * df
  return (aic)
}

#' @rdname deviance.sgdgmf
#' @method BIC sgdgmf
#' @export
BIC.sgdgmf = function (object, ...) {
  dev = deviance(object, normalize = FALSE)
  df = object$npar
  nm = prod(dim(object$Y)) - sum(is.na(object$Y))
  bic = dev + df * log(nm)
  return (bic)
}

#' @title Print the fundamental characteristics of a GMF
#'
#' @description Print some summary information of a GMF model.
#'
#' @param x an object of class \code{sgdgmf}
#' @param ... further arguments passed to or from other methods
#'
#' @returns No return value, called only for printing.
#'
#' @examples
#' \dontshow{
#' Sys.setenv(OPENBLAS_NUM_THREADS = 1)
#' Sys.setenv(MKL_NUM_THREADS = 1)
#' }# Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model with 3 latent factors
#' gmf = sgdgmf.fit(data$Y, ncomp = 3, family = poisson())
#'
#' # Print the GMF object
#' print(gmf)
#'
#' @method print sgdgmf
#' @export
print.sgdgmf = function (x, ...) {
  object = x

  # Percentage of explained deviance
  dev = 100 * (1 - deviance(object, normalize = TRUE))

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
  cat(gettextf("\n Explained deviance: %.2f %%", dev))
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
#' @param ... further arguments passed to or from other methods
#' @param type the type of coefficients which should be returned
#'
#' @return
#' If \code{type="all"}, a list of coefficients containing the fields \code{B}, \code{A}, \code{U} and \code{V}.
#' Otherwise, a matrix of coefficients, corresponding to the selected \code{type}.
#'
#' @examples
#' # Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model with 3 latent factors
#' gmf = sgdgmf.fit(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the estimated coefficients of a GMF model
#' str(coefficients(gmf)) # returns all the coefficients
#' str(coefficients(gmf, type = "scores")) # returns only the scores, say U
#' str(coefficients(gmf, type = "loadings")) # returns only the loadings, say V
#'
#' @method coefficients sgdgmf
#' @export
coefficients.sgdgmf = function (
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

#' @rdname coefficients.sgdgmf
#' @method coef sgdgmf
#' @export
coef.sgdgmf = function (
    object, ...,
    type = c("all", "colreg", "rowreg", "scores", "loadings")
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
#' @examples
#' # Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model with 3 latent factors
#' gmf = sgdgmf.fit(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the deviance residuals of a GMF model
#' str(residuals(gmf)) # returns the overall deviance residuals
#' str(residuals(gmf, partial = TRUE)) # returns the partial residuals
#' str(residuals(gmf, spectrum = TRUE)) # returns the eigenvalues of the residual var-cov matrix
#'
#' @method residuals sgdgmf
#' @export
residuals.sgdgmf = function (
    object, ...,
    type = c("deviance", "pearson", "working", "response", "link"),
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
    rcov = cov(res)
    ncomp = max(1, min(ncomp, ncol(res)))
    pca = NULL
    if (ncomp == ncol(res)) {
      pca = eigen(rcov)
    } else {
      pca = RSpectra::eigs_sym(rcov, ncomp)
    }

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
    object, ...,
    type = c("deviance", "pearson", "working", "response", "link"),
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
#' @param ... further arguments passed to or from other methods
#' @param type the type of fitted values which should be returned
#' @param partial if \code{TRUE}, returns the partial fitted values
#'
#' @return
#' If \code{type="terms"}, a list of fitted values containing the fields \code{XB},
#' \code{AZ} and \code{UV}. Otherwise, a matrix of fitted values in the link or
#' response scale, depending on the selected \code{type}.
#'
#' @examples
#' # Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model with 3 latent factors
#' gmf = sgdgmf.fit(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the fitted values of a GMF model
#' str(fitted(gmf)) # returns the overall fitted values in link scale
#' str(fitted(gmf, type = "response")) # returns the overall fitted values in response scale
#'
#' @method fitted sgdgmf
#' @export
fitted.sgdgmf = function (
    object, ...,
    type = c("link", "response", "terms"), partial = FALSE
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
#' @param ... further arguments passed to or from other methods
#' @param newY optionally, a matrix of new response variable
#' @param newX optionally, a matrix of new covariate values
#' @param type the type of prediction which should be returned
#' @param parallel if \code{TRUE}, allows for parallel computing using the package \code{foreach}
#' @param nthreads number of cores to be used in parallel (only if \code{parallel=TRUE})
#'
#' @return
#' If \code{type="link"} or \code{typr="response"}, a matrix of predictions.
#' If \code{type="terms"}, a list of matrices containing the fields \code{XB}, \code{AZ} and \code{UV}.
#' If \code{type="coef"}, a list of matrices containing the field \code{B}, \code{A}, \code{U} and \code{V}.
#'
#' @details
#' If \code{newY} and \code{newX} are omitted, the predictions are based on the data
#' used for the fit. In that case, the predictions corresponds to the fitted values.
#' If \code{newY} and \code{newX} are provided, a corresponding set of \code{A} and
#' \code{U} are estimated via maximum likelihood using the \code{glm.fit} function.
#' By doing so, \code{B} and \code{V} are kept fixed.
#'
#' @examples
#' # Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 120, m = 20, ncomp = 5, family = poisson())
#' train = sample(1:120, size = 100)
#' test = setdiff(1:120, train)
#'
#' Y = data$Y[train, ]
#' newY = data$Y[test, ]
#'
#' # Fit a GMF model with 3 latent factors
#' gmf = sgdgmf.fit(Y, ncomp = 3, family = poisson())
#'
#' # Get the fitted values of a GMF model
#' str(predict(gmf)) # returns the overall fitted values in link scale
#' str(predict(gmf, type = "response")) # returns the overall fitted values in response scale
#' str(predict(gmf, type = "terms")) # returns the partial fitted values in link scale
#' str(predict(gmf, newY = newY)) # returns the predictions for the new set of responses
#'
#' @method predict sgdgmf
#' @export
predict.sgdgmf = function (
    object, ...,
    newY = NULL, newX = NULL,
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
    if (!is.numeric(newY)) stop("Type error: 'newY' is not numeric.", call. = FALSE)
    if (!is.matrix(newY)) stop("Type error: 'newY' is not a matrix.", call. = FALSE)
    if (!is.null(newX) && !is.numeric(newX)) stop("Type error: 'newX' is not numeric.", call. = FALSE)
    if (!is.null(newX) && !is.matrix(newX)) stop("Type error: 'newX' is not a matrix.", call. = FALSE)
    if (is.null(newX)) {
      if (ncol(object$X) == 1 && stats::sd(object$X[,1]) == 0) {
        newX = matrix(object$X[1,1], nrow = nrow(newY), ncol = 1)
      } else {
        stop("Unspecified input: 'newX' must be provided.", call. = FALSE)
      }
    }

    # Check the dimensions of the input matrices
    ny = nrow(newY); my = ncol(newY)
    nx = nrow(newX); px = ncol(newX)

    if (my != m ) stop("Incompatible dimensions: 'newY' has wrong dimentions.", call. = FALSE)
    if (px != p ) stop("Incompatible dimensions: 'newX' has wrong dimentions.", call. = FALSE)
    if (nx != ny) stop("Incompatible dimensions: 'newX' has wrong dimentions.", call. = FALSE)

    # Check the parallelization settings
    if (!is.logical(parallel)) stop("'parallel' must be a logical values", call. = FALSE)
    if (!is.numeric(nthreads)) stop("'nthreads' mus be a positive integer value", call. = FALSE)
    if (floor(nthreads) < 1) stop("'nthreads' mus be a positive integer value", call. = FALSE)

    i = NULL
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
      # Estimate the latent scores independently with parallel GLM fitting strategy
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

  if (is.list(pred)) {
    pred = lapply(pred, function(p) {dimnames(p) = NULL; p})
  } else {
    dimnames(pred) = NULL
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
#' @param ... further arguments passed to or from other methods
#' @param nsim number of samples
#'
#' @returns An 3-fold array containing the simulated data.
#'
#' @examples
#' # Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model
#' gmf = sgdgmf.fit(data$Y, ncomp = 3, family = poisson())
#'
#' # Simulate new data from a GMF model
#' str(simulate(gmf))
#'
#' @method simulate sgdgmf
#' @export
simulate.sgdgmf = function (
    object, ..., nsim = 1
) {

  check = object$family$family %in% c("gaussian", "poisson", "binomial", "gamma")
  if (!check) {
    stop("Only the following families allow for data simulation: Gaussian, Poisson, Binomial, Gamma.")
  }

  # Data dimensions
  n = nrow(object$Y)
  m = ncol(object$Y)
  N = nsim * n * m

  # Estimated parameters
  mu = object$mu
  phi = object$phi
  wts = object$weights

  # Simulated values
  sim = switch(object$family$family,
    "gaussian" = stats::rnorm(N, mean = mu, sd = sqrt(phi / wts)),
    "poisson" = stats::rpois(N, lambda = mu),
    "binomial" = stats::rbinom(N, size = 1, prob = mu),
    "gamma" = stats::rgamma(N, shape = wts / phi, rate = wts / (phi * mu)),
    "invgaussian" = SuppDists::rinvGauss(N, nu = mu, lambda = wts / phi),
    "negbinom" = MASS::rnegbin(N, mu = mu, theta = phi / wts))

  # Reshaping
  sim = array(sim, dim = c(nsim, n, m))

  # Output
  return(sim)
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
#' @param x an object of class \code{sgdgmf}
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
#' @examples
#' \donttest{# Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model
#' gmf = sgdgmf.fit(data$Y, ncomp = 3, family = poisson())
#'
#' # Plot the residual-based GMF diagnostics
#' plot(gmf, type = "res-fit") # Residuals vs fitted values
#' plot(gmf, type = "std-fit") # Abs-sqrt-transformed residuals vs fitted values
#' plot(gmf, type = "qq") # Residual QQ-plot
#' plot(gmf, type = "hist") # Residual histogram
#' }
#' @method plot sgdgmf
#' @export
plot.sgdgmf = function (
    x, ...,
    type = c("res-idx", "res-fit", "std-fit", "hist", "qq", "ecdf"),
    resid = c("deviance", "pearson", "working", "response", "link"),
    subsample = FALSE, sample.size = 500, partial = FALSE,
    normalize = FALSE, fillna = FALSE
) {
  type = match.arg(type)
  resid = match.arg(resid)
  object = x

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
      ggplot(data = df, mapping = aes_string(x = "index", y = "residuals")) +
      geom_point(alpha = 0.5) + geom_hline(yintercept = 0, col = 2, lty = 2) +
        labs(x = "Index", y = "Residuals", title = "Residuals vs Indicies")
    },
    "res-fit" = {
      df = data.frame(residuals = c(res), fitted = c(fit))
      ggplot(data = df, mapping = aes_string(x = "fitted", y = "residuals")) +
        geom_point(alpha = 0.5) + geom_hline(yintercept = 0, col = 2, lty = 2) +
        labs(x = "Fitted values", y = "Residuals", title = "Residuals vs Fitted values")
    },
    "std-fit" = {
      res = sqrt(abs(res))
      df = data.frame(residuals = c(res), fitted = c(fit))
      ggplot(data = df, mapping = aes_string(x = "fitted", y = "residuals")) +
        geom_point(alpha = 0.5) + geom_hline(yintercept = mean(res), col = 2, lty = 2) +
        labs(x = "Fitted values", y = "|Residuals|", title = "Residuals vs Fitted values")
    },
    "hist" = {
      df = data.frame(residuals = c(res))
      ggplot(data = df, mapping = aes_string(x = "residuals", y = "after_stat(density)")) +
        geom_histogram(bins = 30) + geom_vline(xintercept = 0, col = 2, lty = 2) +
        labs(x = "Residuals", y = "Frequency", title = "Histogram of the residuals")
    },
    "qq" = {
      df = list2DF(stats::qqnorm(scale(c(res)), plot.it = FALSE))
      ggplot(data = df, mapping = aes_string(x = "x", y = "y")) +
        geom_abline(intercept = 0, slope = 1, color = 2, lty = 2) + geom_point(alpha = 0.5) +
        labs(x = "Theoretical quantiles", y = "Empirical quantiles", title = "Residual QQ-plot")
    },
    "ecdf" = {
      zn = scale(c(res))
      zz = seq(from = min(zn), to = max(zn), length = 100)
      df1 = data.frame(x = zn, y = stats::ecdf(zn)(zn))
      df2 = data.frame(x = zz, y = stats::pnorm(zz))
      ggplot() +
        geom_line(data = df2, mapping = aes_string(x = "x", y = "y"), color = 2) +
        geom_point(data = df1, mapping = aes_string(x = "x", y = "y"), alpha = 0.5) +
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
#' @examples
#' \donttest{# Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model
#' gmf = sgdgmf.fit(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the partial residual spectrum of a GMF model
#' screeplot(gmf) # screeplot of the var-cov matrix of the deviance residuals
#' screeplot(gmf, partial = TRUE) # screeplot of the partial residuals
#' screeplot(gmf, cumulative = TRUE) # cumulative screeplot
#' screeplot(gmf, proportion = TRUE) # proportion of explained residual variance
#' }
#' @method screeplot sgdgmf
#' @export
screeplot.sgdgmf = function (
    x, ...,
    ncomp = 20,
    type = c("deviance", "pearson", "working", "response", "link"),
    partial = FALSE, normalize = FALSE,
    cumulative = FALSE, proportion = FALSE
) {
  object = x

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
  plt = ggplot(data = df, mapping = aes_string(x = "components", y = "lambdas")) +
    geom_col() + labs(x = "Components", y = "Eigenvalues", title = "Residual screeplot")

  # Return the ggplot object
  return (plt)
}

#' @title Biplot of a GMF model
#'
#' @description
#' Plot the observations on a two-dimensional projection determined by the
#' estimated score matrix
#'
#' @param x an object of class \code{sgdgmf}
#' @param ... further arguments passed to or from other methods
#' @param choices a length 2 vector specifying the components to plot
#' @param arrange if \code{TRUE}, return a single plot with two panels
#' @param byrow if \code{TRUE}, the panels are arranged row-wise (if \code{arrange=TRUE})
#' @param normalize if \code{TRUE}, orthogonalizes the scores using SVD
#' @param labels a vector of labels which should be plotted
#' @param palette the color-palette which should be used
#' @param titles a 2-dimensional string vector containing the plot titles
#'
#' @return
#' If \code{arrange=TRUE}, a single ggplot object with the selected biplots,
#' otherwise, a list of two ggplot objects showing the row and column latent variables.
#'
#' @examples
#' \donttest{# Load the sgdGMF package
#' library(sgdGMF)
#'
#' # Generate data from a Poisson model
#' data = sim.gmf.data(n = 100, m = 20, ncomp = 5, family = poisson())
#'
#' # Fit a GMF model
#' gmf = sgdgmf.fit(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the biplot of a GMF model
#' biplot(gmf)
#' }
#' @method biplot sgdgmf
#' @export
biplot.sgdgmf = function (
    x, ...,
    choices = 1:2, arrange = TRUE, byrow = FALSE,
    normalize = FALSE, labels = NULL, palette = NULL,
    titles = c(NULL, NULL)
) {
  object = x

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

  colnames(scores) = c("idx", "pc1", "pc2")
  colnames(loadings) = c("idx", "pc1", "pc2")

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
  title.scores = ifelse(!is.null(titles[1]), titles[1], "Scores")
  plt.scores =
    ggplot(data = scores, mapping = aes_string(x = "pc1", y = "pc2", color = "idx", label = "idx")) +
    geom_hline(yintercept = 0, lty = 2, color = "grey40") +
    geom_vline(xintercept = 0, lty = 2, color = "grey40") +
    geom_point() + geom_text(color = 1, size = 2.5, nudge_x = -0.1, nudge_y = +0.1) +
    scale_colour_gradientn(colours = palette$scores) + theme(legend.position = "bottom") +
    labs(x = paste("PC", i), y = paste("PC", j), color = "Index", title = title.scores)

  # Create the loading biplot
  title.loadings = ifelse(!is.null(titles[1]), titles[1], "Loadings")
  plt.loadings =
    ggplot(data = loadings, mapping = aes_string(x = "pc1", y = "pc2", color = "idx", label = "idx")) +
    geom_hline(yintercept = 0, lty = 2, color = "grey40") +
    geom_vline(xintercept = 0, lty = 2, color = "grey40") +
    geom_point() + geom_text(color = 1, size = 2.5, nudge_x = -0.1, nudge_y = +0.1) +
    scale_colour_gradientn(colours = palette$loadings) + theme(legend.position = "bottom") +
    labs(x = paste("PC", i), y = paste("PC", j), color = "Index", title = title.loadings)

  plt = NULL
  if (arrange && byrow)
    plt = ggpubr::ggarrange(plt.scores, plt.loadings, nrow = 2, align = "v", legend = "right")
  if (arrange && !byrow)
    plt = ggpubr::ggarrange(plt.scores, plt.loadings, ncol = 2, align = "h", legend = "bottom")
  if (!arrange)
    plt = list(scores = plt.scores, loadings = plt.loadings)

  # Return the ggplot objects
  return(plt)
}

#' @title Heatmap of a GMF model
#'
#' @description
#' Plots a heatmap of either the data, the fitted values, or the residual values
#' of a GMF model allowing for different types of transformations and normalizations.
#' Moreover, it also permits to plot the latent score and loading matrices.
#'
#' @param x an object of class \code{sgdgmf}
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
#'  # Fit a GMF model
#'  gmf = sgdgmf.fit(data$Y, ncomp = 3, family = poisson())
#'
#' # Get the heatmap of a GMF model
#' image(gmf, type = "data") # original data
#' image(gmf, type = "response") # fitted values in response scale
#' image(gmf, type = "scores") # estimated score matrix
#' image(gmf, type = "loadings") # estimated loading matrix
#' image(gmf, type = "deviance", resid = TRUE) # deviance residual matrix
#' }
#' @method image sgdgmf
#' @export
image.sgdgmf = function (
    x, ...,
    type = c("data", "response", "link", "scores", "loadings", "deviance", "pearson", "working"),
    resid = FALSE, symmetric = FALSE, transpose = FALSE, limits = NULL, palette = NULL
) {
  type = match.arg(type)
  object = x

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

  plt = ggplot(data = df, mapping = aes_string(x = "variable", y = "sample", fill = "value")) +
    geom_raster() + scale_fill_gradientn(colours = palette, limits = limits) +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    theme(legend.position = "bottom", panel.grid = element_blank()) +
    labs(x = "Variables", y = "Samples", fill = "Intensity")

  return (plt)
}


