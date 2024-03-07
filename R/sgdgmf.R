

#' @title Factorize a matrix of non-Gaussian observations
#'
#' @description
#' ...
#'
#' @param Y matrix of responses (\eqn{n \times m})
#' @param X matrix of row fixed effects (\eqn{n \times p})
#' @param Z matrix of column fixed effects (\eqn{q \times m})
#' @param family a family as in the \code{\link{glm}} interface
#' @param ncomp number of random effects to estimate (default 2)
#' @param selection ...
#' @param penalty list of penalty parameters corresponding to u, v, a and b
#' @param init list of initialization options
#' @param control list of optimization options
#'
#' @return
#' ...
#'
#' @details
#' ...
#'
#' @references
#' ...
#'
#' @importFrom svd propack.svd
#'
#' @examples
#' ...
#'
#' @export
sgdgmf = function (
    Y,
    X = NULL,
    Z = NULL,
    family = poisson(),
    ncomp = 2,
    method = c("airwls", "newton", "msgd", "csgd", "rsgd", "csgd"),
    selection = FALSE,
    penalty = list(),
    init = list(),
    control = list()
) {

  # Sanity check for the selection parameter
  if (!is.logical(selection)) {
    stop("Selection must be a logical parameter.")
  }

  # Data dimensions
  n = nrow(Y)
  m = ncol(Y)
  p = ifelse(is.null(X), 0, ncol(X))
  q = ifelse(is.null(Z), 0, ncol(Z))
  f = 0.3

  # Perform model selection minimizing the test deviance
  if (selection) {
    maxcomp = ncomp
    data = partition(Y, f)
    dev = numeric(maxcomp)
    aic = numeric(maxcomp)
    bic = numeric(maxcomp)
    dfs = numeric(maxcomp)

    # Estimation loop
    for (ncomp in 1:maxcomp) {
      cat("Matrix rank:", ncomp, "\n")

      # Effective number of parameters
      df = m * p + n * q + (n + m) * ncomp

      # Estimated mean matrix
      mu = sgdgmf.fit(Y = data$train, X = X, Z = Z, family = family, ncomp = ncomp,
        method = method, penalty = penalty, init = init, control = control)$mu

      # Train and test deviances
      dev.train = matrix.deviance(mu = mu, y = data$train, family = family)
      dev.test = matrix.deviance(mu = mu, y = data$test, family = family)

      # Out-of-sample error and information criteria
      dfs[ncomp] = df
      dev[ncomp] = dev.test / (f*n*m)
      aic[ncomp] = (dev.train + 2 * df) / ((1-f)*n*m)
      bic[ncomp] = (dev.train + 2 * df * log((1-f)*n*m)) / ((1-f)*n*m)
    }

    cat("Final estimation \n")

    # Rank selection
    ncomp = which.min(dev)+1
  }

  # Fit the model using the chosen optimizer
  fit = sgdgmf.fit(Y = Y, X = X, Z = Z, family = family, ncomp = ncomp,
    method = method, penalty = penalty, init = init, control = control)

  if (selection) {
    fit$selection = data.frame(ncomp = 1:maxcomp, dev = dev, aic = aic, bic = bic, df = dfs)
  }

  # Convert the result in a S3 class
  class(fit) = "sgdgmf"

  # Return the fitted model
  return (fit)
}


