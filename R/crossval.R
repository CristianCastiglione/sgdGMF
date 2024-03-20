
#' @title Model selection via cross-validation for generalized matrix factorization models
#'
#' @description
#' ...
#'
#' @param Y matrix of responses (\eqn{n \times m})
#' @param X matrix of row fixed effects (\eqn{n \times p})
#' @param Z matrix of column fixed effects (\eqn{q \times m})
#' @param family a \code{glm} family (see \code{\link{family}} for more details)
#' @param ncomps rank of the latent matrix factorization (default 2)
#' @param weights an optional matrix of weights (\eqn{n \times m})
#' @param offset an optional matrix of offset values (\eqn{n \times m}), that specify a known component to be included in the linear predictor.
#' @param method estimation method to minimize the negative penalized log-likelihood
#' @param penalty list of penalty parameters (see \code{\link{set.penalty}} for more details)
#' @param control.init list of control parameters for the initialization (see \code{\link{set.control.init}} for more details)
#' @param control.alg list of control parameters for the optimization (see \code{\link{set.control.alg}} for more details)
#' @param control.cv list of control parameters for the cross-validation (see \code{\link{set.control.cv}} for more details)
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
#' @examples
#' ...
#'
#' @export
sgdgmf.cv = function (
    Y,
    X = NULL,
    Z = NULL,
    family = poisson(),
    ncomps = seq(from = 1, to = 10, by = 1),
    weights = NULL,
    offset = NULL,
    method = c("airwls", "newton", "msgd", "csgd", "rsgd", "csgd"),
    penalty = list(),
    control.init = list(),
    control.alg = list(),
    control.cv = list()
) {

  # Sanity check for the cross-validation options
  control.cv = do.call("set.control.cv", control.cv)

  # Data dimensions
  n = nrow(Y)
  m = ncol(Y)
  p = ifelse(is.null(X), 0, ncol(X))
  q = ifelse(is.null(Z), 0, ncol(Z))
  f = 0.3

  # Maximum rank
  ncomps = floor(ncomps)
  maxcomp = max(ncomps)
  nfolds = control.cv$nfolds
  parallel = control.cv$parallel
  nthreads = control.cv$nthreads

  # Perform model selection minimizing the test deviance
  {
    # Create a k-fold partition
    data = list()
    for (fold in 1:nfolds) {
      data[[fold]] = partition(Y, f)
    }

    # Initialize the summary data-frame
    cv = data.frame(ncomp = c(), fold = c(), df = c(),
                    aic = c(), bic = c(), cbic = c(), dev = c())

    # Estimation loop
    for (ncomp in ncomps) {

      # Effective number of parameters
      df = m * p + n * q + (n + m) * ncomp

      for (fold in 1:nfolds) {

        # Print the eestimation status
        cat(" Rank:", paste(ncomp, maxcomp, sep = "/"),
            " Fold:", paste(fold, nfolds, sep = "/"), "\n")

        # Estimated mean matrix
        mu = sgdgmf.fit(Y = data[[fold]]$train, X = X, Z = Z, family = family,
                        ncomp = ncomp, method = method, penalty = penalty,
                        control.init = control.init, control.alg = control.alg)$mu

        # Train and test sample sizes
        n.train = (1-f)*n*m
        n.test = f*n*m

        # Train and test goodness-of-fit measures
        dev.train = matrix.deviance(mu = mu, y = data[[fold]]$train, family = family)
        dev.test = matrix.deviance(mu = mu, y = data[[fold]]$test, family = family)
        aic.train = dev.train + 2 * df
        bic.train = dev.train + 2 * df * log(n.train)
        cbic.train = dev.train + 2 * df * log(log(n.train))

        # Summary data-frame
        cv = rbind(cv, data.frame(
            ncomp = ncomp, fold = fold, df = df,
            aic = aic.train / n.train, bic = bic.train / n.train,
            cbic = cbic.train / n.train, dev = dev.test / n.test
          )
        )
      }
    }

    # Rank selection
    if (nfolds > 1) {
      avgcv = data.frame(ncomp = c(), fold = c(), df = c(), aic = c(), bic = c(), dev = c())
      for (ncomp in 1:maxcomp) {
        for (fold in 1:nfolds) {
          idx = (cv$ncomp == ncomp) & (cv$fold == fold)
          df = mean(cv$df[idx])
          aic = mean(cv$aic[idx])
          bic = mean(cv$bic[idx])
          cbic = mean(cv$cbic[idx])
          dev = mean(cv$dev[idx])
          avgstat = data.frame(ncomp = ncomp, df = df, aic = aic, bic = bic, cbic, cbic, dev = dev)
          avgcv = rbind(avgcv, avgstat)
        }
      }
      ncomp = avgcv$ncomp[which.min(avgcv$dev)]
    } else {
      ncomp = cv$ncomp[which.min(cv$dev)]
    }
  }

  cat("Final refit with rank = ", ncomp, " \n")

  # Fit the model using the chosen optimizer
  fit = sgdgmf.fit(Y = Y, X = X, Z = Z, family = family,
                   ncomp = ncomp, method = method, penalty = penalty,
                   control.init = control.init, control.alg = control.alg)

  # Return the cross-validation parameters
  fit$control.cv = control.cv

  # Return the cross-validation statistics
  fit$summary.cv = cv

  # Return the fitted model
  return (fit)
}
