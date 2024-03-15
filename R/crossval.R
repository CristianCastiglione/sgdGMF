
#' @title Model selection via cross-validation for generalized matrix factorization models
#'
#' @description
#' ...
#'
#' @param Y matrix of responses (\eqn{n \times m})
#' @param X matrix of row fixed effects (\eqn{n \times p})
#' @param Z matrix of column fixed effects (\eqn{q \times m})
#' @param family a family as in the \code{\link{glm}} interface
#' @param ncomp number of random effects to estimate (default 2)
#' @param method optimization method: \code{"airwls"} (default),
#' \code{"newton"}, \code{"msgd"}, \code{"csgd"}, \code{"rsgd"}, \code{"bsgd"}
#' @param nfolds number of cross-validation folds (default 5)
#' @param parallel if \code{TRUE}, use parallel \code{foreach} to fit the models over different folds
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
sgdgmf.cv = function (
    Y,
    X = NULL,
    Z = NULL,
    family = poisson(),
    ncomp = 2,
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
  maxcomp = ncomp
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
                    aic = c(), bic = c(), dev = c())

    # Estimation loop
    for (ncomp in 1:maxcomp) {

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

        # Summary data-frame
        cv = rbind(cv,
          data.frame(
            ncomp = ncomp,
            fold = fold,
            df = df,
            aic = aic.train / n.train,
            bic = bic.train / n.train,
            dev = dev.test / n.test
          )
        )
      }
    }

    # Rank selection
    if (nfolds > 1) {
      avgcv = data.frame(ncomp = c(), fold = c(), df = c(), aic = c(), bic = c(), dev = c())
      for (ncomp in 1:maxcomp) {
        for (fold in 1:nfolds) {
          df = mean(cv$df[(cv$ncomp == ncomp) & (cv$fold == fold)])
          aic = mean(cv$aic[(cv$ncomp == ncomp) & (cv$fold == fold)])
          bic = mean(cv$bic[(cv$ncomp == ncomp) & (cv$fold == fold)])
          dev = mean(cv$dev[(cv$ncomp == ncomp) & (cv$fold == fold)])
          avgstat = data.frame(ncomp = ncomp, df = df, aic = aic, bic = bic, dev = dev)
          avgcv = rbind(avgcv, avgstat)
        }
      }
      ncomp = avgcv$ncomp[which.min(avgcv$dev)]
    } else {
      ncomp = vc$ncomp[which.min(cv$dev)]
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
  fit$cv.stat = cv

  # Return the fitted model
  return (fit)
}
