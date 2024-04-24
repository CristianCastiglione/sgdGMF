
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
#' @export sgdgmf.cv
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
  control.init = do.call("set.control.init", control.init)
  control.alg = switch(method,
    "airwls" = do.call("set.control.airwls", control.alg),
    "newton" = do.call("set.control.newton", control.alg),
    "msgd" = do.call("set.control.msgd", control.alg),
    "csgd" = do.call("set.control.csgd", control.alg),
    "rsgd" = do.call("set.control.rsgd", control.alg),
    "bsgd" = do.call("set.control.bsgd", control.alg))
  control.cv = do.call("set.control.cv", control.cv)

  # Data dimensions
  n = nrow(Y)
  m = ncol(Y)
  p = ifelse(is.null(X), 0, ncol(X))
  q = ifelse(is.null(Z), 0, ncol(Z))

  # Maximum rank
  ncomps = floor(ncomps)
  maxcomp = max(ncomps)

  criterion = control.cv$criterion
  refit = control.cv$refit
  nfolds = control.cv$nfolds
  common = control.cv$init == "common"
  parallel = control.cv$parallel
  nthreads = control.cv$nthreads

  if (parallel) stop("Parallelization is still to be implemented.", call. = FALSE)

  # Check whether the algorithm and cross-validation controls are consistent
  if (method %in% c("airwls", "newton")) {
    if (control.cv$parallel && (control.init$parallel || control.alg$parallel)) {
      stop("Nested parallelization is not allowed.")
    }
  }

  if (control.alg$verbose) {
    control.cv$verbose = TRUE
  }

  if (control.cv$parallel && (control.init$verbose || control.alg$verbose || control.cv$verbose)) {
    control.init$verbose = FALSE
    control.alg$verbose = FALSE
    control.cv$verbose = FALSE
    warning("With 'parallel=TRUE', no output is printed, that is 'verbose=FALSE'.",
            call. = FALSE, immediate. = TRUE, domain = NULL)
  }

  # Common initialization
  if (common) {
    time.init = proc.time()
    control.init$values = init.gmf.param(
      Y = Y, X = X, Z = Z, ncomp = maxcomp,
      family = family, method = control.init$method,
      type = control.init$type, niter = control.init$niter,
      values = control.init$values,verbose = control.init$verbose,
      parallel = control.init$parallel, nthreads = control.init$threads)
    control.init$method = "values"
    time.init = as.numeric(proc.time() - time.init)[3]
  }

  # Create a k-fold partition
  data = list()
  f = control.cv$proportion
  for (fold in 1:nfolds) {
    data[[fold]] = partition(Y, f)
  }

  # Cartesian product of groups and folds
  groups = expand.grid(fold = 1:nfolds, ncomp = ncomps)
  niter = nrow(groups)

  # Register the clusters
  if (parallel) {
    ncores = parallel::detectCores() - 1
    ncores = max(1, min(nthreads, ncores))
    clust = parallel::makeCluster(ncores)
    parallel::makeExport(c("sgdGMF"))
    doParallel::registerDoParallel(clust)
  }

  if (!parallel) {
    # Sequential estimation loop
    cv = foreach (iter = 1:niter, .combine = "rbind") %do% {

      # Set the number of components and the group
      ncomp = groups$ncomp[iter]
      fold = groups$fold[iter]

      # Set the train and test data
      train = data[[fold]]$train
      test = data[[fold]]$test

      # Run the estimation and compute the GoF statistics
      sgdgmf.cv.step(
        train = train, test = test, X = X, Z = Z, family = family,
        ncomp = ncomp, maxcomp = maxcomp, fold = fold, nfolds = nfolds,
        method = method, penalty = penalty, control.init = control.init,
        control.alg = control.alg, control.cv = control.cv)
    }
  } else {
    # Parallel estimation loop
    # WARNING:
    # ... this "foreach" loop does not work yet!
    # ... The problem seems to be that the function "sgdgmf.cv.step"
    # ... is not visible from the "foreach" environment, as well as
    # ... all the other functions that are called by "sgdgmf.cv.step".
    cv = foreach (iter = 1:niter, .packages = c("sgdGMF"), .combine = "rbind") %dopar% {

      # Set the number of components and the group
      ncomp = groups$ncomp[iter]
      fold = groups$fold[iter]

      # Set the train and test data
      train = data[[fold]]$train
      test = data[[fold]]$test

      # Run the estimation and compute the GoF statistics
      sgdgmf.cv.step(
        train = train, test = test, X = X, Z = Z, family = family,
        ncomp = ncomp, maxcomp = maxcomp, fold = fold, nfolds = nfolds,
        method = method, penalty = penalty, control.init = control.init,
        control.alg = control.alg, control.cv = control.cv)
    }
  }

  # Close the connection to the clusters
  if (parallel) {
    parallel::stopCluster(clust)
  }

  # Rank selection
  if (nfolds > 1) {
    avgcv = data.frame()
    for (ncomp in 1:maxcomp) {
      for (fold in 1:nfolds) {
        idx = (cv$ncomp == ncomp) & (cv$fold == fold)
        avgstat = data.frame(
          ncomp = ncomp,
          df = mean(cv$df[idx]),
          dev = mean(cv$dev[idx], na.rm = TRUE),
          aic = mean(cv$aic[idx], na.rm = TRUE),
          bic = mean(cv$bic[idx], na.rm = TRUE),
          sic = mean(cv$sic[idx], na.rm = TRUE)
        )
        avgcv = rbind(avgcv, avgstat)
      }
    }
    ncomp = switch(criterion,
      "dev" = avgcv$ncomp[which.min(avgcv$dev)],
      "aic" = avgcv$ncomp[which.min(avgcv$aic)],
      "bic" = avgcv$ncomp[which.min(avgcv$bic)],
      "sic" = avgcv$ncomp[which.min(avgcv$sic)])
  } else {
    ncomp = switch(criterion,
      "dev" = cv$ncomp[which.min(cv$dev)],
      "aic" = cv$ncomp[which.min(cv$aic)],
      "bic" = cv$ncomp[which.min(cv$bic)],
      "sic" = cv$ncomp[which.min(cv$sic)])
  }

  if (refit) {
    # Re-fit the model using the chosen optimizer
    cat("Final refit with rank =", ncomp, "\n")
    fit = sgdgmf.fit(
      Y = Y, X = X, Z = Z, family = family,
      ncomp = ncomp, method = method, penalty = penalty,
      control.init = control.init, control.alg = control.alg)
    if (common) fit$exe.time[1] = time.init
  } else {
    # Do not re-fit the model, just return the summary statistics
    fit = list()
    fit$control.init = control.init
    fit$control.alg = control.alg
    fit$control.cv = control.cv
  }

  # Return the cross-validation parameters
  fit$control.cv = control.cv

  # Return the cross-validation statistics
  fit$summary.cv = cv

  # Return the fitted model
  return (fit)
}


#' @keywords internal
sgdgmf.cv.step = function (
    train, test, X, Z, family,
    ncomp, maxcomp, fold, nfolds, method, penalty,
    control.init, control.alg, control.cv
) {
  # Set the model dimensions
  n = nrow(train)
  m = ncol(train)
  p = ncol(X)
  q = ncol(Z)

  # Fraction of the data to use as a test set
  f = control.cv$proportion

  # Effective number of parameters
  df = m * p + n * q + (n + m) * ncomp

  # Print the estimation status
  if (control.cv$verbose)
    cat(" Rank:", paste(ncomp, maxcomp, sep = "/"),
        " Fold:", paste(fold, nfolds, sep = "/"), "\n")

  # Select the correct number of columns of the initial U and V matrices
  if (control.init$method == "values") {
    control.init$values$U = control.init$values$U[, 1:ncomp, drop = FALSE]
    control.init$values$V = control.init$values$V[, 1:ncomp, drop = FALSE]
  }

  # Estimated mean matrix
  mu = sgdgmf.fit(
    Y = train, X = X, Z = Z, family = family,
    ncomp = ncomp, method = method, penalty = penalty,
    control.init = control.init, control.alg = control.alg)$mu

  # Train and test sample sizes
  n.train = (1-f)*n*m
  n.test = f*n*m

  # Train and test goodness-of-fit measures
  dev.train = matrix.deviance(mu = mu, y = train, family = family)
  dev.test = matrix.deviance(mu = mu, y = test, family = family)
  aic.train = dev.train + 2 * df
  bic.train = dev.train + df * log(n.train)
  sic.train = dev.train + df * log(n.train) / n.train

  # Return a data-frame with all the GoF statistics
  data.frame(
    ncomp = ncomp, fold = fold, df = df,
    aic = aic.train / n.train, bic = bic.train / n.train,
    sic = sic.train / n.train, dev = dev.test / n.test
  )
}


