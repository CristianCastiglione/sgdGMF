
#' @title Model selection via cross-validation for generalized matrix factorization models
#'
#' @description
#' K-fold cross-validation for generalized matrix factorization (GMF) models.
#'
#' @param Y matrix of responses (\eqn{n \times m})
#' @param X matrix of row fixed effects (\eqn{n \times p})
#' @param Z matrix of column fixed effects (\eqn{q \times m})
#' @param family a \code{glm} family (see \code{\link{family}} for more details)
#' @param ncomps rank of the latent matrix factorization (default 2)
#' @param weights an optional matrix of weights (\eqn{n \times m})
#' @param offset an optional matrix of offset values (\eqn{n \times m}), that specify a known component to be included in the linear predictor.
#' @param method estimation method to minimize the negative penalized log-likelihood
#' @param sampling sub-sampling strategy to use if \code{method = "sgd"}
#' @param penalty list of penalty parameters (see \code{\link{set.penalty}} for more details)
#' @param control.init list of control parameters for the initialization (see \code{\link{set.control.init}} for more details)
#' @param control.alg list of control parameters for the optimization (see \code{\link{set.control.alg}} for more details)
#' @param control.cv list of control parameters for the cross-validation (see \code{\link{set.control.cv}} for more details)
#'
#' @return
#' If \code{refit = FALSE}, the function returns a list containing \code{control.init},
#' \code{control.alg}, \code{control.cv} and \code{summary.cv}. The latter is a matrix
#' collecting the cross-validation results for each combination of fold and latent
#' dimension.
#'
#' If \code{refit = TRUE}, the function returns an object of class \code{sgdgmf},
#' obtained by refitting the model on the whole data matrix using the latent dimension
#' selected via cross-validation. The returned object also contains the \code{summary.cv}
#' information along with the other standard output of the \code{\link{sgdgmf.fit}} function.
#'
#' @details
#' Cross-validation is performed by minimizing the estimated out-of-sample error, which
#' can be measured in terms of averaged deviance, AIC or BIC calculated on fold-specific
#' test sets. Within each fold, the test set is defined as a fixed proportion of entries
#' in the response matrix which are held out from the estimation process.
#' To this end, the test set entries are hidden by \code{NA} values when training the
#' model. Then, the predicted, i.e. imputed, values are used to compute the fold-specific
#' out-of-sample error.
#'
#' @example examples/example-crossval.R
#'
#' @export sgdgmf.cv
sgdgmf.cv = function (
    Y,
    X = NULL,
    Z = NULL,
    family = gaussian(),
    ncomps = seq(from = 1, to = 10, by = 1),
    weights = NULL,
    offset = NULL,
    method = c("airwls", "newton", "sgd"),
    sampling = c("block", "coord"),
    penalty = list(),
    control.init = list(),
    control.alg = list(),
    control.cv = list()
) {

  method = match.arg(method)

  # Sanity check for the cross-validation options
  control.cv = do.call("set.control.cv", control.cv)
  control.init = do.call("set.control.init", control.init)
  control.alg = set.control.alg(method, sampling, control.alg)

  # Data dimensions
  n = nrow(Y)
  m = ncol(Y)
  p = ifelse(is.null(X), 1, ncol(X)) # if X = NULL, sgdgmf.fit set X = rep(1, n)
  q = ifelse(is.null(Z), 1, ncol(Z)) # if X = NULL, sgdgmf.fit set Z = rep(1, m)

  # Maximum rank
  ncomps = floor(ncomps)
  maxcomp = max(ncomps)

  criterion = control.cv$criterion
  refit = control.cv$refit
  nfolds = control.cv$nfolds
  common = isTRUE(control.cv$init == "common")
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
    control.init$values = sgdgmf.init(
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
        train = train, test = test, X = X, Z = Z, family = family, ncomp = ncomp,
        maxcomp = maxcomp, fold = fold, nfolds = nfolds, method = method,
        sampling = sampling, penalty = penalty, control.init = control.init,
        control.alg = control.alg, control.cv = control.cv)
    }
  } else {
    # Parallel estimation loop
    # WARNING:
    # ... this "foreach" loop does not work yet!
    # ... The problem seems to be that the function "sgdgmf.cv.step"
    # ... is not visible from the "foreach" environment, as well as
    # ... all the other functions that are called by "sgdgmf.cv.step".
    # SOLUTION:
    # ... I guess that the problem can be solved by just building and installing
    # ... the "sgdGMF" package. By doing this, the "foreach" option .package = c("sgdGMF")
    # ... should work properly, including in the "foreach" environment
    # ... all the functions and utilities of the "sgdGMF" package
    cv = foreach (iter = 1:niter, .packages = c("sgdGMF"), .combine = "rbind") %dopar% {

      # Set the number of components and the group
      ncomp = groups$ncomp[iter]
      fold = groups$fold[iter]

      # Set the train and test data
      train = data[[fold]]$train
      test = data[[fold]]$test

      # Run the estimation and compute the GoF statistics
      sgdgmf.cv.step(
        train = train, test = test, X = X, Z = Z, family = family, ncomp = ncomp,
        maxcomp = maxcomp, fold = fold, nfolds = nfolds, method = method,
        sampling = sampling, penalty = penalty, control.init = control.init,
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
    for (ncomp in ncomps) {
      for (fold in 1:nfolds) {
        idx = (cv$ncomp == ncomp) & (cv$fold == fold)
        avgstat = data.frame(
          ncomp = ncomp,
          df = mean(cv$df[idx]),
          dev = mean(cv$dev[idx], na.rm = TRUE),
          mae = mean(cv$mae[idx], na.rm = TRUE),
          mse = mean(cv$mse[idx], na.rm = TRUE),
          aic = mean(cv$aic[idx], na.rm = TRUE),
          bic = mean(cv$bic[idx], na.rm = TRUE)
        )
        avgcv = rbind(avgcv, avgstat)
      }
    }
    ncomp = switch(criterion,
      "dev" = avgcv$ncomp[which.min(avgcv$dev)],
      "aic" = avgcv$ncomp[which.min(avgcv$aic)],
      "bic" = avgcv$ncomp[which.min(avgcv$bic)],
      "mae" = avgcv$ncomp[which.min(avgcv$mae)],
      "mse" = avgcv$ncomp[which.min(avgcv$mse)])
  } else {
    ncomp = switch(criterion,
      "dev" = cv$ncomp[which.min(cv$dev)],
      "aic" = cv$ncomp[which.min(cv$aic)],
      "bic" = cv$ncomp[which.min(cv$bic)],
      "mae" = cv$ncomp[which.min(cv$mae)],
      "mse" = cv$ncomp[which.min(cv$mse)])
  }

  if (refit) {
    # Re-fit the model using the selected optimizer
    cat("Final refit with rank =", ncomp, "\n")
    control.init$values$U = control.init$values$U[, 1:ncomp, drop = FALSE]
    control.init$values$V = control.init$values$V[, 1:ncomp, drop = FALSE]

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

#' @title Single step of cross-validation for generalized matrix factorization models
#'
#' @description
#' Internal function running a single step of cross-validation for generalized
#' matrix factorization (GMF) models and calculating some goodness-of-fit measures
#' on the train and test sets.
#'
#' @return
#' Returns a \code{data.frame}  containing the current number of latent factors
#' in the model (\code{ncomp}), the fold identifier (\code{fold}), the degrees of
#' freedom, i.e. the number of parameters, of the model (\code{df}), the AIC, BIC
#' and deviance (respectively, \code{aic}, \code{bic}, \code{dev})
#' calculated on the train and test sets.
#'
#' @keywords internal
sgdgmf.cv.step = function (
    train, test, X, Z, family, ncomp, maxcomp,
    fold, nfolds, method, sampling, penalty,
    control.init, control.alg, control.cv
) {
  # Set the model dimensions
  n = nrow(train)
  m = ncol(train)
  p = ifelse(is.null(X), 1, ncol(X))
  q = ifelse(is.null(Z), 1, ncol(Z))

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
    ncomp = ncomp, method = method, sampling = sampling, penalty = penalty,
    control.init = control.init, control.alg = control.alg)$mu

  # Train and test sample sizes
  n.train = n * m * (1-f)
  n.test = n * m * f

  # Train and test goodness-of-fit measures
  dev.train = sum(family$dev.resids(train, mu, 1), na.rm = TRUE)
  dev.test = sum(family$dev.resids(test, mu, 1), na.rm = TRUE)
  aic.train = dev.train + 2 * df
  bic.train = dev.train + df * log(n.train)
  mae.train = sum(abs(train - mu), na.rm = TRUE)
  mae.test = sum(abs(test - mu), na.rm = TRUE)
  mse.train = sum((train - mu)^2, na.rm = TRUE)
  mse.test = sum((test - mu)^2, na.rm = TRUE)

  # Return a data-frame with all the GoF statistics
  data.frame(
    ncomp = ncomp,
    fold = fold,
    df = df,
    train = n.train,
    test = n.test,
    aic = aic.train / n.train,
    bic = bic.train / n.train,
    dev = dev.test / n.test,
    mae = mae.test / n.test,
    mse = mse.test / n.test
  )
}


