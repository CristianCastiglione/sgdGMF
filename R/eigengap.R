

#' @title Rank selection via eigenvalue-gap methods
#'
#' @description
#' Select the number of significant principal components of a matrix via
#' exploitation of eigenvalue-gap methods
#'
#' @param Y matrix of responses (\eqn{n \times m})
#' @param X matrix of row-specific fixed effects (\eqn{n \times p})
#' @param Z matrix of column-specific fixed effects (\eqn{q \times m})
#' @param maxcomp maximum number of eigenvalues to compute
#' @param family a family as in the \code{\link{glm}} interface (default \code{gaussian()})
#' @param weights matrix of optional weights (\eqn{n \times m})
#' @param offset matrix of optional offsets (\eqn{n \times m})
#' @param method rank selection method
#' @param type.reg regression method to be used to profile out the covariate effects
#' @param type.res residual type to be decomposed
#' @param maxiter maximum number of iterations
#' @param parallel if \code{TRUE}, allows for parallel computing using \code{foreach}
#' @param nthreads number of cores to be used in parallel (only if \code{parallel=TRUE})
#'
#' @references
#' Onatski, A. (2010).
#' Determining the number of factors from empirical distribution of eigenvalues.
#' Review of Economics and Statistics, 92(4): 1004-1016
#'
#' Fan, J., Guo, j. and Zheng, S. (2020)
#' Estimating number of factors by adjusted eigenvalues thresholding
#' Journal of the American Statistical Association, 117(538): 852--861
#'
#' Wang, L. and Carvalho, L. (2023)
#' Deviance matrix factorization
#' Electronic Journal of Statistics, 17(2): 3762-3810
#'
#' @export select.rank
select.rank = function (
    Y,
    X = NULL,
    Z = NULL,
    maxcomp = ncol(Y),
    family = gaussian(),
    weights = NULL,
    offset = NULL,
    method = c("onatski", "act"),
    type.reg = c("ols", "glm"),
    type.res = c("deviance", "pearson", "working", "link"),
    maxiter = 10,
    parallel = FALSE,
    nthreads = 1
) {
  # Set the selection method
  method = match.arg(method)
  type.reg = match.arg(type.reg)
  type.res = match.arg(type.res)

  # Set the family-specific data transformation
  family = set.family(family)

  # Fill the missing values
  Y[] = apply(Y, 2, function (x) {
    x[is.na(x)] = mean(x, na.rm = TRUE)
    return (x)
  })

  # Compute the transformed data
  gY = family$transform(Y)

  # Set the covariate matrices
  if (is.null(X)) X = matrix(1, nrow = n, ncol = 1)
  if (is.null(Z)) Z = matrix(1, nrow = n, ncol = 1)

  # Register and open the connection to the clusters
  clust = NULL
  if (parallel) {
    ncores = parallel::detectCores() - 1
    ncores = max(1, min(nthreads, ncores))
    clust = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(clust)
  }

  # Set the required GLM matrices
  eta = matrix(NA, nrow = n, ncol = m)
  mu  = matrix(NA, nrow = n, ncol = m)
  res = matrix(NA, nrow = n, ncol = m)

  # Initialize the regression coefficients
  B = switch(type.reg,
    "ols" = ols.fit.coef(gY, X, family = family, offset = NULL),
    "glm" = vglm.fit.coef(Y, X, family = family, offset = NULL,
                          parallel = parallel, nthreads = nthreads, clust = clust))
  #
  eta[] = tcrossprod(X, B)
  A = switch(type.reg,
    "ols" = ols.fit.coef(t(gY), Z, family = family, offset = t(eta)),
    "glm" = vglm.fit.coef(t(Y), Z, family = family, offset = t(eta),
                          parallel = parallel, nthreads = nthreads, clust = clust))

  # Close the connection to the clusters
  if (parallel) parallel::stopCluster(clust)

  ## # Initialize the regression coefficients
  ## eta = matrix(NA, nrow = n, ncol = m)
  ## mu  = matrix(NA, nrow = n, ncol = m)
  ## res = matrix(NA, nrow = n, ncol = m)
  ## B = foreach(j = 1:m, .combine = "rbind") %do% {
  ##   yj = as.vector(Y[,j])
  ##   fit = stats::glm.fit(x = X, y = yj, family = family)
  ##   t(fit$coefficients)
  ## }
  ## #
  ## eta[] = tcrossprod(X, B)
  ## A = foreach(i = 1:n, .combine = "rbind") %do% {
  ##   yi = as.vector(Y[i,])
  ##   oi = as.vector(eta[i,])
  ##   fit = stats::glm.fit(x = Z, y = yi, family = family, offset = oi)
  ##   fit$coefficients
  ## }

  # Initialize the linear predictor and the conditional mean matrix
  eta[] = eta + tcrossprod(A, Z)
  mu[] = family$linkinv(eta)
  res[] = switch(type.res,
    "deviance" = sign(Y - mu) * sqrt(abs(family$dev.resids(Y, mu, 1))),
    "pearson" = (Y - mu) / sqrt(abs(family$variance(mu))),
    "working" = (Y - mu) / abs(family$mu.eta(eta)),
    "link" = (gY - eta))

  # Select the optimal rank
  ncomp = switch(method,
    "onatski" = select.rank.onatski(res, maxcomp, maxiter)$ncomp,
    "act" = select.rank.act(res, maxcomp))

  # Return the selected rank
  return (ncomp)
}


#' @title Rank selection via the Onatski method
#'
#' @description
#' Select the number of significant principal components of a matrix via the
#' Onatski method
#'
#' @param Y matrix to be decomposed
#' @param maxcomp maximum number of eigenvalues to compute
#' @param maxiter maximum number of iterations
#'
#' @references
#' Onatski, A. (2010).
#' Determining the number of factors from empirical distribution of eigenvalues.
#' Review of Economics and Statistics, 92(4): 1004-1016
#'
#' @keywords internal
select.rank.onatski = function (Y, maxcomp = 50, maxiter = 100) {

  # Set the matrix dimension
  n = nrow(Y)
  m = ncol(Y)

  # Safety check for the number of maximum components
  if (maxcomp > m - 5) {
    maxcomp = m - 6
    warning("Rank selection: 'maxcomp' set to default value.",
            call. = FALSE, immediate. = TRUE, domain = NULL)
  }

  # Compute the spectrum of the covariance matrix of Y
  lambdas = eigen(cov(Y))$values

  # Initialize the loop parameters
  tol = 1e+03
  iter = 0
  ju = maxcomp + 1
  ncomp = maxcomp

  # Onatski selection loop
  while (tol >= 1) {
    j = ju
    yreg = lambdas[j:(j+4)]
    xreg = (j + (-1:3))^(2/3)
    delta = 2 * abs(cov(yreg, xreg)) / var(xreg)
    flag = which(-diff(lambdas) >= delta)
    ncomp = ifelse(length(flag) == 0, 0, tail(flag[flag <= maxcomp], 1))
    ju = ncomp + 1
    tol = abs(ju - j)
    iter = iter + 1
    if (iter > maxiter) break
  }

  # Check if the search was successful
  success = iter < maxiter

  # Return the selected rank
  list(ncomp = ncomp, delta = delta,
       niter = iter, success = success)
}


#' @title Rank selection via adjust correlation thresholding
#'
#' @description
#' Select the number of significant principal components of a matrix via adjust
#' correlation threshold (ACT)
#'
#' @param Y matrix to be decomposed
#' @param maxcomp maximum number of eigenvalues to compute
#'
#' @references
#' Fan, J., Guo, j. and Zheng, S. (2020)
#' Estimating number of factors by adjusted eigenvalues thresholding
#' Journal of the American Statistical Association, 117(538): 852--861
#'
#' @keywords internal
select.rank.act = function (Y, maxcomp = NULL) {
  # Set the data dimensions
  n = nrow(Y); p = ncol(Y); d = maxcomp
  d = ifelse(is.null(maxcomp), p, d)

  # Compute the spectrum of the correlation matrix of Y
  if (p == d) {
    lambdas = eigen(cor(Y))$values
  } else {
    lambdas = RSpectra::eigs_sym(cor(Y), d)$values
    lambda0 = p - sum(lambdas)
    lambdas = c(lambdas, rep(lambda0, p - d))
  }

  # ACT selection loop
  m = rep(0, times = d-1)
  for (j in 1:(d-1)) {
    delta1 = 1 / (lambdas[(j+1):d] - lambdas[j])
    delta2 = 1 / ((3 * lambdas[j] + lambdas[j+1]) / 4 - lambdas[j])
    m[j] = 1 / (d-j) * (sum(delta1) + delta2)
  }

  # Rank selection
  rho = (d - 1:(d-1)) / (n-1)
  m1 = rho * m - (1 - rho) / lambdas[1:(d-1)]
  adj.lambdas = - 1 / m1
  ncomp = sum(adj.lambdas > 1 + sqrt(d / n))
  ncomp = max(1, ncomp)

  # Return the selected rank
  return (ncomp)
}

