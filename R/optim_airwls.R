
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
#' @param intercept ...
#' @param dispersion ...
#' @param penalty ...
#' @param init ...
#' @param control ...
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
gmf.airwls = function (
    Y,
    X = NULL,
    Z = NULL,
    family = poisson(),
    ncomp = 2,
    intercept = c(FALSE, FALSE),
    dispersion = FALSE,
    penalty = list(),
    init = list(),
    control = list()) {

  # Save the initial time
  init.time = proc.time()

  # Set control, penalty and initialization parameters
  control = set.control("airwls", control)
  penalty = set.penalty(penalty)
  init = set.init(init)

  maxiter = control$maxiter
  stepsize = control$stepsize
  nsteps = control$nsteps
  tol = control$tol
  damping = control$damping
  verbose = control$verbose
  frequency = control$frequency
  parallel = 1


  # Derive dimensions
  n = nrow(Y)
  m = ncol(Y)
  nm = n * m

  # Range of the data
  ymin = min(Y, na.rm = TRUE)
  ymax = max(Y, na.rm = TRUE)
  yrng = ymax - ymin
  yup = ymax + control$eps * yrng
  ylo = ymin - control$eps * yrng

  # Add the intercept as a column of ones in X or Z
  if (intercept[1]) X = cbind(matrix(1, nrow = n, ncol = 1), X)
  if (intercept[2]) Z = rbind(matrix(1, nrow = m, ncol = 1), Z)

  # X still may be NULL if no intercept and empty input X. Then p = 0
  p = 0
  q = 0
  d = ncomp
  if (!is.null(X)) p = ncol(X)
  if (!is.null(Z)) q = ncol(Z)

  # Initialize U, V and beta using the selected method
  .init = gmf.init(Y, X, Z, d, family, init$method,
                   init$niter, init$values, init$verbose)

  # Build the left and right decomposition matrices
  U = cbind(X, .init$bz, .init$u)
  V = cbind(.init$bx, Z, .init$v)

  # Build the penalization vectors for the U and V columns
  pu = c(rep(0, p), rep(penalty$b, q), rep(penalty$u, d))
  pv = c(rep(penalty$b, p), rep(0, q), rep(penalty$v, d))

  # Get the column indices of (A,U) and (B,V)
  idu = c((p+1):(p+q), (p+q+1):(p+q+d))
  idv = c(1:p, (p+q+1):(p+q+d))

  # Save the optimization history
  trace = list(
    iter = numeric(maxiter),
    change = numeric(maxiter),
    deviance = numeric(maxiter),
    exe.time = numeric(maxiter))

  # Remember the last deviance
  llast = Inf

  # remember where is missing data so that we can keep replacing it
  # with more and more accurate models
  isna = is.na(Y)
  anyna = any(isna)

  # Get the linear predictor
  eta = tcrossprod(U, V)

  # Get the predicted means
  mu = family$linkinv(eta)

  # Fill the missing values using the initial values
  if (anyna) Y[isna] = mu[isna]

  # Check whether there are extreme predictions
  above = which(mu > yup)
  below = which(mu < ylo)
  check = length(above) > 0 | length(below) > 0
  if (check) {
    # Correct extreme predictions
    mu[above] = yup
    mu[below] = ylo
    # Becktransform the linear predictor
    eta[above] = family$linkfun(yup)
    eta[below] = family$linkfun(ylo)
  }

  # Get the starting penalty and deviance
  penaltyu = matrix.penalty(U[,idu], pu[idu])
  penaltyv = matrix.penalty(V[,idv], pv[idv])
  deviance = matrix.deviance(mu, Y, family)
  objective = deviance + 0.5 * (penaltyu + penaltyv)
  change = Inf

  # Save the current time
  mid.time = proc.time()
  exe.time = as.numeric(mid.time - init.time)[3]

  # Standard output
  if (verbose) {
    cat(rep("-", 44), "\n", sep = "")
    print.status(1, deviance / nm, 1, exe.time)
  }

  for (t in 1:maxiter) {
    # After the first iteration, replace NAs with model values
    if (t > 1) Y[isna] = mu[isna]

    # Perform airwls_internal_steps internal steps of the regularized IRWLS
    for (j in 1:nsteps) {
      # Get column coefficients, correct for regression using the offset
      offset = NULL
      if (q) offset = tcrossprod(U[, (p+1):(p+q), drop = FALSE], Z)
      coefs = slice.glm(U[,idv], Y, slices = m, coefs = t(V[,idv]),
                        offset = offset, penalized = penalty$v, parallel = 1,
                        family = family, method = "step", stepsize = stepsize)
      V[,idv] = t(coefs)
    }
    for (j in 1:nsteps) {
      # Get row coefficients, correct for regression using the offset
      offset = NULL
      if (p) offset = t(tcrossprod(X, V[, 1:p, drop = FALSE]))
      coefs = t(slice.glm(V[,idu], t(Y), slices = n, coefs = t(U[,idu]),
                          offset = offset, penalized = penalty$u, parallel = 1,
                          family = family, method = "step", stepsize = stepsize))
      U[,idu] = coefs
    }

    # Get the linear predictor
    eta = tcrossprod(U, V)

    # Get the predicted means
    mu = family$linkinv(eta)

    # Correct extremely wrong values
    above = which(mu > yup)
    below = which(mu < ylo)
    mu[above] = yup
    mu[below] = ylo

    # Becktransform the linear predictor
    eta[above] = family$linkfun(mu[above])
    eta[below] = family$linkfun(mu[below])

    # Get the predicted variances
    var = family$variance(mu)

    # Get the dispersion parameters
    phi = if (dispersion) get.dispersion(Y, mu, family) else rep(1, m)

    # Compute penalties and mean deviance
    penaltyu = matrix.penalty(U[,idu], pu[idu])
    penaltyv = matrix.penalty(V[,idv], pv[idv])
    deviance = matrix.deviance(mu, Y, family)
    objectivet = objective
    objective = deviance + 0.5 * (penaltyu + penaltyv)
    change = abs(objective - objectivet) / (abs(objectivet) + 1e-04)

    # Save the current time
    mid.time = proc.time()

    # Save the partial execution time
    exe.time = as.numeric(mid.time - init.time)[3]

    # Diagnostic output: iteration, deviance, relative change
    if (t %% frequency == 0 && verbose) {
      print.status(t, deviance / nm, change, exe.time)
    }

    # Store the current optimization state
    trace$iter[t] = t
    trace$change[t] = change
    trace$deviance[t] = deviance
    trace$exe.time[t] = exe.time

    # Check the stopping criteria
    if (change < tol){
      break
    } else{
      llast = deviance
    }
  }

  # Prune the unused records in the trace
  if (t < maxiter) {
    trace$iter = trace$iter[1:t]
    trace$change = trace$change[1:t]
    trace$deviance = trace$deviance[1:t]
    trace$exe.time = trace$exe.time[1:t]
  }

  # Get the final dispersion parameter
  phi = if (dispersion) get.dispersion(Y, mu, family) else rep(1, m)

  # Get the estimated coefficients
  Bx = Bz = NULL
  if (q) Bz = U[, (p+1):(p+q), drop = FALSE]
  if (p) Bx = V[, 1:p, drop = FALSE]
  U = U[, (p+q+1):(p+q+d), drop = FALSE]
  V = V[, (p+q+1):(p+q+d), drop = FALSE]

  # Orthogonalize U and V
  if (control$normalize) {
    s <- svd::propack.svd(tcrossprod(U, V), neig = d)
    U <- s$u %*% diag(sqrt(s$d))
    V <- s$v %*% diag(sqrt(s$d))
    D <- s$d
  }

  # Save the ending time
  end.time = proc.time()

  # Final execution time
  exe.time = as.numeric(end.time - init.time)[3]

  if (verbose) {
    cat(rep("-", 44), "\n", sep = "")
    print.status(t, deviance / nm, change, exe.time)
  }

  # Return the optimization result
  list(
    data = list(Y = Y, X = X, Z = Z, isna = isna),
    dims = list(n = n, m = m, d = d, p = p, q = q),
    coef = list(betaX = Bx, betaZ = Bz, U = U, V = V, delta = D),
    pred = list(eta = eta, mu = mu, var = var),
    dispersion = phi,
    deviance = deviance,
    exe.time = exe.time,
    method = "airwls",
    penalty = penalty,
    control = control,
    init = init,
    trace = trace)
}

## list(v = V,
##      u = U,
##      d = D,
##      beta.x = Bx,
##      beta.z = Bz,
##      phi = phi,
##      family = family,
##      linear.predictor = eta,
##      fitted.values = mu,
##      fitted.variances = var,
##      deviance = deviance,
##      penalty = penaltyu + penaltyv,
##      execution.time = exe.time,
##      trace = trace)
