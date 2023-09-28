
#' @title Factorize a matrix of non-Gaussian observations using averaged SGD
#'
#' @description
#' Fit a penalized generalized matrix factorization model using a coordinate
#' averaged stochastic gradient descent
#'
#' @param Y ...
#' @param X ...
#' @param Z ...
#' @param family ...
#' @param ncomp ...
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
gmf.coord.sgd = function (
    Y,
    X = NULL,
    Z = NULL,
    family = gaussian(),
    ncomp = 2,
    intercept = c(FALSE, FALSE),
    dispersion = FALSE,
    penalty = list(),
    init = list(),
    control = list()) {

  # Save the initial time
  init.time = proc.time()

  # Set control, penalty and initialization parameters
  control = set.control("csgd", control)
  penalty = set.penalty(penalty)
  init = set.init(init)

  # normalize = control$normalize
  nafill = control$nafill
  maxiter = control$maxiter
  # size = control$size
  # eps = control$eps
  rate0 = control$rate0
  decay = control$decay
  damping = control$damping
  rate_g = control$rate1
  rate_h = control$rate2
  burn = control$burn
  verbose = control$verbose
  frequency = control$frequency

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
  if (!is.null(X) & is.matrix(X)) p = ncol(X)
  if (!is.null(Z) & is.matrix(Z)) q = ncol(Z)

  # Normalization factors for minibatch estimate corrections
  sizer = min(n, control$size[1])
  sizec = min(m, control$size[2])

  scaler = n / sizer
  scalec = m / sizec

  # Column and row minibatch samples
  mbr = sample.minibatch(n, sizer, randomize = TRUE)
  mbc = sample.minibatch(m, sizec, randomize = TRUE)

  # Initialize U, V and beta using the selected method
  .init = gmf.init(Y, X, Z, d, family, init$method,
                   init$niter, init$values, init$verbose)

  # Build the left and right decomposition matrices
  U = Ut = cbind(X, .init$bz, .init$u)
  V = Vt = cbind(.init$bx, Z, .init$v)

  # Build the penalization vectors for the U and V columns
  pu = c(rep(0, p), rep(penalty$b, q), rep(penalty$u, d))
  pv = c(rep(penalty$b, p), rep(0, q), rep(penalty$v, d))

  # Get the column indices of (A,U) and (B,V)
  idu = c((p+1):(p+q), (p+q+1):(p+q+d))
  idv = c(1:p, (p+q+1):(p+q+d))

  # Initialize the gradients and Hessians of U and V
  du = matrix(0, nrow = n, ncol = q+d)
  dv = matrix(0, nrow = m, ncol = p+d)
  ddu = matrix(0, nrow = n, ncol = q+d)
  ddv = matrix(0, nrow = m, ncol = p+d)

  # Save the optimization history
  trace = list(
    iter = numeric(),
    change = numeric(),
    deviance = numeric(),
    exe.time = numeric())

  # remember where there are missing data so that we can keep
  # replacing them with more and more accurate models
  isna = is.na(Y)
  anyna = any(isna)

  # Get the linear predictor
  eta = tcrossprod(U, V)

  # Get the predicted means
  mu = family$linkinv(eta)

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

  # Fill the missing values using the initial values
  if (anyna) Y[isna] = mu[isna]

  # Get the starting penalty and deviance
  penaltyu = matrix.penalty(U[,idu], pu[idu])
  penaltyv = matrix.penalty(V[,idv], pv[idv])
  deviance = matrix.deviance(mu, Y, family)
  objective = deviance + 0.5 * penaltyu + 0.5 * penaltyv
  change = Inf

  # Save the current time
  mid.time = proc.time()
  exe.time = as.numeric(mid.time - init.time)[3]

  # Store the current optimization state
  trace$iter = c(trace$iter, 0)
  trace$change = c(trace$change, Inf)
  trace$deviance = c(trace$deviance, deviance)
  trace$exe.time = c(trace$exe.time, exe.time)

  # Print the optimization status
  if (verbose) {
    cat(rep("-", 72), "\n", sep = "")
    print.status(1, deviance / nm, 1, exe.time, 0, 0)
  }

  # Optimization loop
  for (t in 1:maxiter) {

    # After the first iteration, replace NAs with model values
    if (t > 1){
      # To save time we do not fill the NA cells at every iteration
      if (t %% nafill == 0) {

        # Check whether there are extreme predictions
        above = which(mu > yup)
        below = which(mu < ylo)
        check = length(above) > 0 | length(below) > 0

        if (check) {
          # Correct extreme predictions
          mu[above] = yup
          mu[below] = ylo

          # Beck-transform the linear predictor
          eta[above] = family$linkfun(mu[above])
          eta[below] = family$linkfun(mu[below])
        }

        # Impute the missing values
        if (anyna) Y[isna] = mu[isna]
      }
    }

    # Update the learning rate
    ratet = update.sgd.rate(t, rate0, decay)

    # Sample the mini-batch indices
    idxr = mbr[[select.minibatch(t, length(mbr))]]
    idxc = mbc[[select.minibatch(t, length(mbc))]]

    # Normalization factors for minibatch estimate corrections
    scaler = n / length(idxr)
    scalec = m / length(idxc)

    # Linear predictor
    etar = tcrossprod(U[idxr, , drop = FALSE], V)
    etac = tcrossprod(U, V[idxc, , drop = FALSE])

    yr = Y[idxr, , drop = FALSE]
    yc = Y[, idxc, drop = FALSE]

    eta[idxr, ] = etar
    eta[, idxc] = etac

    mu[idxr, ] = family$linkinv(etar)
    mu[, idxc] = family$linkinv(etac)

    # Row and column minibatch derivatives of the likelihood wrt eta
    derivr = get.eta.deriv(yr, etar, family, trasp = TRUE)
    derivc = get.eta.deriv(yc, etac, family, trasp = FALSE)

    # U and V updates via averaged stochastic gradient
    dv = update.sgd.grad(V[, idv, drop = FALSE], U[idxr, idv, drop = FALSE], derivr$deta, dv, rate_g, scaler, penalty$v)
    du = update.sgd.grad(U[, idu, drop = FALSE], V[idxc, idu, drop = FALSE], derivc$deta, du, rate_g, scalec, penalty$u)
    ddv = update.sgd.hess(V[, idv, drop = FALSE], U[idxr, idv, drop = FALSE], derivr$ddeta, ddv, rate_h, scaler, penalty$v, damping)
    ddu = update.sgd.hess(U[, idu, drop = FALSE], V[idxc, idu, drop = FALSE], derivc$ddeta, ddu, rate_h, scalec, penalty$u, damping)

    U[, idu] = update.sgd.params(U[, idu, drop = FALSE], du, ddu, ratet)
    V[, idv] = update.sgd.params(V[, idv, drop = FALSE], dv, ddv, ratet)
    Ut[, idu] = smooth.sgd.params(U[, idu, drop = FALSE], Ut[, idu, drop = FALSE], burn, t)
    Vt[, idv] = smooth.sgd.params(V[, idv, drop = FALSE], Vt[, idv, drop = FALSE], burn, t)

    # Diagnostic output: iteration, deviance, and execution time
    if (t %% frequency == 0) {

      # Compute the deviance and the penalization terms
      penaltyu = matrix.penalty(U[,idu], pu[idu])
      penaltyv = matrix.penalty(V[,idv], pv[idv])
      deviance = matrix.deviance(mu, Y, family)
      objectivet = objective
      objective = deviance + 0.5 * penaltyu + 0.5 * penaltyv
      change = abs(objective - objectivet) / (abs(objectivet) + 1e-04)

      # Save the current time
      mid.time = proc.time()

      # Save the partial execution time
      exe.time = as.numeric(mid.time - init.time)[3]

      # Print the optimization status
      rowfrac = floor(100 * t * sizer / n)
      colfrac = floor(100 * t * sizec / m)
      if (verbose) print.status(t, deviance / nm, change, exe.time, rowfrac, colfrac)

      # Store the current optimization state
      trace$iter = c(trace$iter, t)
      trace$change = c(trace$change, change)
      trace$deviance = c(trace$deviance, deviance)
      trace$exe.time = c(trace$exe.time, exe.time)
    }
  }

  # Return the smoothed estimates
  Bx = Bz = NULL
  if (p) Bx = Vt[, 1:p, drop = FALSE]
  if (q) Bz = Ut[, (p+1):(p+q), drop = FALSE]
  U = Ut[, (p+q+1):(p+q+d), drop = FALSE]
  V = Vt[, (p+q+1):(p+q+d), drop = FALSE]

  # Final eta, mean and variance values
  eta = tcrossprod(Ut, Vt)
  mu  = matrix(family$linkinv(eta), nrow = n, ncol = m)
  var = matrix(family$variance(mu), nrow = n, ncol = m)

  # Check whether there are extreme predictions
  above = which(mu > yup)
  below = which(mu < ylo)
  check = length(above) > 0 | length(below) > 0
  if (check) {
    # Correct extreme predictions
    mu[above] = yup
    mu[below] = ylo
    # Correct extreme variances
    var[above] = family$variance(yup)
    var[below] = family$variance(ylo)
    # Becktransform the linear predictor
    eta[above] = family$linkfun(yup)
    eta[below] = family$linkfun(ylo)
  }

  # Final penalty and deviance
  penaltyu = matrix.penalty(Ut[, idu], pu[idu])
  penaltyv = matrix.penalty(Vt[, idv], pv[idv])
  deviance = matrix.deviance(mu, Y, family)

  # Final dispersion parameter
  phi = if (dispersion) get.dispersion(Y, mu, family) else rep(1, m)

  # Normalize U and V so that they are the left and right singualar vectors of U*Vt
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

  # Print the optimization status
  if (verbose) {
    cat(rep("-", 72), "\n", sep = "")
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
    method = "c-sgd",
    penalty = penalty,
    control = control,
    init = init,
    trace = trace)
}

## list(u = U,
##      v = V,
##      d = D,
##      x = X,
##      z = Z,
##      beta.x = Bx,
##      beta.z = Bz,
##      phi = phi,
##      n = n,
##      m = m,
##      ncomp = d,
##      family = family,
##      linear.predictor = eta,
##      fitted.values = mu,
##      fitted.variances = var,
##      deviance = deviance,
##      penalty = penaltyu + penaltyv,
##      execution.time = exe.time,
##      trace = trace)
