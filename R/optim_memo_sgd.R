
#' @title Factorize a matrix of non-Gaussian observations using memoized SGD
#'
#' @description
#' Fit a penalized generalized matrix factorization model using a coordinate
#' averaged stochastic gradient descent with memoized updates
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
gmf.memo.sgd = function (
    Y,
    X = NULL,
    Z = NULL,
    family = gaussian(),
    ncomp = 2,
    intercept = FALSE,
    dispersion = FALSE,
    penalty = list(),
    init = list(),
    control = list()) {

  # Save the initial time
  init.time = proc.time()

  # Set control, penalty and initialization parameters
  control = set.control("msgd", control)
  penalty = set.penalty(penalty)
  init = set.init(init)

  # normalize = control$normalize
  nafill = control$nafill
  maxiter = control$maxiter
  epochs = control$epochs
  rate0 = control$rate0
  decay = control$decay
  damping = control$damping
  rate_g = control$rate1
  rate_h = control$rate2
  burn = control$burn
  verbose = control$verbose
  frequency = control$frequency
  progress = control$progress

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
  size = min(n, control$size)
  scale = n / size

  # Column and row minibatch samples
  mb = sample.minibatch(n, size, randomize = TRUE)
  nmb = length(mb)

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

  # Initialized the differentials of V to be accumulated along one epoch
  dV = array(0, dim = c(m, p+d, nmb))
  ddV = array(0, dim = c(m, p+d, nmb))

  # S = array(0, dim = c(ncol(Y), length(idv), length(mb)))

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
    print.status.2(1, deviance / nm, 1, exe.time, 0)
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

    # Cycle over the minibatch subsamples
    if (verbose & progress) progress.bar(0, 0, 50)
    for (b in 1:nmb) {
      # Sample the mini-batch indices
      idx = mb[[select.minibatch(b, nmb)]]

      # Normalization factors for minibatch estimate corrections
      scale = n / length(idx)

      # Linear predictor
      etat = tcrossprod(U[idx, , drop = FALSE], V)
      yt = Y[idx, , drop = FALSE]

      eta[idx, ] = etat
      mu[idx, ] = family$linkinv(etat)

      # Row and column minibatch derivatives of the likelihood wrt eta
      deriv = get.eta.deriv(yt, etat, family, trasp = TRUE)

      # U averaged stochastic gradient estimation
      du[idx, ]  = update.sgd.grad(U[idx, idu], V[, idu, drop = FALSE], t(deriv$deta), du[idx, ], rate_g, 1, penalty$u)
      ddu[idx, ] = update.sgd.hess(U[idx, idu], V[, idu, drop = FALSE], t(deriv$ddeta), ddu[idx, ], rate_h, 1, penalty$u, damping)

      # V averaged stochastic gradient estimation
      dV[, , b]  = update.sgd.grad(V[, idv], U[idx, idv, drop = FALSE], deriv$deta, dv, rate_g, scale, penalty$v)
      ddV[, , b] = update.sgd.hess(V[, idv], U[idx, idv, drop = FALSE], deriv$ddeta, ddv, rate_h, scale, penalty$v, damping)

      # U update via averaged stochastic gradient
      U[idx, idu] = update.sgd.params(U[idx, idu], du[idx, ], ddu[idx, ], ratet)
      Ut[, idu] = smooth.sgd.params(U[, idu], Ut[, idu], burn, t)

      # Update the progress bar
      if (verbose & progress) {
        progress.bar(floor(50 * b / nmb), 0, 50)
        cat("\b\r")
      }
    }

    # V gradient accumulation
    dv = rowMeans(dV, dims = 2)
    ddv = rowMeans(ddV, dims = 2)

    # V update via averaged stochastic gradient
    V[, idv] = update.sgd.params(V[, idv], dv, ddv, ratet)
    Vt[, idv] = smooth.sgd.params(V[, idv], Vt[, idv], burn, t)

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
      scanned = floor(100 * t * size / n)
      if (verbose) print.status.2(t, deviance / nm, change, exe.time, scanned)

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
    method = "m-sgd",
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
