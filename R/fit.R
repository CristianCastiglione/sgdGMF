
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
sgdgmf.fit = function (
    Y,
    X = NULL,
    Z = NULL,
    family = poisson(),
    ncomp = 2,
    method = "airwls",
    penalty = list(),
    init = list(),
    control = list()) {

  # Check and set the covariate matrices
  Y = set.mat.Y(Y)
  X = set.mat.X(X, nrow(Y), "X")
  Z = set.mat.X(Z, ncol(Y), "Z")

  # Check and set the control parameters
  ctr = set.control(method, control)
  lambda = set.penalty(penalty)
  init = set.init(init)
  familyname = family$family
  linkname = family$link

  # Select the correct estimation method
  if (method == "airwls") {
    # AIRWLS algorithm
    fit = cpp.fit.airwls(
      Y = Y, X = X, B = init$B, A = init$A, Z = Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = ctr$maxiter, nstep = ctr$nstep, stepsize = ctr$stepsize,
      eps = ctr$eps, nafill = ctr$nafill, tol = control$tol,
      damping = ctr$damping, verbose = ctr$verbose,
      frequency = ctr$frequency, parallel = ctr$parallel,
      nthreads = ctr$nthreads
    )
  }
  if (method == "newton") {
    # Quasi-Newton algorithm
    fit = cpp.fit.newton(
      Y = Y, X = X, B = init$B, A = init$A, Z = Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = ctr$maxiter, stepsize = ctr$stepsize, eps = ctr$eps,
      nafill = ctr$nafill, tol = ctr$tol, damping = ctr$damping,
      verbose = ctr$verbose, frequency = ctr$frequency,
      parallel = ctr$parallel, nthreads = ctr$nthreads
    )
  }
  if (method == "msgd") {
    # Memoized SGD algorithm
    fit = cpp.fit.msgd(
      Y = Y, X = X, B = init$B, A = init$A, Z = Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = ctr$maxiter, eps = ctr$eps, nafill = ctr$nafill, tol = ctr$tol,
      size = ctr$size, burn = ctr$burn, rate0 = ctr$rate0, decay = ctr$decay,
      damping = ctr$damping, rate1 = control$rate[1], rate2 = ctr$rate[2],
      parallel = ctr$parallel, nthreads = ctr$nthreads, verbose = ctr$verbose,
      frequency = ctr$frequency, progress = ctr$progress
    )
  }
  if (method == "csgd") {
    # Coordinatewise SGD algorithm
    fit = cpp.fit.csgd(
      Y = Y, X = X, B = init$B, A = init$A, Z = Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = ctr$maxiter, eps = ctr$eps, nafill = ctr$nafill, tol = ctr$tol,
      size1 = ctr$size[1], size2 = ctr$size[2], burn = ctr$burn, rate0 = ctr$rate0,
      decay = ctr$decay, damping = ctr$damping, rate1 = ctr$rate[1], rate2 = ctr$rate[2],
      parallel = ctr$parallel, nthreads = ctr$nthreads, verbose = ctr$verbose,
      frequency = ctr$frequency, progress = ctr$progress
    )
  }
  if (method == "rsgd") {
    # Rowwise SGD algorithm
    fit = cpp.fit.rsgd(
      Y = Y, X = X, B = init$B, A = init$A, Z = Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = ctr$maxiter, eps = ctr$eps, nafill = ctr$nafill, tol = ctr$tol,
      size1 = ctr$size[1], size2 = ctr$size[2], burn = ctr$burn, rate0 = ctr$rate0,
      decay = ctr$decay, damping = ctr$damping, rate1 = ctr$rate[1], rate2 = ctr$rate[2],
      parallel = ctr$parallel, nthreads = ctr$nthreads, verbose = ctr$verbose,
      frequency = ctr$frequency, progress = ctr$progress
    )
  }
  if (method == "bsgd") {
    # Blockwise SGD algorithm
    fit = cpp.fit.bsgd(
      Y = Y, X = X, B = init$B, A = init$A, Z = Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = ctr$maxiter, eps = ctr$eps, nafill = ctr$nafill, tol = ctr$tol,
      size1 = ctr$size[1], size2 = ctr$size[2], burn = ctr$burn, rate0 = ctr$rate0,
      decay = ctr$decay, damping = ctr$damping, rate1 = ctr$rate[1], rate2 = ctr$rate[2],
      parallel = ctr$parallel, nthreads = ctr$nthreads, verbose = ctr$verbose,
      frequency = ctr$frequency, progress = ctr$progress
    )
  }

  return (fit)
}
