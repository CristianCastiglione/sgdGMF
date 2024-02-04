
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

  # Set and check the estimation method
  methods = c("airwls", "newton", "csgd", "rsgd", "bsgd")

  if (!(method %in% methods)) {
    stop("Not allowed estimation method.")
  }

  # Set and check the control parameters
  control = set.control(method, control)
  lambda = set.penalty(penalty)
  init = set.init(init)
  familyname = family$family
  linkname = family$link

  # Select the correct estimation method
  if (method == "airwls") {
    fit = cpp.fit.airwls(
      Y = Y, X = init$X, B = init$B, A = init$A, Z = init$Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = control$maxiter, nstep = control$nstep, stepsize = control$stepsize,
      eps = control$eps, nafill = control$nafill, tol = control$tol,
      damping = control$damping, verbose = control$verbose,
      frequency = control$frequency, parallel = control$parallel,
      nthreads = control$nthreads
    )
  }
  if (method == "newton") {
    fit = cpp.fit.newton(
      Y = Y, X = init$X, B = init$B, A = init$A, Z = init$Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = control$maxiter, stepsize = control$stepsize, eps = control$eps,
      nafill = control$nafill, tol = control$tol, damping = control$damping,
      verbose = control$verbose, frequency = control$frequency,
      parallel = control$parallel, nthreads = control$nthreads
    )
  }
  if (method == "msgd") {
    fit = cpp.fit.msgd(
      Y = Y, X = init$X, B = init$B, A = init$A, Z = init$Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = control$maxiter, eps = control$eps, nafill = control$nafill,
      tol = control$tol, size = control$size,
      burn = control$burn, rate0 = control$rate0, decay = control$decay,
      damping = control$damping, rate1 = control$rate[1], rate2 = control$rate[2],
      parallel = control$parallel, nthreads = control$nthreads,
      verbose = control$verbose, frequency = control$frequency,
      progress = control$progress
    )
  }
  if (method == "csgd") {
    fit = cpp.fit.csgd(
      Y = Y, X = init$X, B = init$B, A = init$A, Z = init$Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = control$maxiter, eps = control$eps, nafill = control$nafill,
      tol = control$tol, size1 = control$size[1], size2 = control$size[2],
      burn = control$burn, rate0 = control$rate0, decay = control$decay,
      damping = control$damping, rate1 = control$rate[1], rate2 = control$rate[2],
      parallel = control$parallel, nthreads = control$nthreads,
      verbose = control$verbose, frequency = control$frequency,
      progress = control$progress
    )
  }
  if (method == "rsgd") {
    fit = cpp.fit.rsgd(
      Y = Y, X = init$X, B = init$B, A = init$A, Z = init$Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = control$maxiter, eps = control$eps, nafill = control$nafill,
      tol = control$tol, size1 = control$size[1], size2 = control$size[2],
      burn = control$burn, rate0 = control$rate0, decay = control$decay,
      damping = control$damping, rate1 = control$rate[1], rate2 = control$rate[2],
      parallel = control$parallel, nthreads = control$nthreads,
      verbose = control$verbose, frequency = control$frequency,
      progress = control$progress
    )
  }
  if (method == "bsgd") {
    fit = cpp.fit.bsgd(
      Y = Y, X = init$X, B = init$B, A = init$A, Z = init$Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = control$maxiter, eps = control$eps, nafill = control$nafill,
      tol = control$tol, size1 = control$size[1], size2 = control$size[2],
      burn = control$burn, rate0 = control$rate0, decay = control$decay,
      damping = control$damping, rate1 = control$rate[1], rate2 = control$rate[2],
      parallel = control$parallel, nthreads = control$nthreads,
      verbose = control$verbose, frequency = control$frequency,
      progress = control$progress
    )
  }

  return (fit)
}
