
#' @title Factorize a matrix of non-Gaussian observations via SGD
#'
#' @description
#' Fit a generalized matrix factorization (GMF) model for non-Gaussian data using
#' either deterministic or stochastic optimization methods.
#' It is an alternative to PCA when the observed data are binary, counts, positive
#' scores or, more generally, when the conditional distribution of the observations
#' can be appropriately described using a dispersion exponential family
#' or a quasi-likelihood model.
#' Some examples are Gaussian, Gamma, Binomial and Poisson laws.
#'
#' The dependence among the observations and the variables in the sample can be
#' taken into account through appropriate row- and column-specific regression effects.
#' The residual variability is then modeled through a low-rank latent matrix factorization.
#'
#' For the estimation, the package implements two deterministic optimization methods
#' and four stochastic optimization algorithms.
#'
#' @param Y matrix of responses (\eqn{n \times m})
#' @param X matrix of row fixed effects (\eqn{n \times p})
#' @param Z matrix of column fixed effects (\eqn{q \times m})
#' @param family a family as in the \code{\link{glm}} interface (default \code{gaussian()})
#' @param ncomp rank of the latent matrix factorization (default 2)
#' @param method optimization method: \code{"airwls"} (default),
#' \code{"newton"}, \code{"msgd"}, \code{"csgd"}, \code{"rsgd"}, \code{"bsgd"}
#' @param penalty list of penalty parameters
#' @param init list of initialization options
#' @param control list of optimization options
#'
#' @return
#' An \code{sgdgmf} object, namely a list, containing the estimated parameters of the GMF model.
#' In particular, the returned object collects the following information:
#' \itemize{
#'   \item \code{method}: the selected estimation method
#'   \item \code{family}: the model family name
#'   \item \code{link}: the link function name
#'   \item \code{idu}: a vector of indices identify the columns of the augmented
#'     matrix \eqn{U_* = [X, \Gamma, U]} corresponding to \eqn{\Gamma} and \eqn{U}
#'   \item \code{idv}: a vector of indices identify the columns of the augmented
#'     matrix \eqn{V_* = [B, Z, V]} corresponding to \eqn{B} and \eqn{V}
#'   \item \code{U}: the estimated factor matrix \eqn{\hat{U}_* = [X, \hat{\Gamma}, \hat{U}]}
#'   \item \code{V}: the estimated loading matrix \eqn{\hat{V}_* = [\hat{B}, Z, \hat{V}]}
#'   \item \code{eta}: the estimated linear predictor
#'   \item \code{mu}: the estimated mean matrix
#'   \item \code{var}: the estimated variance matrix
#'   \item \code{phi}: the estimated dispersion parameter
#'   \item \code{penalty}: the penalty value at the end of the optimization
#'   \item \code{deviance}: the deviance value at the end of the optimization
#'   \item \code{objective}: the penalized objective function at the end of the optimization
#'   \item \code{exe.time}: the total execution time in seconds
#'   \item \code{trace}: a trace matrix recording the optimization history
#' }
#'
#'
#' @details
#' The model we consider is defined as follows.
#' Let \eqn{Y = (y_{ij})} be a matrix of observed data of dimension \eqn{n \times m}.
#' We assume for the \eqn{(i,j)}th observation in the matrix a dispersion exponential family law
#' \eqn{(y_{ij} \mid \theta_{ij}) \sim EF(\theta_{ij}, \phi)}, where \eqn{\theta_{ij}} is the
#' natural parameter and \eqn{\phi} is the dispersion parameter.
#' Recall that the conditional probability density function of \eqn{y_{ij}} is given by
#' \deqn{f (y_{ij}; \psi) = \exp \big[ \{(y_{ij} \theta_{ij} - b(\theta_{ij})\} / \phi - c(y_{ij}, \phi) \big],}
#' where \eqn{\psi} is the vector of unknown parameters to be estimated,
#' \eqn{b(\cdot)} is a convex twice differentiable log-partition function,
#' and \eqn{c(\cdot,\cdot)} is the cumulant function of the family.
#'
#' The conditional mean of \eqn{y_{ij}}, say \eqn{\mu_{ij}}, is then modeled as
#' \deqn{g(\mu_{ij}) = \eta_{ij} = x_i^\top \beta_j + \gamma_i^\top z_j + u_i^\top v_j,}
#' where \eqn{g(\cdot)} is a bijective twice differentiable link function, \eqn{\eta_{ij}} is
#' a linear predictor, \eqn{x_i \in \mathbb{R}^p} and \eqn{z_j \in \mathbb{R}^q} are
#' observed covariate vectors, \eqn{\beta_j \in \mathbb{R}^p} and \eqn{\gamma_j \in \mathbb{R}^q}
#' are unknown regression parameters and, finally, \eqn{u_i \in \mathbb{R}^d} and
#' \eqn{v_j \in \mathbb{R}^d} are latent vector explaining the residual varibility
#' not captured by the regression effects.
#' Equivalently, in matrix form, we have
#' \eqn{g(\mu) = \eta = X B^\top + \Gamma Z^\top + U V^\top.}
#'
#' The natural parameter \eqn{\theta_{ij}} is linked to the conditional mean of \eqn{y_{ij}}
#' through the equation \eqn{E(y_{ij}) = \mu_{ij} = b'(\theta_{ij})}.
#' Similarly, the variance of \eqn{y_{ij}} is given by
#' \eqn{\text{Var}(y_{ij}) = \phi \,\nu(\mu_{ij}) = \phi \,b''(\mu_{ij})}, where \eqn{\nu(\cdot)}
#' is the so-called variance function of the family.
#'
#' The estimation of the model parameters is performed by minimizing the penalized negative log-likelihood function
#' \deqn{\ell_\lambda (\psi; y) = - 2 \sum_{i,j = 1}^{n,m} \log f(y_{ij}; \psi) + \lambda \| U \|_F^2 + \lambda \| V \|_F^2,}
#' where \eqn{\lambda > 0} is a regularization parameter and \eqn{\|\cdot\|} is the Frobenious norm.
#' Quasi-likelihood models, where \eqn{f(y; \psi)} is substituted by
#' \eqn{Q(y, \mu) = \exp \big[ \int_y^\mu \{(y - t) / \phi \nu(t)\} \,dt \big] / \phi},
#' can also be considered.
#' The regularization parameters for \eqn{U} and \eqn{V} do not have necessary the same value.
#' Additional \eqn{\ell_2} penalization terms can be introduced to regularize \eqn{B} and \eqn{\Gamma}.
#'
#' To obtain the penalized maximum likelihood estimate, we here employs
#' six different algorithms
#' \itemize{
#'   \item AIRWLS: alternated iterative re-weighted least squares
#'   \item Newton: quasi-Newton algorithm with diagonal Hessian
#'   \item M-SGD: memoized stochastic gradient descent
#'   \item C-SGD: coordinate-wise stochastic gradient descent
#'   \item R-SGD: row-wise stochastic gradient descent
#'   \item B-SGD: block-wise stochastic gradient descent
#' }
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
    method = c("airwls", "newton", "msgd", "csgd", "rsgd", "bsgd"),
    penalty = list(),
    init = list(),
    control = list()
) {

  # Check and set the model matrices
  Y = set.mat.Y(Y)
  X = set.mat.X(X, nrow(Y), "X")
  Z = set.mat.X(Z, ncol(Y), "Z")

  # Check and set the control parameters
  method = match.arg(method)
  ctr = set.control(method, control)
  lambda = set.penalty(penalty)
  init = set.init(init)
  familyname = family$family
  linkname = family$link

  # Initialize the parameters
  init = init.param(
    Y = Y, X = X, Z = Z, ncomp = ncomp, family = family,
    method = init$method, niter = init$niter,
    values = init$values, verbose = init$verbose)

  # Select the correct estimation method
  if (method == "airwls") {
    # AIRWLS algorithm
    fit = cpp.fit.airwls(
      Y = Y, X = X, B = init$B, A = init$A, Z = Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = ctr$maxiter, nstep = ctr$nstep, stepsize = ctr$stepsize,
      eps = ctr$eps, nafill = ctr$nafill, tol = ctr$tol,
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
      damping = ctr$damping, rate1 = control$rate1, rate2 = ctr$rate2,
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
      decay = ctr$decay, damping = ctr$damping, rate1 = ctr$rate1, rate2 = ctr$rate2,
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
      decay = ctr$decay, damping = ctr$damping, rate1 = ctr$rate1, rate2 = ctr$rate2,
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
      decay = ctr$decay, damping = ctr$damping, rate1 = ctr$rate1, rate2 = ctr$rate2,
      parallel = ctr$parallel, nthreads = ctr$nthreads, verbose = ctr$verbose,
      frequency = ctr$frequency, progress = ctr$progress
    )
  }

  # Model dimensions
  n = nrow(Y); m = ncol(Y)
  p = ncol(X); q = ncol(Z)

  df = p * m + q * n + (n + m) * ncomp
  nm = n * m - sum(is.na(Y))

  # Set the indices of A, B, U and V
  idxA = seq(from = p+1, to = p+q)
  idxB = seq(from = 1, to = p)
  idxU = seq(from = p+q+1, to = p+q+ncomp)
  idxV = seq(from = p+q+1, to = p+q+ncomp)

  # Output list
  out = list()
  out$method = fit$method
  out$family = family
  out$ncomp = ncomp
  out$control = ctr
  out$A = fit$U[, idxA]
  out$B = fit$V[, idxB]
  out$U = fit$U[, idxU]
  out$V = fit$V[, idxV]
  out$eta = fit$eta
  out$mu = fit$mu
  out$var = fit$var
  out$phi = fit$phi
  out$penalty = fit$penalty
  out$deviance = fit$deviance
  out$objective = fit$objective
  out$aic = fit$deviance + 2 * df
  out$bic = fit$deviance + 2 * df * log(nm)
  out$exe.time = fit$exe.time
  out$trace = as.data.frame(fit$trace)
  colnames(out$trace) = c("iter", "dev", "pen", "pdev", "change", "time")

  # Normalize the latent factors
  if (ctr$normalize) {
    uv = normalize.uv(out$U, out$V, method = "qr")
    out$U = uv$U
    out$V = uv$V
  }

  return (out)
}
