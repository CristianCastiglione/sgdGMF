
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
#' @param family a \code{glm} family (see \code{\link{family}} for more details)
#' @param ncomp rank of the latent matrix factorization (default 2)
#' @param weights an optional matrix of weights (\eqn{n \times m})
#' @param offset an optional matrix of offset values (\eqn{n \times m}), that specify a known component to be included in the linear predictor.
#' @param method estimation method to minimize the negative penalized log-likelihood
#' @param penalty list of penalty parameters (see \code{\link{set.penalty}} for more details)
#' @param control.init list of control parameters for the initialization (see \code{\link{set.control.init}} for more details)
#' @param control.alg list of control parameters for the optimization (see \code{\link{set.control.alg}} for more details)
#'
#' @return
#' An \code{sgdgmf} object, namely a list, containing the estimated parameters of the GMF model.
#' In particular, the returned object collects the following information:
#' \itemize{
#'   \item \code{method}: the selected estimation method
#'   \item \code{family}: the model family
#'   \item \code{ncomp}:
#'   \item \code{npar}:
#'   \item \code{control.init}:
#'   \item \code{control.alg}:
#'   \item \code{control.cv}:
#'   \item \code{Y}:
#'   \item \code{X}:
#'   \item \code{Z}:
#'   \item \code{B}: the estimated col-specific coefficient matrix
#'   \item \code{A}: the estimated row-specific coefficient matrix
#'   \item \code{U}: the estimated factor matrix
#'   \item \code{V}: the estimated loading matrix
#'   \item \code{eta}: the estimated linear predictor
#'   \item \code{mu}: the estimated mean matrix
#'   \item \code{var}: the estimated variance matrix
#'   \item \code{phi}: the estimated dispersion parameter
#'   \item \code{penalty}: the penalty value at the end of the optimization
#'   \item \code{deviance}: the deviance value at the end of the optimization
#'   \item \code{objective}: the penalized objective function at the end of the optimization
#'   \item \code{aic}:
#'   \item \code{bic}:
#'   \item \code{cbic}:
#'   \item \code{exe.time}: the total execution time in seconds
#'   \item \code{trace}: a trace matrix recording the optimization history
#'   \item \code{summary.cv}:
#' }
#'
#' @details
#' The model we consider is defined as follows.
#' Let \eqn{Y = (y_{ij})} be a matrix of observed data of dimension \eqn{n \times m}.
#' We assume for the \eqn{(i,j)}th observation in the matrix a dispersion exponential family law
#' \eqn{(y_{ij} \mid \theta_{ij}) \sim EF(\theta_{ij}, \phi)}, where \eqn{\theta_{ij}} is the
#' natural parameter and \eqn{\phi} is the dispersion parameter.
#' Recall that the conditional probability density function of \eqn{y_{ij}} is given by
#' \deqn{f (y_{ij}; \psi) = \exp \big[ w_{ij} \{(y_{ij} \theta_{ij} - b(\theta_{ij})\} / \phi - c(y_{ij}, \phi) \big],}
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
#' where \eqn{\lambda > 0} is a regularization parameter and \eqn{\|\cdot\|_F} is the Frobenious norm.
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
#' Kidzinnski, L., Hui, F.K.C., Warton, D.I. and Hastie, J.H. (2022).
#' \emph{Generalized Matrix Factorization: efficient algorithms for fitting generalized linear latent variable models to large data arrays.}
#' Journal of Machine Learning Research, 23: 1-29.
#'
#' Wang, L. and Carvalho, L. (2023).
#' \emph{Deviance matrix factorization.}
#' Electronic Journal of Statistics, 17(2): 3762-3810.
#'
#' @examples
#' ...
#'
#' @export sgdgmf.fit
sgdgmf.fit = function (
    Y,
    X = NULL,
    Z = NULL,
    family = poisson(),
    ncomp = 2,
    weights = NULL,
    offset = NULL,
    method = c("airwls", "newton", "msgd", "csgd", "rsgd", "bsgd"),
    penalty = list(),
    control.init = list(),
    control.alg = list()
) {

  # Check and set the model matrices
  Y = set.mat.Y(Y)
  X = set.mat.X(X, nrow(Y), "X")
  Z = set.mat.X(Z, ncol(Y), "Z")

  # Set the model dimensions
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  q = ncol(Z)

  # Check the model family
  family = set.family(family)
  familyname = family$family
  linkname = family$link

  # Check the optimization method
  method = match.arg(method)

  # Check the penalty terms
  lambda = do.call("set.penalty", penalty)

  # Check the control parameters for the initialization
  control.init = do.call("set.control.init", control.init)

  # Check the control parameters for the optimization
  if (method == "airwls") control.alg = do.call("set.control.airwls", control.alg)
  if (method == "newton") control.alg = do.call("set.control.newton", control.alg)
  if (method == "msgd") control.alg = do.call("set.control.msgd", control.alg)
  if (method == "csgd") control.alg = do.call("set.control.csgd", control.alg)
  if (method == "rsgd") control.alg = do.call("set.control.rsgd", control.alg)
  if (method == "bsgd") control.alg = do.call("set.control.bsgd", control.alg)

  alg = control.alg

  # Initialize the parameters
  time.init = proc.time()
  init = init.gmf.param(
    Y = Y, X = X, Z = Z, ncomp = ncomp,
    family = family, method = control.init$method, type = control.init$type,
    niter = control.init$niter, values = control.init$values,
    verbose = control.init$verbose, parallel = control.init$parallel,
    nthreads = control.init$threads, savedata = FALSE)
  time.init = as.numeric(proc.time() - time.init)[3]

  # Select the correct estimation method
  time.optim = proc.time()
  if (method == "airwls") {
    # AIRWLS algorithm
    fit = cpp.fit.airwls(
      Y = Y, X = X, B = init$B, A = init$A, Z = Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = alg$maxiter, nstep = alg$nstep, stepsize = alg$stepsize,
      eps = alg$eps, nafill = alg$nafill, tol = alg$tol,
      damping = alg$damping, verbose = alg$verbose,
      frequency = alg$frequency, parallel = alg$parallel,
      nthreads = alg$nthreads
    )
  }
  if (method == "newton") {
    # Quasi-Newton algorithm
    fit = cpp.fit.newton(
      Y = Y, X = X, B = init$B, A = init$A, Z = Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = alg$maxiter, stepsize = alg$stepsize, eps = alg$eps,
      nafill = alg$nafill, tol = alg$tol, damping = alg$damping,
      verbose = alg$verbose, frequency = alg$frequency,
      parallel = alg$parallel, nthreads = alg$nthreads
    )
  }
  if (method == "msgd") {
    # Memoized SGD algorithm
    fit = cpp.fit.msgd(
      Y = Y, X = X, B = init$B, A = init$A, Z = Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = alg$maxiter, eps = alg$eps, nafill = alg$nafill, tol = alg$tol,
      size = alg$size, burn = alg$burn, rate0 = alg$rate0, decay = alg$decay,
      damping = alg$damping, rate1 = control$rate1, rate2 = alg$rate2,
      parallel = FALSE, nthreads = 1, verbose = alg$verbose,
      frequency = alg$frequency, progress = alg$progress
    )
  }
  if (method == "csgd") {
    # Coordinatewise SGD algorithm
    fit = cpp.fit.csgd(
      Y = Y, X = X, B = init$B, A = init$A, Z = Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = alg$maxiter, eps = alg$eps, nafill = alg$nafill, tol = alg$tol,
      size1 = alg$size[1], size2 = alg$size[2], burn = alg$burn, rate0 = alg$rate0,
      decay = alg$decay, damping = alg$damping, rate1 = alg$rate1, rate2 = alg$rate2,
      parallel = FALSE, nthreads = 1, verbose = alg$verbose,
      frequency = alg$frequency, progress = alg$progress
    )
  }
  if (method == "rsgd") {
    # Rowwise SGD algorithm
    fit = cpp.fit.rsgd(
      Y = Y, X = X, B = init$B, A = init$A, Z = Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = alg$maxiter, eps = alg$eps, nafill = alg$nafill, tol = alg$tol,
      size1 = alg$size[1], size2 = alg$size[2], burn = alg$burn, rate0 = alg$rate0,
      decay = alg$decay, damping = alg$damping, rate1 = alg$rate1, rate2 = alg$rate2,
      parallel = FALSE, nthreads = 1, verbose = alg$verbose,
      frequency = alg$frequency, progress = alg$progress
    )
  }
  if (method == "bsgd") {
    # Blockwise SGD algorithm
    fit = cpp.fit.bsgd(
      Y = Y, X = X, B = init$B, A = init$A, Z = Z, U = init$U, V = init$V,
      familyname = familyname, linkname = linkname, ncomp = ncomp, lambda = lambda,
      maxiter = alg$maxiter, eps = alg$eps, nafill = alg$nafill, tol = alg$tol,
      size1 = alg$size[1], size2 = alg$size[2], burn = alg$burn, rate0 = alg$rate0,
      decay = alg$decay, damping = alg$damping, rate1 = alg$rate1, rate2 = alg$rate2,
      parallel = FALSE, nthreads = 1, verbose = alg$verbose,
      frequency = alg$frequency, progress = alg$progress
    )
  }
  time.optim = as.numeric(proc.time() - time.optim)[3]
  time.tot = time.init + time.optim
  exe.time = c(init = time.init, optim = time.optim, tot = time.tot)

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

  # Set the optimization history
  colnames(fit$trace) = c("iter", "dev", "pen", "pdev", "change", "time")

  # Set the overdispersion parameter of a Negative Binomial family
  if (family$family == "negbinom" | substring(family$family, 1, 17) == "Negative Binomial") {
    family = MASS::negative.binomial(theta = fit$phi)
    family = set.family(family)
    family$theta = fit$phi
  }

  # Output list
  out = list()
  out$method = method
  out$family = family
  out$ncomp = ncomp
  out$npar = m*p + n*q + (n+m)*ncomp
  out$control.init = control.init
  out$control.alg = control.alg
  out$control.cv = list()
  out$Y = Y
  out$X = X
  out$Z = Z
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
  out$cbic = fit$deviance + 2 * df * log(log(nm))
  out$exe.time = exe.time
  out$trace = as.data.frame(fit$trace)
  out$summary.cv = data.frame()

  # Normalize the latent factors
  if (alg$normalize) {
    uv = normalize.uv(out$U, out$V, method = "qr")
    out$U = uv$U
    out$V = uv$V
  }

  # Set the S3 class of the output model
  class(out) = "sgdgmf"

  # Return the estimated model
  return (out)
}
