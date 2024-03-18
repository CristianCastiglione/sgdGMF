
#' @title Check and set the response matrix Y
#'
#' @description
#' Check if the input response matrix is well-defined and return the same
#' matrix without attributes such as row and column names.
#'
#' @keywords internal
set.mat.Y = function (Y) {
  if (!is.numeric(Y)) stop("Y is not numeric.")
  if (!is.matrix(Y)) stop("Y is not a matrix.")
  Y = matrix(c(Y), nrow = nrow(Y), ncol = ncol(Y))
  return (Y)
}

#' @title Check and set the covariate matrix X
#'
#' @description
#' Check if the input covariate matrix X is well-defined and return the same
#' matrix without attributes such as row and column names.
#'
#' @keywords internal
set.mat.X = function (X, n, mat = "X") {
  if (!is.null(X)) {
    if (!is.numeric(X)) stop(paste(mat, "is not numeric."))
    if (!is.matrix(X)) stop(paste(mat, "is not a matrix."))
    if (nrow(X) != n) stop(paste("The dimensions of", mat, "are not compatible with Y."))
    if (anyNA(X)) stop(paste(mat, "contains some NA."))
    if (sum(apply(X, 2, sd) == 0) > 1) stop(paste(mat, "has too many constant columns."))
    X = matrix(c(X), nrow = nrow(X), ncol = ncol(X))
  } else {
    X = matrix(1, nrow = n, ncol = 1)
  }
  return (X)
}

#' @title Check and set the model family
#'
#' @description
#' Check if the model family is allowed and return it eventually with a
#' different family name for compatibility with the \code{C++} implementation
#'
#' @param family a \code{glm} family (see \code{\link{family}} for more details)
#'
#' @export
set.family = function (family) {

  flag = TRUE
  if (check.class(family, "family")) {
    # Gaussian family
    if (family$family == "gaussian") {
      if (family$link == "identity") {
        family$transform = function (y) y
        flag = FALSE
      }
    }
    # Binomial family
    if (family$family %in% c("binomial", "quasibinomial")) {
      if (family$link %in% c("logit", "probit", "cauchit", "cloglog")) {
        family$transform = function (y) jitter(2 * y - 1, amount = 0.25)
        flag = FALSE
      }
    }
    # Poisson family
    if (family$family %in% c("poisson", "quasipoisson")) {
      if (family$link == "log") {
        family$transform = function (y) family$linkfun(y + 0.1)
        flag = FALSE
      }
    }
    # Gamma family
    if (family$family %in% c("gamma", "Gamma")) {
      if (family$link %in% c("inverse", "log", "sqrt")) {
        family$family = "gamma"
        family$transform = function (y) family$linkfun(y)
        flag = FALSE
      }
    }
    # Negative Binomial family
    if (family$family == "negbinom" | substring(family$family, 1, 17) == "Negative Binomial") {
      if (family$link %in% c("inverse", "log", "sqrt")) {
        family$family = "negbinom"
        family$transform = function (y) family$linkfun(y + (y == 0) / 6)
        flag = FALSE
      }
    }
    attr(family$transform, "srcref") = NULL
  }
  if (flag) {
    stop("Family not available")
  }

  return (family)
}


#' @title Set the jittering function for mapping the data (DEPRECATED!)
#'
#' @description
#' Return a function f which maps the original data y into the perturbed data
#' z = f(y) such that z is a real number having the same scale of eta = linkinv(mu).
#' This is useful for initialization purposes.
#'
#' @keywords internal
set.jitter = function (family) {
  f = NULL
  if (family$family %in% c("gaussian", "quasi")) {
    f = function(x) x
  }
  if (family$family %in% c("binomial", "quasibinomial")) {
    f = function(x) jitter(2 * x - 1, amount = 0.25)
  }
  if (family$family %in% c("poisson", "quasipoisson")) {
    f = function(x) family$linkfun(x + 0.1)
  }
  if (family$family %in% c("Gamma", "inverse.gaussian")) {
    f = function(x) family$linkfun(x)
  }
  if (family$family == "negbinom" |
      family$family == "Negative Binomial" |
      substring(family$family, first = 1, last = 17) == "Negative Binomial") {
    f = function(x) family$linkfun(x + (x == 0) / 6)
  }
  return (f)
}

#' @title Check and set the penalty parameters
#'
#' @description
#' Check if the input penalty parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined penalty parameters.
#'
#' @param B penalty parameter of \code{B}
#' @param A penalty parameter of \code{A}
#' @param U penalty parameter of \code{U}
#' @param V penalty parameter of \code{V}
#'
#' @export
set.penalty = function (B = 0, A = 0, U = 1, V = 0) {
  penalty = c(0, 0, 1, 0)
  if (is.numeric(B) && B >= .0) penalty[1] = B
  if (is.numeric(A) && A >= .0) penalty[2] = A
  if (is.numeric(U) && U >= .0) penalty[3] = U
  if (is.numeric(V) && V >= .0) penalty[4] = V
  if (sum(abs(penalty[3:4])) == 0) penalty[3] = 1.
  return (penalty)
}

#' @title Check and set the initialization parameters
#'
#' @description
#' Check if the input initialization parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined initialization parameters.
#'
#' @param method initialization method
#' @param type residual type to be decomposed in the \code{"glm"} method
#' @param values list of custom initialization parameters fixed by the user
#' @param niter number if refinement iterations in the \code{"svd"} method
#' @param normalize if \code{TRUE}, normalize \code{U} and \code{V} to orthogonal \code{U} and lower triangular \code{V}
#' @param verbose if \code{TRUE}, print the initialization state
#' @param parallel if \code{TRUE}, use parallel computing for the \code{"glm"} method
#' @param nthreads number of cores to be used in the \code{"glm"} method
#'
#' @details
#' It takes in input a list of options which define how to initialize the GMF model
#' parameters. It returns a list of safe initialization parameters containing the
#' following elements:
#'
#' @export
set.control.init = function (
    method = c("svd", "glm", "random", "values"),
    type = c("deviance", "pearson", "working"),
    values = list(),
    niter = 5,
    normalize = TRUE,
    verbose = FALSE,
    parallel = FALSE,
    nthreads = 1
) {
  ctr = list()
  ctr$method = match.arg(method)
  ctr$type = match.arg(type)
  ctr$values = list()
  ctr$niter = 5
  ctr$normalize = TRUE
  ctr$verbose = FALSE
  ctr$parallel = FALSE
  ctr$nthreads = 1

  if (is.numeric(niter) && niter > 0) ctr$niter = floor(niter)
  if (is.logical(normalize)) ctr$normalize = normalize
  if (is.logical(verbose)) ctr$verbose = verbose
  if (is.logical(parallel)) ctr$parallel = parallel
  if (is.numeric(nthreads) && nthreads > 0) ctr$nthreads = floor(nthreads)

  if (is.list(values)) {
    if (length(values) > 0) {
      if (c("B", "A", "U", "V") %in% names(values)) {
        if (is.numeric(values$B) && is.matrix(values$B)) ctr$values$B = values$B
        if (is.numeric(values$A) && is.matrix(values$A)) ctr$values$A = values$A
        if (is.numeric(values$U) && is.matrix(values$U)) ctr$values$U = values$U
        if (is.numeric(values$V) && is.matrix(values$V)) ctr$values$V = values$V
      } else {
        values = list()
      }
    }
  }

  return (ctr)
}

#' @title Check and set the control parameters for the AIRWLS algorithm
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @param normalize if \code{TRUE}, normalize \code{U} and \code{V} to uncorrelated Gaussian \code{U} and upper triangular \code{V} with positive diagonal
#' @param maxiter maximum number of iterations
#' @param nstep number of IRWLS steps in each inner loop of AIRWLS
#' @param stepsize step-size parameter scaling each IRWLS step
#' @param eps how much shrinkage has to be introduced on extreme predictions lying outside of the data range
#' @param nafill how frequently the \code{NA} values are filled, by default \code{NA} values are filled at each iteration of the algorithm
#' @param tol tolerance threshold for the stopping criterion
#' @param damping regularization parameter which is added to the diagonal of the Hessian to ensure numerical stability
#' @param verbose if \code{TRUE}, print the optimization status (default \code{TRUE})
#' @param frequency how often the optimization status is printed (only if \code{verbose=TRUE})
#' @param parallel if \code{TRUE}, allows for parallel computing using the \code{C++} library \code{OpenMP}
#' @param nthreads number of cores to be used in parallel (only if \code{parallel=TRUE})
#'
#' @export
set.control.airwls = function (
    normalize = TRUE,
    maxiter = 100,
    nstep = 1,
    stepsize = 0.1,
    eps = 1e-08,
    nafill = 1,
    tol = 1e-05,
    damping = 1e-03,
    verbose = TRUE,
    frequency = 10,
    parallel = FALSE,
    nthreads = 1
) {
  # Set the default control parameters
  ctr = list()
  ctr$normalize = TRUE
  ctr$maxiter = 100
  ctr$nstep = 1
  ctr$stepsize = 0.1
  ctr$eps = 1e-08
  ctr$nafill = 1
  ctr$tol = 1e-05
  ctr$damping = 1e-03
  ctr$verbose = TRUE
  ctr$frequency = 10
  ctr$parallel = FALSE
  ctr$nthreads = 1

  # Standard safety checks
  if (is.logical(normalize)) ctr$normalize = normalize
  if (is.numeric(maxiter) && maxiter >= 1) ctr$maxiter = floor(maxiter)
  if (is.numeric(nstep) && nstep >= 1) ctr$nstep = floor(nstep)
  if (is.numeric(stepsize) && stepsize > 0) ctr$stepsize = stepsize
  if (is.numeric(eps) && eps > 0) ctr$eps = eps
  if (is.numeric(nafill) && nafill >= 1) ctr$nafill = floor(nafill)
  if (is.numeric(tol) && tol > 0) ctr$tol = tol
  if (is.numeric(damping) && damping >= 0) ctr$damping = damping
  if (is.logical(verbose)) ctr$verbose = verbose
  if (is.numeric(frequency) && frequency >= 1) ctr$frequency = floor(frequency)
  if (is.logical(parallel)) ctr$parallel = parallel
  if (is.numeric(nthreads) && nthreads >= 1) ctr$nthreads = floor(nthreads)

  # Additional safety checks
  if (ctr$stepsize > 1) ctr$stepsize = 1
  if (ctr$eps > 1e-01) ctr$eps = 1e-01
  if (ctr$frequency > ctr$maxiter) ctr$frequency = ctr$maxiter

  # Return the checked parameters
  return (ctr)
}


#' @title Check and set the control parameters for the Newton algorithm
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @param normalize if \code{TRUE}, normalize \code{U} and \code{V} to uncorrelated Gaussian \code{U} and upper triangular \code{V} with positive diagonal
#' @param maxiter maximum number of iterations
#' @param stepsize step-size parameter scaling each IRWLS step
#' @param eps how much shrinkage has to be introduced on extreme predictions lying outside of the data range
#' @param nafill how frequently the \code{NA} values are filled, by default \code{NA} values are filled at each iteration of the algorithm
#' @param tol tolerance threshold for the stopping criterion
#' @param damping regularization parameter which is added to the Hessian to ensure numerical stability
#' @param verbose if \code{TRUE}, print the optimization status
#' @param frequency how often the optimization status is printed (only if \code{verbose=TRUE}
#' @param parallel if \code{TRUE}, allows for parallel computing using the \code{C++} library \code{OpenMP}
#' @param nthreads number of cores to be used in parallel (only if \code{parallel=TTUE})
#'
#' @export
set.control.newton = function (
    normalize = TRUE,
    maxiter = 500,
    stepsize = 0.01,
    eps = 1e-08,
    nafill = 1,
    tol = 1e-05,
    damping = 1e-03,
    verbose = TRUE,
    frequency = 50,
    parallel = FALSE,
    nthreads = 1
) {
  # Set the default control parameters
  ctr = list()
  ctr$normalize = TRUE
  ctr$maxiter = 500
  ctr$stepsize = 0.01
  ctr$eps = 1e-08
  ctr$nafill = 1
  ctr$tol = 1e-05
  ctr$damping = 1e-03
  ctr$verbose = TRUE
  ctr$frequency = 50
  ctr$parallel = FALSE
  ctr$nthreads = 1

  # Standard safety checks
  if (is.logical(normalize)) ctr$normalize = normalize
  if (is.numeric(maxiter) && maxiter >= 1) ctr$maxiter = floor(maxiter)
  if (is.numeric(stepsize) && stepsize > 0) ctr$stepsize = stepsize
  if (is.numeric(eps) && eps > 0) ctr$eps = eps
  if (is.numeric(nafill) && nafill >= 1) ctr$nafill = floor(nafill)
  if (is.numeric(tol) && tol > 0) ctr$tol = tol
  if (is.numeric(damping) && damping > 0) ctr$damping = damping
  if (is.logical(verbose)) ctr$verbose = verbose
  if (is.numeric(frequency) && frequency >= 1) ctr$frequency = floor(frequency)
  if (is.logical(parallel)) ctr$parallel = parallel
  if (is.numeric(nthreads) && nthreads >= 1) ctr$nthreads = floor(nthreads)

  # Additional consistency checks
  if (ctr$stepsize > 1) ctr$stepsize = 1
  if (ctr$eps > 1e-01) ctr$eps = 1e-01
  if (ctr$frequency > ctr$maxiter) ctr$frequency = ctr$maxiter

  # Return the checked parameters
  return (ctr)
}


#' @title Check and set the control parameters for the memoized-SGD algorithm
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @export
set.control.msgd = function (
    normalize = TRUE,
    maxiter = 100,
    eps = 1e-08,
    nafill = 10,
    tol = 1e-05,
    size = 100,
    burn = 90,
    rate0 = 0.01,
    decay = 1.0,
    damping = 1e-04,
    rate1 = 0.05,
    rate2 = 0.01,
    parallel = FALSE,
    nthreads = 1,
    verbose = TRUE,
    frequency = 10,
    progress = FALSE
) {
  # Set the default control parameters
  ctr = list()
  ctr$normalize = TRUE
  ctr$maxiter = 100
  ctr$eps = 1e-08
  ctr$nafill = 10
  ctr$tol = 1e-05
  ctr$size = 100
  ctr$burn = 90
  ctr$rate0 = 0.01
  ctr$decay = 1.0
  ctr$damping = 1e-04
  ctr$rate1 = 0.05
  ctr$rate2 = 0.01
  ctr$parallel = FALSE
  ctr$nthreads = 1
  ctr$verbose = TRUE
  ctr$frequency = 10
  ctr$progress = FALSE

  # Standard safety checks
  if (is.logical(normalize)) ctr$normalize = normalize
  if (is.numeric(maxiter) && maxiter >= 1) ctr$maxiter = floor(maxiter)
  if (is.numeric(eps) && eps > 0) ctr$eps = eps
  if (is.numeric(nafill) && nafill >= 1) ctr$nafill = floor(nafill)
  if (is.numeric(tol) && tol > 0) ctr$tol = tol
  if (is.numeric(size) && size >= 1) ctr$size = floor(size)
  if (is.numeric(burn) && burn > 0 && burn <= 1) ctr$burn = burn
  if (is.numeric(rate0) && rate0 > 0) ctr$rate0 = rate0
  if (is.numeric(decay) && decay > 0) ctr$decay = decay
  if (is.numeric(damping) && damping > 0) ctr$damping = damping
  if (is.numeric(rate1) && rate1 > 0) ctr$rate1 = rate1
  if (is.numeric(rate2) && rate2 > 0) ctr$rate2 = rate2
  if (is.logical(parallel)) ctr$parallel = parallel
  if (is.numeric(nthreads) && nthreads >= 1) ctr$nthreads = floor(nthreads)
  if (is.logical(verbose)) ctr$verbose = verbose
  if (is.numeric(frequency) && frequency >= 1) ctr$frequency = floor(frequency)
  if (is.logical(progress)) ctr$progress = progress

  # Addional consistency checks
  if (ctr$nafill > ctr$maxiter) ctr$nafill = ctr$maxiter
  if (ctr$eps > 1e-01) ctr$eps = 1e-01
  if (ctr$rate1 > 1 - 1e-08) ctr$rate1 = 1 - 1e-08
  if (ctr$rate2 > 1 - 1e-08) ctr$rate2 = 1 - 1e-08
  if (ctr$frequency > ctr$maxiter) ctr$frequency = ctr$maxiter

  # Return the checked parameters
  return (ctr)
}

#' @title Check and set the control parameters for the coordinate-SGD algorithm
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @param normalize if \code{TRUE}, normalize \code{U} and \code{V} to uncorrelated Gaussian \code{U} and upper triangular \code{V} with positive diagonal
#' @param maxiter maximum number of iterations
#' @param eps how much shrinkage has to be introduced on extreme predictions lying outside of the data range
#' @param nafill how frequently the \code{NA} values are filled, by default \code{NA} values are filled at each iteration of the algorithm
#' @param tol tolerance threshold for the stopping criterion
#' @param size mini-batch size, the first value is for row sub-sample, the second value is for column sub-sample
#' @param burn percentage of iterations to ignore before performing Polyak averaging
#' @param rate0 initial learning rate
#' @param decay learning rate decay
#' @param damping regularization parameter which is added to the Hessian to ensure numerical stability
#' @param rate1 exponential decay rate for the moment estimate of the gradient
#' @param rate2 exponential decay rate for the moment estimate of the Hessian
#' @param verbose if \code{TRUE}, print the optimization status
#' @param frequency how often the optimization status is printed (only if \code{verbose=TRUE})
#' @param progress if \code{TRUE}, print a compact progress-bar instead of a full-report of the optimization status (only if \code{verbose=TRUE})
#'
#' @export
set.control.csgd = function (
    normalize = TRUE,
    maxiter = 1000,
    eps = 1e-08,
    nafill = 10,
    tol = 1e-08,
    size = c(100, 100),
    burn = 1,
    rate0 = 0.01,
    decay = 0.01,
    damping = 1e-03,
    rate1 = 0.1,
    rate2 = 0.01,
    verbose = TRUE,
    frequency = 250,
    progress = FALSE
) {
  # Set the default control parameters
  ctr = list()
  ctr$normalize = TRUE
  ctr$maxiter = 1000
  ctr$eps = 1e-08
  ctr$nafill = 10
  ctr$tol = 1e-08
  ctr$size = c(100, 100)
  ctr$burn = 1
  ctr$rate0 = 0.01
  ctr$decay = 0.01
  ctr$damping = 1e-03
  ctr$rate1 = 0.1
  ctr$rate2 = 0.01
  ctr$verbose = TRUE
  ctr$frequency = 250
  ctr$progress = FALSE

  # Standard safety checks
  if (is.logical(normalize)) ctr$normalize = normalize
  if (is.numeric(maxiter) && maxiter >= 1) ctr$maxiter = floor(maxiter)
  if (is.numeric(eps) && eps > 0) ctr$eps = eps
  if (is.numeric(nafill) && nafill >= 1) ctr$nafill = floor(nafill)
  if (is.numeric(tol) && tol > 0) ctr$tol = tol
  if (is.numeric(size) && all(size >= 1)) ctr$size = floor(size[1:2])
  if (is.numeric(burn) && burn > 0 && burn <= 1) ctr$burn = burn
  if (is.numeric(rate0) && rate0 > 0) ctr$rate0 = rate0
  if (is.numeric(decay) && decay > 0) ctr$decay = decay
  if (is.numeric(damping) && damping > 0) ctr$damping = damping
  if (is.numeric(rate1) && rate1 > 0) ctr$rate1 = rate1
  if (is.numeric(rate2) && rate2 > 0) ctr$rate2 = rate2
  if (is.logical(verbose)) ctr$verbose = verbose
  if (is.numeric(frequency) && frequency >= 1) ctr$frequency = floor(frequency)
  if (is.logical(progress)) ctr$progress = progress

  # Additional consistency checks
  if (ctr$nafill > ctr$maxiter) ctr$nafill = ctr$maxiter
  if (ctr$eps > 1e-01) ctr$eps = 1e-01
  if (ctr$rate1 > 1 - 1e-08) ctr$rate1 = 1 - 1e-08
  if (ctr$rate2 > 1 - 1e-08) ctr$rate2 = 1 - 1e-08
  if (ctr$frequency > ctr$maxiter) ctr$frequency = ctr$maxiter

  # Return the corrected control parameters
  return (ctr)
}

#' @title Check and set the control parameters for the rowwise-SGD algorithm
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @export
set.control.rsgd = function (
    normalize = TRUE,
    maxiter = 100,
    eps = 1e-08,
    nafill = 10,
    tol = 1e-05,
    size = 100,
    burn = 90,
    rate0 = 0.01,
    decay = 1.0,
    damping = 1e-04,
    rate1 = 0.05,
    rate2 = 0.01,
    verbose = TRUE,
    frequency = 10,
    progress = FALSE
) {
  # Set the default control parameters
  ctr = list()
  ctr$normalize = TRUE
  ctr$maxiter = 100
  ctr$eps = 1e-08
  ctr$nafill = 10
  ctr$tol = 1e-05
  ctr$size = 100
  ctr$burn = 90
  ctr$rate0 = 0.01
  ctr$decay = 1.0
  ctr$damping = 1e-04
  ctr$rate1 = 0.05
  ctr$rate2 = 0.01
  ctr$verbose = TRUE
  ctr$frequency = 10
  ctr$progress = FALSE

  # Standard safety checks
  if (is.logical(normalize)) ctr$normalize = normalize
  if (is.numeric(maxiter) && maxiter >= 1) ctr$maxiter = floor(maxiter)
  if (is.numeric(eps) && eps > 0) ctr$eps = eps
  if (is.numeric(nafill) && nafill >= 1) ctr$nafill = floor(nafill)
  if (is.numeric(tol) && tol > 0) ctr$tol = tol
  if (is.numeric(size) && size >= 1) ctr$size = floor(size)
  if (is.numeric(burn) && burn > 0 && burn <= 1) ctr$burn = burn
  if (is.numeric(rate0) && rate0 > 0) ctr$rate0 = rate0
  if (is.numeric(decay) && decay > 0) ctr$decay = decay
  if (is.numeric(damping) && damping > 0) ctr$damping = damping
  if (is.numeric(rate1) && rate1 > 0) ctr$rate = rate1
  if (is.numeric(rate2) && rate2 > 0) ctr$rate = rate2
  if (is.logical(verbose)) ctr$verbose = verbose
  if (is.numeric(frequency) && frequency >= 1) ctr$frequency = floor(frequency)
  if (is.logical(progress)) ctr$progress = progress

  # Additional consistency checks
  if (ctr$nafill > ctr$maxiter) ctr$nafill = ctr$maxiter
  if (ctr$eps > 1e-01) ctr$eps = 1e-01
  if (ctr$rate1 > 1 - 1e-08) ctr$rate1 = 1 - 1e-08
  if (ctr$rate2 > 1 - 1e-08) ctr$rate2 = 1 - 1e-08
  if (ctr$frequency > ctr$maxiter) ctr$frequency = ctr$maxiter

  # Return the corrected control parameters
  return (ctr)
}

#' @title Check and set the control parameters for the blockwise-SGD algorithm
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @param normalize if \code{TRUE}, normalize \code{U} and \code{V} to uncorrelated Gaussian \code{U} and upper triangular \code{V} with positive diagonal
#' @param maxiter maximum number of iterations
#' @param eps how much shrinkage has to be introduced on extreme predictions lying outside of the data range
#' @param nafill how frequently the \code{NA} values are filled, by default \code{NA} values are filled at each iteration of the algorithm
#' @param tol tolerance threshold for the stopping criterion
#' @param size mini-batch size, the first value is for row sub-sample, the second value is for column sub-sample
#' @param burn percentage of iterations to ignore before performing Polyak averaging
#' @param rate0 initial learning rate
#' @param decay learning rate decay
#' @param damping regularization parameter which is added to the Hessian to ensure numerical stability
#' @param rate1 exponential decay rate for the moment estimate of the gradient
#' @param rate2 exponential decay rate for the moment estimate of the Hessian
#' @param verbose if \code{TRUE}, print the optimization status
#' @param frequency how often the optimization status is printed (only if \code{verbose=TRUE})
#' @param progress if \code{TRUE}, print a compact progress-bar instead of a full-report of the optimization status (only if \code{verbose=TRUE})
#'
#' @export
set.control.bsgd = function (
    normalize = TRUE,
    maxiter = 1000,
    eps = 1e-08,
    nafill = 10,
    tol = 1e-08,
    size = c(100, 100),
    burn = 1,
    rate0 = 0.01,
    decay = 0.01,
    damping = 1e-03,
    rate1 = 0.1,
    rate2 = 0.01,
    verbose = TRUE,
    frequency = 250,
    progress = FALSE
) {
  # Set the default control parameters
  ctr = list()
  ctr$normalize = TRUE
  ctr$maxiter = 1000
  ctr$eps = 1e-08
  ctr$nafill = 10
  ctr$tol = 1e-08
  ctr$size = c(100, 100)
  ctr$burn = 1
  ctr$rate0 = 0.01
  ctr$decay = 0.01
  ctr$damping = 1e-03
  ctr$rate1 = 0.1
  ctr$rate2 = 0.01
  ctr$verbose = TRUE
  ctr$frequency = 250
  ctr$progress = FALSE

  # Standard safety checks
  if (is.logical(normalize)) ctr$normalize = normalize
  if (is.numeric(maxiter) && maxiter >= 1) ctr$maxiter = floor(maxiter)
  if (is.numeric(eps) && eps > 0) ctr$eps = eps
  if (is.numeric(nafill) && nafill >= 1) ctr$nafill = floor(nafill)
  if (is.numeric(tol) && tol > 0) ctr$tol = tol
  if (is.numeric(size) & all(size >= 1)) ctr$size = floor(size[1:2])
  if (is.numeric(burn) && burn > 0 && burn <= 1) ctr$burn = burn
  if (is.numeric(rate0) && rate0 > 0) ctr$rate0 = rate0
  if (is.numeric(decay) && decay > 0) ctr$decay = decay
  if (is.numeric(damping) && damping > 0) ctr$damping = damping
  if (is.numeric(rate1) && rate1 > 0) ctr$rate1 = rate1
  if (is.numeric(rate2) && rate2 > 0) ctr$rate2 = rate2
  if (is.logical(verbose)) ctr$verbose = verbose
  if (is.numeric(frequency) && frequency >= 1) ctr$frequency = floor(frequency)
  if (is.logical(progress)) ctr$progress = progress

  # Additional consistency checks
  if (ctr$nafill > ctr$maxiter) ctr$nafill = ctr$maxiter
  if (ctr$eps > 1e-01) ctr$eps = 1e-01
  if (ctr$rate1 > 1 - 1e-08) ctr$rate1 = 1 - 1e-08
  if (ctr$rate2 > 1 - 1e-08) ctr$rate2 = 1 - 1e-08
  if (ctr$frequency > ctr$maxiter) ctr$frequency = ctr$maxiter

  # Return the check control parameters
  return (ctr)
}

#' @title Check and set the control parameters for the select optimization algorithm
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @param method optimization method: \code{"airwls"}, \code{"newton"}, \code{"msgd"}, \code{"csgd"}, \code{"rsgd"}, \code{"bsgd"}
#' @param control list of algorithm-specific control parameters
#'
#' @details
#' It is not necessary to provide a complete list of control parameters, one can
#' just specify a list containing the parameters he/she needs to change from the
#' default values. Wrongly specified parameters are ignored or set to default values.
#' For a detailed description of all the algorithm-specific control parameters,
#' please refer to
#' \code{\link{set.control.airwls}} (\code{method="airwls"}),
#' \code{\link{set.control.newton}} (\code{method="newton"}),
#' \code{\link{set.control.csgd}} (\code{method="csgd"}),
#' \code{\link{set.control.bsgd}} (\code{method="bsgd"}).
#'
#' @export
set.control.alg = function (
    method = c("airwls", "newton", "msgd", "csgd", "rsgd", "bsgd"),
    control = list()
) {
  # Check the optimization method
  method = match.arg(method)

  # Set the input values for the control parameters, check whether the
  # they are allowed and, if they aren't, set them to default values
  if (method == "airwls") control = do.call("set.control.airwls", control)
  if (method == "newton") control = do.call("set.control.newton", control)
  if (method == "msgd") control = do.call("set.control.msgd", control)
  if (method == "csgd") control = do.call("set.control.csgd", control)
  if (method == "rsgd") control = do.call("set.control.rsgd", control)
  if (method == "bsgd") control = do.call("set.control.bsgd", control)

  # Return the corrected control parameters
  return (control)
}


#' @title Check and set the cross-validation parameters
#'
#' @description
#' Check if the input cross-validation parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined cross-validation parameters.
#'
#' @param nfolds number of cross-validation folds
#' @param parallel if \code{TRUE}, allows for parallel computing
#' @param nthreads number of cores to use in parallel (only if \code{parallel=TRUE})
#'
#' @export
set.control.cv = function (nfolds = 5, parallel = FALSE, nthreads = 1) {
  ctr = list(nfolds = 5, parallel = FALSE, nthreads = 1)
  if (is.numeric(nfolds) && nfolds >= 1) ctr$nfolds = floor(nfolds)
  if (is.logical(parallel)) ctr$parallel = parallel
  if (is.numeric(nthreads) && nthreads >= 1) ctr$nthreads = floor(nthreads)
  return (ctr)
}

