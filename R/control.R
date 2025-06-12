
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
set.mat.X = function (X, n, m) {
  if (!is.null(X)) {
    if (!is.numeric(X)) stop("X is not numeric.")
    if (!is.matrix(X)) stop("X is not a matrix.")
    if (nrow(X) != n) stop("The dimensions of X are not compatible with Y.")
    if (anyNA(X)) stop("X contains some NA.")
    if (sum(apply(X, 2, stats::sd) == 0) > 1) stop("X has too many constant columns.")
    X = matrix(c(X), nrow = nrow(X), ncol = ncol(X))
  } else {
    X = matrix(1, nrow = n, ncol = 1)
  }
  return (X)
}

#' @title Check and set the covariate matrix X
#'
#' @description
#' Check if the input covariate matrix X is well-defined and return the same
#' matrix without attributes such as row and column names.
#'
#' @keywords internal
set.mat.Z = function (Z, n, m) {
  if (!is.null(Z)) {
    if (!is.numeric(Z)) stop("Z is not numeric.")
    if (!is.matrix(Z)) stop("Z is not a matrix.")
    if (nrow(Z) != m) stop("The dimensions of Z are not compatible with Y.")
    if (anyNA(Z)) stop("Z contains some NA.")
    if (sum(apply(Z, 2, stats::sd) == 0) > 1) stop("Z has too many constant columns.")
    Z = matrix(c(Z), nrow = nrow(Z), ncol = ncol(Z))
  } else {
    Z = matrix(1, nrow = m, ncol = 1)
  }
  return (Z)
}

#' @title Check and set the weighting matrix
#'
#' @description
#' Check if the input weighting matrix is well-defined and return the same
#' matrix without attributes such as row and column names.
#'
#' @keywords internal
set.mat.weights = function (W, n, m) {
  if (!is.null(W)) {
    if (!is.numeric(W)) stop("weights is not numeric.")
    if (!is.matrix(W)) stop("weights is not a matrix.")
    if (nrow(W) != n) stop("The dimensions of weights are not compatible with Y.")
    if (ncol(W) != m) stop("The dimensions of weights are not compatible with Y.")
    if (anyNA(W)) stop("weights contains some NA.")
    if (any(W < 0)) stop("weights contains negative values.")
    if (all(W == 0)) stop("weights contains only 0 values.")
    W = matrix(c(W), nrow = nrow(W), ncol = ncol(W))
  } else {
    W = matrix(1, nrow = n, ncol = m)
  }
  return (W)
}

#' @title Check and set the offset matrix
#'
#' @description
#' Check if the input offset matrix is well-defined and return the same
#' matrix without attributes such as row and column names.
#'
#' @keywords internal
set.mat.offset = function (O, n, m) {
  if (!is.null(O)) {
    if (!is.numeric(O)) stop("weights is not numeric.")
    if (!is.matrix(O)) stop("weights is not a matrix.")
    if (nrow(O) != n) stop("The dimensions of weights are not compatible with Y.")
    if (ncol(O) != m) stop("The dimensions of weights are not compatible with Y.")
    if (anyNA(O)) stop("weights contains some NA.")
    O = matrix(c(O), nrow = nrow(O), ncol = ncol(O))
  } else {
    O = matrix(0, nrow = n, ncol = m)
  }
  return (O)
}

#' @title Check and set the model family
#'
#' @description
#' Check if the model family is allowed and return it eventually with a
#' different family name for compatibility with the \code{C++} implementation
#'
#' @param family a \code{glm} family (see \code{\link{family}} for more details)
#'
#' @keywords internal
set.family = function (family) {

  flag = TRUE
  if (is(family, "family")) {
    # Gaussian family
    if (family$family == "gaussian") {
      if (family$link == "identity") {
        family$varfun = "const"
        family$transform = function (y) y
        flag = FALSE
      }
    }
    # Binomial family
    if (family$family %in% c("binomial", "quasibinomial")) {
      if (family$link %in% c("logit", "probit", "cauchit", "cloglog")) {
        family$varfun = "mu(1-mu)"
        family$transform = function (y) jitter(2 * y - 1, amount = 0.25)
        flag = FALSE
      }
    }
    # Poisson family
    if (family$family %in% c("poisson", "quasipoisson")) {
      if (family$link == "log") {
        family$varfun = "mu"
        family$transform = function (y) family$linkfun(y + 0.1)
        flag = FALSE
      }
    }
    # Gamma family
    if (family$family %in% c("gamma", "Gamma")) {
      if (family$link %in% c("inverse", "log", "sqrt")) {
        family$family = "gamma"
        family$varfun = "mu^2"
        family$transform = function (y) family$linkfun(y)
        flag = FALSE
      }
    }
    # Inverse-Gaussian family
    if (family$family %in% c("invgaussian", "inverse.gaussian")) {
      if (family$link %in% c("inverse", "1/mu^2", "log", "sqrt")) {
        family$family = "invgaussian"
        family$varfun = "mu^3"
        family$transform = function (y) family$linkfun(y)
        flag = FALSE
      }
    }
    # Negative Binomial family
    if (family$family == "negbinom" | substring(family$family, 1, 17) == "Negative Binomial") {
      if (family$link %in% c("inverse", "log", "sqrt")) {
        family$family = "negbinom"
        family$varfun = "mu(1+t*mu)"
        family$transform = function (y) family$linkfun(y + (y == 0) / 6)
        flag = FALSE
      }
    }
    if (family$family == "quasi") {
      family$varfun = ifelse(family$varfun == "constant", "const", family$varfun)
      family$transform = switch(family$link,
        "identity" = function(y) y,
        "inverse" = function(y) family$linkfun(y),
        "log" = function(y) family$linkfun(y + 0.1),
        "sqrt" = function(y) family$linkfun(y + 0.1),
        "logit" = function(y) jitter(2 * y - 1, amount = 0.25),
        "probit" = function(y) jitter(2 * y - 1, amount = 0.25),
        "caichit" = function(y) jitter(2 * y - 1, amount = 0.25),
        "cloglog" = function(y) jitter(2 * y - 1, amount = 0.25))
      flag = FALSE
    }
    attr(family$transform, "srcref") = NULL
  }
  if (flag) {
    stop("Family not available")
  }

  return (family)
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
#' @keywords internal
set.penalty = function (B = 0, A = 0, U = 1, V = 0) {
  penalty = c(0, 0, 1, 0)

  message = function (var)
    warning(paste0("Penalty parameter: '", var, "' was set to default value."),
            call. = FALSE, immediate. = TRUE, domain = NULL)

  if (is.numeric(B) && B >= .0) penalty[1] = B else message("B")
  if (is.numeric(A) && A >= .0) penalty[2] = A else message("A")
  if (is.numeric(U) && U >= .0) penalty[3] = U else message("U")
  if (is.numeric(V) && V >= .0) penalty[4] = V else message("V")
  if (sum(abs(penalty[3:4])) == 0) {penalty[3] = 1.; message("U")}

  return (penalty)
}

#' @title Check and set the initialization parameters for a GMF model
#'
#' @description
#' Check if the input initialization parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined options which specify how
#' to initialize a GMF model. See \code{\link{sgdgmf.init}} for more details upon the methods used for initialisation.
#'
#' @param method initialization method (see \code{\link{sgdgmf.init}} for more details upon the initialization methods used)
#' @param type residual type to be decomposed (see \code{\link{sgdgmf.init}} for more details upon the residuals used)
#' @param values list of custom initialization parameters fixed by the user
#' @param niter number if refinement iterations in the \code{"svd"} method
#' @param normalize if \code{TRUE}, normalize \code{U} and \code{V} to orthogonal \code{U} and lower triangular \code{V}
#' @param verbose if \code{TRUE}, print the initialization state
#' @param parallel if \code{TRUE}, use parallel computing for the \code{"glm"} method
#' @param nthreads number of cores to be used in the \code{"glm"} method
#'
#' @returns A \code{list} of control parameters for the initialization
#'
#' @seealso \code{\link{set.control.alg}}, \code{\link{set.control.cv}}, \code{\link{sgdgmf.init}}
#'
#' @examples
#' library(sgdGMF)
#'
#' # Empty call
#' set.control.init()
#'
#' # Parametrized call
#' set.control.init(method = "glm", type = "deviance", niter = 10)
#'
#' @export set.control.init
set.control.init = function (
    method = c("ols", "glm", "random", "values"),
    type = c("deviance", "pearson", "working", "link"),
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

  message = function (var)
    warning(paste0("Init. control: '", var,"' was set to default value."),
            call. = FALSE, immediate. = TRUE, domain = NULL)

  if (is.numeric(niter) && niter > 0) ctr$niter = floor(niter) else message("niter")
  if (is.logical(normalize)) ctr$normalize = normalize else message("normalize")
  if (is.logical(verbose)) ctr$verbose = verbose else message("verbose")
  if (is.logical(parallel)) ctr$parallel = parallel else message("parallel")
  if (is.numeric(nthreads) && nthreads > 0) ctr$nthreads = floor(nthreads) else message("nthreads")

  if (is.list(values)) {
    if (length(values) > 0) {
      if (all(c("B", "A", "U", "V") %in% names(values))) {
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
#' Check if the input control parameters of the AIRWLS algorithm are allowed
#' and set them to default values if they are not. Returns a list of
#' well-defined control parameters.
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
#' @param savedata if \code{TRUE}, saves a copy of the data and fitted values
#'
#' @returns A \code{list} of control parameters for the AIRWLS algorithm
#'
#' @examples
#' library(sgdGMF)
#'
#' # Empty call
#' set.control.airwls()
#'
#' # Parametrized call
#' set.control.airwls(maxiter = 100, nstep = 5, stepsize = 0.3)
#'
#'
#' @export set.control.airwls
set.control.airwls = function (
    normalize = TRUE,
    maxiter = 100,
    nstep = 1,
    stepsize = 0.1,
    eps = 1e-08,
    nafill = 1,
    tol = 1e-05,
    damping = 1e-03,
    verbose = FALSE,
    frequency = 10,
    parallel = FALSE,
    nthreads = 1,
    savedata = TRUE
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
  ctr$verbose = FALSE
  ctr$frequency = 10
  ctr$parallel = FALSE
  ctr$nthreads = 1
  ctr$savedata = TRUE

  message = function (var)
    warning(paste0("AIRWLS control: '", var,"' was set to default value."),
            call. = FALSE, immediate. = TRUE, domain = NULL)

  # Standard safety checks
  if (is.logical(normalize)) ctr$normalize = normalize else message("normalize")
  if (is.numeric(maxiter) && maxiter >= 1) ctr$maxiter = floor(maxiter) else message("maxiter")
  if (is.numeric(nstep) && nstep >= 1) ctr$nstep = floor(nstep) else message("nstep")
  if (is.numeric(stepsize) && stepsize > 0) ctr$stepsize = stepsize else message("stepsize")
  if (is.numeric(eps) && eps > 0) ctr$eps = eps else message("eps")
  if (is.numeric(nafill) && nafill >= 1) ctr$nafill = floor(nafill) else message("nafill")
  if (is.numeric(tol) && tol > 0) ctr$tol = tol else message("tol")
  if (is.numeric(damping) && damping >= 0) ctr$damping = damping else message("damping")
  if (is.logical(verbose)) ctr$verbose = verbose else message("verbose")
  if (is.numeric(frequency) && frequency >= 1) ctr$frequency = floor(frequency) else message("frequency")
  if (is.logical(parallel)) ctr$parallel = parallel else message("parallel")
  if (is.numeric(nthreads) && nthreads >= 1) ctr$nthreads = floor(nthreads) else message("nthreads")
  if (is.logical(savedata)) ctr$savedata = savedata else message("verbose")

  # Additional safety checks
  if (ctr$stepsize > 1) {ctr$stepsize = 1; message("stepsize")}
  if (ctr$eps > 1e-01) {ctr$eps = 1e-01; message("eps")}
  if (ctr$frequency > ctr$maxiter) {ctr$frequency = ctr$maxiter; message("maxiter")}

  # Set the number of threads
  ncores = parallel::detectCores() - 1
  ctr$nthreads = floor(max(1, min(ctr$nthreads, ncores)))

  # OpenMP check
  if (ctr$parallel & !omp.check()) {
    ctr$parallel = FALSE
    ctr$nthreads = 1
    warning("OpenMP not detected. Parallel computing options are unavailable.",
            call. = FALSE, immediate. = TRUE, domain = NULL)
  }

  # Return the checked parameters
  return (ctr)
}


#' @title Check and set the control parameters for the Newton algorithm
#'
#' @description
#' Check if the input control parameters of the quasi-Newton algorithm  are
#' allowed and set them to default values if they are not. Returns a list of
#' well-defined control parameters.
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
#' @param savedata if \code{TRUE}, saves a copy of the data and fitted values
#'
#' @returns A \code{list} of control parameters for the quasi-Newton algorithm
#'
#' @examples
#' library(sgdGMF)
#'
#' # Empty call
#' set.control.newton()
#'
#' # Parametrized call
#' set.control.newton(maxiter = 1000, stepsize = 0.01, tol = 1e-04)
#'
#' @export set.control.newton
set.control.newton = function (
    normalize = TRUE,
    maxiter = 500,
    stepsize = 0.01,
    eps = 1e-08,
    nafill = 1,
    tol = 1e-05,
    damping = 1e-03,
    verbose = FALSE,
    frequency = 50,
    parallel = FALSE,
    nthreads = 1,
    savedata = TRUE
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
  ctr$verbose = FALSE
  ctr$frequency = 50
  ctr$parallel = FALSE
  ctr$nthreads = 1
  ctr$savedata = TRUE

  message = function (var)
    warning(paste0("Newton control: '", var,"' was set to default value."),
            call. = FALSE, immediate. = TRUE, domain = NULL)

  # Standard safety checks
  if (is.logical(normalize)) ctr$normalize = normalize else message("normalize")
  if (is.numeric(maxiter) && maxiter >= 1) ctr$maxiter = floor(maxiter) else message("maxiter")
  if (is.numeric(stepsize) && stepsize > 0) ctr$stepsize = stepsize else message("stepsize")
  if (is.numeric(eps) && eps > 0) ctr$eps = eps else message("eps")
  if (is.numeric(nafill) && nafill >= 1) ctr$nafill = floor(nafill) else message("nafill")
  if (is.numeric(tol) && tol > 0) ctr$tol = tol else message("tol")
  if (is.numeric(damping) && damping > 0) ctr$damping = damping else message("damping")
  if (is.logical(verbose)) ctr$verbose = verbose else message("verbose")
  if (is.numeric(frequency) && frequency >= 1) ctr$frequency = floor(frequency) else message("frequency")
  if (is.logical(parallel)) ctr$parallel = parallel else message("parallel")
  if (is.numeric(nthreads) && nthreads >= 1) ctr$nthreads = floor(nthreads) else message("nthreads")
  if (is.logical(savedata)) ctr$savedata = savedata else message("verbose")

  # Additional consistency checks
  if (ctr$stepsize > 1) {ctr$stepsize = 1; message("stepsize")}
  if (ctr$eps > 1e-01) {ctr$eps = 1e-01; message("eps")}
  if (ctr$frequency > ctr$maxiter) {ctr$frequency = ctr$maxiter; message("maxiter")}

  # Set the number of threads
  ncores = parallel::detectCores() - 1
  ctr$nthreads = floor(max(1, min(ctr$nthreads, ncores)))

  # OpenMP check
  if (ctr$parallel & !omp.check()) {
    ctr$parallel = FALSE
    ctr$nthreads = 1
    warning("OpenMP not detected. Parallel computing options are unavailable.",
            call. = FALSE, immediate. = TRUE, domain = NULL)
  }
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
#' @param savedata if \code{TRUE}, saves a copy of the data and fitted values
#'
#' @returns A \code{list} of control parameters for the adaptive SGD algorithm with coordinate-wise sub-sampling
#'
#' @examples
#' library(sgdGMF)
#'
#' # Empty call
#' set.control.coord.sgd()
#'
#' # Parametrized call
#' set.control.coord.sgd(maxiter = 2000, rate0 = 0.01, decay = 0.01)
#'
#' @export set.control.coord.sgd
set.control.coord.sgd = function (
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
    verbose = FALSE,
    frequency = 250,
    progress = FALSE,
    savedata = TRUE
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
  ctr$verbose = FALSE
  ctr$frequency = 250
  ctr$progress = FALSE
  ctr$savedata = TRUE

  message = function (var)
    warning(paste0("C-SGD control: '", var,"' was set to default value."),
            call. = FALSE, immediate. = TRUE, domain = NULL)

  # Standard safety checks
  if (is.logical(normalize)) ctr$normalize = normalize else message("normalize")
  if (is.numeric(maxiter) && maxiter >= 1) ctr$maxiter = floor(maxiter) else message("maxiter")
  if (is.numeric(eps) && eps > 0) ctr$eps = eps else message("eps")
  if (is.numeric(nafill) && nafill >= 1) ctr$nafill = floor(nafill) else message("nafill")
  if (is.numeric(tol) && tol > 0) ctr$tol = tol else message("tol")
  if (is.numeric(size) && all(size >= 1)) ctr$size = floor(size[1:2]) else message("size")
  if (is.numeric(burn) && burn > 0 && burn <= 1) ctr$burn = burn else message("burn")
  if (is.numeric(rate0) && rate0 > 0) ctr$rate0 = rate0 else message("rate0")
  if (is.numeric(decay) && decay > 0) ctr$decay = decay else message("decay")
  if (is.numeric(damping) && damping > 0) ctr$damping = damping else message("damping")
  if (is.numeric(rate1) && rate1 > 0) ctr$rate1 = rate1 else message("rate1")
  if (is.numeric(rate2) && rate2 > 0) ctr$rate2 = rate2 else message("rate2")
  if (is.logical(verbose)) ctr$verbose = verbose else message("verbose")
  if (is.numeric(frequency) && frequency >= 1) ctr$frequency = floor(frequency) else message("frequency")
  if (is.logical(progress)) ctr$progress = progress else message("progress")
  if (is.logical(savedata)) ctr$savedata = savedata else message("verbose")

  # Additional consistency checks
  if (ctr$nafill > ctr$maxiter) {ctr$nafill = ctr$maxiter; message("maxiter")}
  if (ctr$eps > 1e-01) {ctr$eps = 1e-01; message("eps")}
  if (ctr$rate1 > 1 - 1e-08) {ctr$rate1 = 1 - 1e-08; message("rate1")}
  if (ctr$rate2 > 1 - 1e-08) {ctr$rate2 = 1 - 1e-08; message("rate2")}
  if (ctr$frequency > ctr$maxiter) {ctr$frequency = ctr$maxiter; message("frequency")}

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
#' @param savedata if \code{TRUE}, saves a copy of the data and fitted values
#'
#' @returns A \code{list} of control parameters for the adaptive SGD algorithm with block-wise sub-sampling
#'
#' @examples
#' library(sgdGMF)
#'
#' # Empty call
#' set.control.block.sgd()
#'
#' # Parametrized call
#' set.control.block.sgd(maxiter = 2000, rate0 = 0.01, decay = 0.01)
#'
#'
#' @export set.control.block.sgd
set.control.block.sgd = function (
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
    verbose = FALSE,
    frequency = 250,
    progress = FALSE,
    savedata = TRUE
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
  ctr$verbose = FALSE
  ctr$frequency = 250
  ctr$progress = FALSE
  ctr$savedata = TRUE

  message = function (var)
    warning(paste0("B-SGD control: '", var,"' was set to default value."),
            call. = FALSE, immediate. = TRUE, domain = NULL)

  # Standard safety checks
  if (is.logical(normalize)) ctr$normalize = normalize else message("normalize")
  if (is.numeric(maxiter) && maxiter >= 1) ctr$maxiter = floor(maxiter) else message("maxiter")
  if (is.numeric(eps) && eps > 0) ctr$eps = eps else message("eps")
  if (is.numeric(nafill) && nafill >= 1) ctr$nafill = floor(nafill) else message("nafill")
  if (is.numeric(tol) && tol > 0) ctr$tol = tol else message("tol")
  if (is.numeric(size) & all(size >= 1)) ctr$size = floor(size[1:2]) else message("size")
  if (is.numeric(burn) && burn > 0 && burn <= 1) ctr$burn = burn else message("burn")
  if (is.numeric(rate0) && rate0 > 0) ctr$rate0 = rate0 else message("rate0")
  if (is.numeric(decay) && decay > 0) ctr$decay = decay else message("decay")
  if (is.numeric(damping) && damping > 0) ctr$damping = damping else message("damping")
  if (is.numeric(rate1) && rate1 > 0) ctr$rate1 = rate1 else message("rate1")
  if (is.numeric(rate2) && rate2 > 0) ctr$rate2 = rate2 else message("rate2")
  if (is.logical(verbose)) ctr$verbose = verbose else message("verbose")
  if (is.numeric(frequency) && frequency >= 1) ctr$frequency = floor(frequency) else message("frequency")
  if (is.logical(progress)) ctr$progress = progress else message("progress")
  if (is.logical(savedata)) ctr$savedata = savedata else message("verbose")

  # Additional consistency checks
  if (ctr$nafill > ctr$maxiter) {ctr$nafill = ctr$maxiter; message("nafill")}
  if (ctr$eps > 1e-01) {ctr$eps = 1e-01; message("eps")}
  if (ctr$rate1 > 1 - 1e-08) {ctr$rate1 = 1 - 1e-08; message("rate1")}
  if (ctr$rate2 > 1 - 1e-08) {ctr$rate2 = 1 - 1e-08; message("rate2")}
  if (ctr$frequency > ctr$maxiter) {ctr$frequency = ctr$maxiter; message("frequency")}

  # Return the check control parameters
  return (ctr)
}

#' @title Check and set the control parameters for the select optimization algorithm
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @param method optimization method to use
#' @param sampling sub-sampling method to use
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
#' \code{\link{set.control.block.sgd}} (\code{method="sgd"}, \code{sampling="block"}).
#' \code{\link{set.control.coord.sgd}} (\code{method="sgd"}, \code{sampling="coord"}),
#'
#' @returns A \code{list} of control parameters for the selected estimation algorithm
#'
#' @seealso \code{\link{set.control.init}}, \code{\link{set.control.cv}}, \code{\link{sgdgmf.fit}}
#'
#' @examples
#' library(sgdGMF)
#'
#' # Empty call
#' set.control.alg()
#'
#' # Parametrized call
#' set.control.alg(method = "airwls", control = list(maxiter = 200, stepsize = 0.3))
#'
#'
#' @export set.control.alg
set.control.alg = function (
    method = c("airwls", "newton", "sgd"),
    sampling = c("block", "coord", "rnd-block"),
    control = list()
) {
  # Check the optimization method
  method = match.arg(method)
  sampling = match.arg(sampling)

  # Set the input values for the control parameters, check whether the
  # they are allowed and, if they aren't, set them to default values
  if (method == "airwls") control = do.call("set.control.airwls", control)
  if (method == "newton") control = do.call("set.control.newton", control)
  if (method == "sgd" & sampling == "block") control = do.call("set.control.block.sgd", control)
  if (method == "sgd" & sampling == "coord") control = do.call("set.control.coord.sgd", control)
  if (method == "sgd" & sampling == "rnd-block") control = do.call("set.control.block.sgd", control)

  # Return the corrected control parameters
  return (control)
}


#' @title Check and set the cross-validation parameters
#'
#' @description
#' Check if the input cross-validation parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined cross-validation parameters.
#'
#' @param criterion information criterion to minimize for selecting the matrix rank
#' @param refit if \code{TRUE}, refit the model with the selected rank and return the fitted model
#' @param nfolds number of cross-validation folds
#' @param proportion proportion of the data to be used as test set in each fold
#' @param init initialization approach to use
#' @param verbose if \code{TRUE}, print the cross-validation status
#' @param parallel if \code{TRUE}, allows for parallel computing
#' @param nthreads number of cores to use in parallel (only if \code{parallel=TRUE})
#'
#' @returns A \code{list} of control parameters for the cross-validation algorithm
#'
#' @seealso \code{\link{set.control.init}}, \code{\link{set.control.alg}}, \code{\link{sgdgmf.cv}}
#'
#' @examples
#' library(sgdGMF)
#'
#' # Empty call
#' set.control.cv()
#'
#' # Parametrized call
#' set.control.cv(criterion = "bic", proportion = 0.2)
#'
#' @export set.control.cv
set.control.cv = function (
    criterion = c("dev", "mae", "mse", "aic", "bic"),
    refit = TRUE,
    nfolds = 5,
    proportion = 0.3,
    init = c("common", "separate"),
    verbose = FALSE,
    parallel = FALSE,
    nthreads = 1
) {
  # Set the default parameters
  ctr = list()
  ctr$criterion = match.arg(criterion)
  ctr$refit = TRUE
  ctr$nfolds = 5
  ctr$proportion = 0.3
  ctr$init = match.arg(init)
  ctr$verbose = FALSE
  ctr$parallel = FALSE
  ctr$nthreads = 1

  # Set the warning message
  message = function (var)
    warning(paste0("Cross-validation control: '", var,"' was set to default value."),
            call. = FALSE, immediate. = TRUE, domain = NULL)

  # Check all the input parameters
  if (is.logical(refit))
    ctr$refit = refit else message("refit")
  if (is.numeric(nfolds) && nfolds >= 1)
    ctr$nfolds = floor(nfolds) else message("nfolds")
  if (is.numeric(proportion) && proportion > 0 && proportion < 1)
    ctr$proportion = proportion else message("proportion")
  if (is.logical(verbose))
    ctr$verbose = verbose else message("verbose")
  if (is.logical(parallel))
    ctr$parallel = parallel else message("parallel")
  if (is.numeric(nthreads) && nthreads >= 1)
    ctr$nthreads = floor(nthreads) else message("nthreads")

  # Set the number of threads
  ncores = parallel::detectCores() - 1
  ctr$nthreads = floor(max(1, min(ctr$nthreads, ncores)))

  # Return the checked parameters
  return (ctr)
}
