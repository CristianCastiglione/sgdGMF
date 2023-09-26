
#' @title Check and set the penalty parameters
#'
#' @description
#' Check if the input penalty parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined penalty parameters.
#'
#' @keywords internal
set.penalty = function (penalty) {

  default = list(u = 1, v = 0, b = 0)

  check = function (object) {
    flag = FALSE
    if (!is.null(object)) {
      if (class(object) == "numeric") {
        if (object >= 0) {
          flag = TRUE
        }
      }
    }
    return (flag)
  }

  if (class(penalty) == "list") {
    if (check(penalty$u)) default$u = penalty$u
    if (check(penalty$v)) default$v = penalty$v
    if (check(penalty$b)) default$b = penalty$b
    if (default$u == 0 & default$v == 0) default$u = 1
  }

  return (default)
}

#' @title Check and set the initialization parameters
#'
#' @description
#' Check if the input initialization parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined initialization parameters.
#'
#' @keywords internal
set.init = function (init) {

  default = list()
  default$method = "ls-svd"
  default$values = list()
  default$niter = 5
  default$normalize = TRUE
  default$verbose = FALSE

  check = function (object, clss) {
    flag = FALSE
    if (!is.null(object)) {
      if (class(object) == clss) {
        flag = TRUE
      }
    }
    return (flag)
  }

  if (check(init$method, "character")) {
    if (init$method %in% c("glm-svd", "ls-svd", "svd", "random", "values")) {
      default$method = init$method
    }
  }
  if (check(init$niter, "numeric")) {
    if (floor(init$niter) >= 0) {
      default$niter = floor(init$niter)
    }
  }
  if (check(init$normalize, "logical")) {
    default$normalize = init$normalize
  }
  if (check(init$verbose, "logical")) {
    default$verbose = init$verbose
  }
  if (check(init$values, "list")) {
    default$values = init$values
  }

  return (default)
}

#' @title Check and set the control parameters for the optimization
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @keywords internal
set.control = function (method, control) {

  default = list()

  # Check whether an object exists and if it belongs to a particular class
  check.base = function (object, clss) {
    flag = FALSE
    if (!is.null(object)) {
      if (class(object) == clss) {
        flag = TRUE
      }
    }
    return (flag)
  }

  # Check whether an object is boolean
  check.bool = function (object) {
    flag = check.base(object, "logical")
    return (flag)
  }

  # Check whether an object is a positive number
  check.pos = function (object) {
    flag = FALSE
    if (check.base(object, "numeric")) {
      if (object > 0) {
        flag = TRUE
      }
    }
    return (flag)
  }

  # Check whether an object is a positive integer
  check.int = function (object) {
    flag = FALSE
    if (check.base(object, "numeric")) {
      if (all(floor(object) > 0)) {
        flag = TRUE
      }
    }
    return (flag)
  }

  if (!(method %in% c("airwls", "newton", "csgd", "asgd", "msgd", "svd"))) {
    stop("The specified optimization method is not implemented.")
  }

  # Set the input values for the control parameters, check whether the
  # they are allowed and, if they aren't, set them to default values
  if (method == "airwls") {
    default$normalize = TRUE
    default$maxiter = 100
    default$stepsize = 0.1
    default$nsteps = 1
    default$eps = 1e-04
    default$tol = 1e-03
    default$damping = 1e-08
    default$verbose = TRUE
    default$frequency = 10

    if (!is.null(control)) {
      # Standard safety checks
      if (check.bool(control$normalize)) default$normalize = control$normalize
      if (check.int (control$maxiter)) default$maxiter = floor(control$maxiter)
      if (check.pos (control$stepsize)) default$stepsize = control$stepsize
      if (check.int (control$nsteps)) default$nsteps = floor(control$nsteps)
      if (check.pos (control$eps)) default$eps = control$eps
      if (check.pos (control$tol)) default$tol = control$tol
      if (check.pos (control$damping)) default$damping = control$damping
      if (check.bool(control$verbose)) default$verbose = control$verbose
      if (check.int (control$frequency)) default$frequency = floor(control$frequency)

      # Addional consistency checks
      if (default$stepsize > 1) default$stepsize = 1
      if (default$eps > 1e-01) default$eps = 1e-01
      if (default$frequency > default$maxiter) default$frequency = default$maxiter
    }
  }
  if (method == "newton") {
    default$normalize = TRUE
    default$maxiter = 1000
    default$stepsize = 0.01
    default$eps = 1e-04
    default$tol = 1e-08
    default$damping = 1e-03
    default$verbose = TRUE
    default$frequency = 100

    if (!is.null(control)) {
      # Standard safety checks
      if (check.bool(control$normalize)) default$normalize = control$normalize
      if (check.int (control$maxiter)) default$maxiter = floor(control$maxiter)
      if (check.pos (control$stepsize)) default$stepsize = control$stepsize
      if (check.pos (control$eps)) default$eps = control$eps
      if (check.pos (control$tol)) default$tol = control$tol
      if (check.pos (control$damping)) default$damping = control$damping
      if (check.bool(control$verbose)) default$verbose = control$verbose
      if (check.int (control$frequency)) default$frequency = control$frequency

      # Addional consistency checks
      if (default$stepsize > 1) default$stepsize = 1
      if (default$eps > 1e-01) default$eps = 1e-01
      if (default$frequency > default$maxiter) default$frequency = default$maxiter
    }
  }
  if (method == "svd") {
    default$normalize = TRUE
    default$maxiter = 100
    default$stepsize = 0.1
    default$eps = 1e-04
    default$tol = 1e-03
    default$damping = 1e-08
    default$verbose = TRUE
    default$frequency = 10

    if (!is.null(control)) {
      # Standard safety checks
      if (check.bool(control$normalize)) default$normalize = control$normalize
      if (check.int (control$maxiter)) default$maxiter = floor(control$maxiter)
      if (check.pos (control$stepsize)) default$stepsize = control$stepsize
      if (check.pos (control$eps)) default$eps = control$eps
      if (check.pos (control$tol)) default$tol = control$tol
      if (check.pos (control$damping)) default$damping = control$damping
      if (check.bool(control$verbose)) default$verbose = control$verbose
      if (check.int (control$frequency)) default$frequency = floor(control$frequency)

      # Addional consistency checks
      if (default$stepsize > 1) default$stepsize = 1
      if (default$eps > 1e-01) default$eps = 1e-01
      if (default$frequency > default$maxiter) default$frequency = default$maxiter
    }
  }
  if (method == "csgd") {
    default$normalize = TRUE
    default$nafill = 10
    default$maxiter = 5000
    default$size = c(10, 10)
    default$eps = 1e-04
    default$rate0 = 0.01
    default$decay = 0.1
    default$damping = 1e-04
    default$rate1 = 0.05
    default$rate2 = 0.01
    default$burn = 2500
    default$verbose = TRUE
    default$frequency = 250

    if (!is.null(control)) {
      # Standard safety checks
      if (check.bool(control$normalize)) default$normalize = control$normalize
      if (check.int (control$nafill)) default$nafill = floor(control$nafill)
      if (check.int (control$maxiter)) default$maxiter = floor(control$maxiter)
      if (check.int (control$size)) default$size = floor(control$size)
      if (check.pos (control$eps)) default$eps = control$eps
      if (check.pos (control$rate0)) default$rate0 = control$rate0
      if (check.pos (control$decay)) default$decay = control$decay
      if (check.pos (control$damping)) default$damping = control$damping
      if (check.pos (control$rate1)) default$rate1 = control$rate1
      if (check.pos (control$rate2)) default$rate2 = control$rate2
      if (check.int (control$burn)) default$burn = floor(control$burn)
      if (check.bool(control$verbose)) default$verbose = control$verbose
      if (check.int (control$frequency)) default$frequency = floor(control$frequency)

      # Addional consistency checks
      if (default$nafill > default$maxiter) default$nafill = default$maxiter
      if (default$eps > 1e-01) default$eps = 1e-01
      if (default$rate1 > 1 - 1e-08) default$rate1 = 1 - 1e-08
      if (default$rate2 > 1 - 1e-08) default$rate1 = 1 - 1e-08
      if (default$burn > default$maxiter) default$burn = default$maxiter
      if (default$frequency > default$maxiter) default$frequency = default$maxiter
    }
  }
  if (method == "asgd") {
    default$normalize = TRUE
    default$nafill = 10
    default$maxiter = 5000
    default$eps = 1e-04
    default$size = c(10, 10)
    default$epochs = 10
    default$rate0 = 0.01
    default$decay = 0.1
    default$damping = 1e-04
    default$rate1 = 0.05
    default$rate2 = 0.01
    default$burn = 2500
    default$verbose = TRUE
    default$frequency = 250

    if (!is.null(control)) {
      # Standard safety checks
      if (check.bool(control$normalize)) default$normalize = control$normalize
      if (check.int (control$nafill)) default$nafill = floor(control$nafill)
      if (check.int (control$maxiter)) default$maxiter = floor(control$maxiter)
      if (check.int (control$size)) default$size = floor(control$size)
      if (check.pos (control$eps)) default$eps = control$eps
      if (check.pos (control$rate0)) default$rate0 = control$rate0
      if (check.pos (control$decay)) default$decay = control$decay
      if (check.pos (control$damping)) default$damping = control$damping
      if (check.pos (control$rate1)) default$rate1 = control$rate1
      if (check.pos (control$rate2)) default$rate2 = control$rate2
      if (check.int (control$burn)) default$burn = floor(control$burn)
      if (check.bool(control$verbose)) default$verbose = control$verbose
      if (check.int (control$frequency)) default$frequency = floor(control$frequency)

      # Addional consistency checks
      if (default$nafill > default$maxiter) default$nafill = default$maxiter
      if (default$eps > 1e-01) default$eps = 1e-01
      if (default$rate1 > 1 - 1e-08) default$rate1 = 1 - 1e-08
      if (default$rate2 > 1 - 1e-08) default$rate1 = 1 - 1e-08
      if (default$burn > default$maxiter) default$burn = default$maxiter
      if (default$frequency > default$maxiter) default$frequency = default$maxiter
    }
  }
  if (method == "msgd") {
    default$normalize = TRUE
    default$nafill = 10
    default$maxiter = 100
    default$size = 100
    default$eps = 1e-04
    default$rate0 = 0.01
    default$decay = 1.0
    default$damping = 1e-04
    default$rate1 = 0.05
    default$rate2 = 0.01
    default$burn = 90
    default$verbose = TRUE
    default$frequency = 10
    default$progress = FALSE

    if (!is.null(control)) {
      # Standard safety checks
      if (check.bool(control$normalize)) default$normalize = control$normalize
      if (check.int (control$nafill)) default$nafill = floor(control$nafill)
      if (check.int (control$maxiter)) default$maxiter = floor(control$maxiter)
      if (check.int (control$size)) default$size = floor(control$size)
      if (check.pos (control$eps)) default$eps = control$eps
      if (check.pos (control$rate0)) default$rate0 = control$rate0
      if (check.pos (control$decay)) default$decay = control$decay
      if (check.pos (control$damping)) default$damping = control$damping
      if (check.pos (control$rate1)) default$rate1 = control$rate1
      if (check.pos (control$rate2)) default$rate2 = control$rate2
      if (check.int (control$burn)) default$burn = floor(control$burn)
      if (check.bool(control$verbose)) default$verbose = control$verbose
      if (check.int (control$frequency)) default$frequency = floor(control$frequency)
      if (check.bool(control$progress)) default$progress = control$progress

      # Addional consistency checks
      if (default$nafill > default$maxiter) default$nafill = default$maxiter
      if (default$eps > 1e-01) default$eps = 1e-01
      if (default$rate1 > 1 - 1e-08) default$rate1 = 1 - 1e-08
      if (default$rate2 > 1 - 1e-08) default$rate1 = 1 - 1e-08
      if (default$burn > default$maxiter) default$burn = default$maxiter
      if (default$frequency > default$maxiter) default$frequency = default$maxiter
    }
  }

  return (default)
}

#' @title Set the jittering function for mapping the data
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
  if (substring(family$family, first = 1, last = 17) == "Negative Binomial") {
    f = function(x) family$linkfun(x + (x == 0) / 6)
  }
  return (f)
}
