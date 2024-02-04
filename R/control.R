
#' @title Class check
#' @description Check whether an object exists and if it belongs to a particular class
#' @keywords internal
check.class = function (object, clss) {
  flag = FALSE
  if (!is.null(object)) {
    if (class(object) == clss) {
      flag = TRUE
    }
  }
  return (flag)
}

#' @title Boolean check
#' @description Check whether an object is boolean
#' @keywords internal
check.bool = function (object) {
  flag = check.class(object, "logical")
  return (flag)
}

#' @title Positive real check
#' @description Check whether an object is a positive number
#' @keywords internal
check.pos = function (object) {
  flag = FALSE
  if (check.class(object, "numeric")) {
    if (all(object > 0)) {
      flag = TRUE
    }
  }
  return (flag)
}

#' @title Non-negative real check
#' @description Check whether an object is a positive number
#' @keywords internal
check.nneg = function (object) {
  flag = FALSE
  if (check.class(object, "numeric")) {
    if (all(object >= 0)) {
      flag = TRUE
    }
  }
  return (flag)
}

#' @title Integer check
#'
#' @description
#' Check whether an object is a positive integer
#'
#' @keywords internal
check.int = function (object) {
  flag = FALSE
  if (check.class(object, "numeric")) {
    if (all(floor(object) > 0)) {
      flag = TRUE
    }
  }
  return (flag)
}


#' @title Check and set the penalty parameters
#'
#' @description
#' Check if the input penalty parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined penalty parameters.
#'
#' @keywords internal
set.penalty = function (penalty) {

  default = list(u = 1, v = 0, b = 0, a = 0)

  if (check.class(penalty, "list")) {
    if (check.nneg(penalty$u)) default$u = penalty$u
    if (check.nneg(penalty$v)) default$v = penalty$v
    if (check.nneg(penalty$b)) default$b = penalty$b
    if (check.nneg(penalty$a)) default$a = penalty$a
    if (default$u == 0 & default$v == 0) default$u = 1
  }

  penalty = c(default$b, default$a, default$u, default$v)

  return (penalty)
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
  default$method = "svd"
  default$values = list()
  default$niter = 5
  default$normalize = TRUE
  default$verbose = FALSE

  if (check(init$method, "character")) {
    if (init$method %in% c("glm", "svd", "random", "values")) {
      default$method = init$method
    }
  }
  if (check.class(init$niter, "numeric")) {
    if (floor(init$niter) >= 0) {
      default$niter = floor(init$niter)
    }
  }
  if (check.class(init$normalize, "logical")) {
    default$normalize = init$normalize
  }
  if (check.class(init$verbose, "logical")) {
    default$verbose = init$verbose
  }
  if (check.class(init$values, "list")) {
    default$values = init$values
  }

  return (default)
}

#' @title Check and set the control parameters for the AIRWLS algorithm
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @keywords internal
set.control.airwls = function (control) {
  # Set the default control parameters
  default = list()
  default$normalize = TRUE
  default$maxiter = 100
  default$nstep = 1
  default$stepsize = 0.1
  default$eps = 1e-04
  default$nafill = 1
  default$tol = 1e-05
  default$damping = 1e-03
  default$verbose = TRUE
  default$frequency = 10
  default$parallel = FALSE
  default$nthreads = 1

  if (!is.null(control)) {
    # Standard safety checks
    if (check.bool(control$normalize)) default$normalize = control$normalize
    if (check.int(control$maxiter)) default$maxiter = floor(control$maxiter)
    if (check.int(control$nsteps)) default$nsteps = floor(control$nsteps)
    if (check.pos(control$stepsize)) default$stepsize = control$stepsize
    if (check.pos(control$eps)) default$eps = control$eps
    if (check.int(control$nafill)) default$nafill = floor(control$nafill)
    if (check.pos(control$tol)) default$tol = control$tol
    if (check.pos(control$damping)) default$damping = control$damping
    if (check.bool(control$verbose)) default$verbose = control$verbose
    if (check.int(control$frequency)) default$frequency = floor(control$frequency)
    if (check.bool(control$parallel)) default$parallel = control$parallel
    if (check.int(control$nthreads)) default$nthreads = floor(control$nthreads)

    # Addional consistency checks
    if (default$stepsize > 1) default$stepsize = 1
    if (default$eps > 1e-01) default$eps = 1e-01
    if (default$frequency > default$maxiter) default$frequency = default$maxiter
  }

  return (default)
}

#' @title Check and set the control parameters for the Newton algorithm
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @keywords internal
set.control.newton = function (control) {
  # Set the default control parameters
  default = list()
  default$normalize = TRUE
  default$maxiter = 500
  default$stepsize = 0.01
  default$eps = 1e-04
  default$nafill = 1
  default$tol = 1e-05
  default$damping = 1e-03
  default$verbose = TRUE
  default$frequency = 50
  default$parallel = FALSE
  default$nthreads = 1

  if (!is.null(control)) {
    # Standard safety checks
    if (check.bool(control$normalize)) default$normalize = control$normalize
    if (check.int(control$maxiter)) default$maxiter = floor(control$maxiter)
    if (check.pos(control$stepsize)) default$stepsize = control$stepsize
    if (check.pos(control$eps)) default$eps = control$eps
    if (check.int(control$nafill)) default$nafill = floor(control$nafill)
    if (check.pos(control$tol)) default$tol = control$tol
    if (check.pos(control$damping)) default$damping = control$damping
    if (check.bool(control$verbose)) default$verbose = control$verbose
    if (check.int(control$frequency)) default$frequency = control$frequency
    if (check.bool(control$parallel)) default$parallel = control$parallel
    if (check.int(control$nthreads)) default$nthreads = floor(control$nthreads)

    # Addional consistency checks
    if (default$stepsize > 1) default$stepsize = 1
    if (default$eps > 1e-01) default$eps = 1e-01
    if (default$frequency > default$maxiter) default$frequency = default$maxiter
  }

  return (default)
}

#' @title Check and set the control parameters for the memoized-SGD algorithm
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @keywords internal
set.control.msgd = function (control) {
  # Set the default control parameters
  defualt = list()
  default$normalize = TRUE
  default$maxiter = 100
  default$eps = 1e-04
  default$nafill = 10
  default$tol = 1e-05
  default$size = 100
  default$burn = 90
  default$rate0 = 0.01
  default$decay = 1.0
  default$damping = 1e-04
  default$rate = c(0.05, 0.01)
  default$parallel = FALSE
  default$nthreads = 1
  default$verbose = TRUE
  default$frequency = 10
  default$progress = FALSE

  if (!is.null(control)) {
    # Standard safety checks
    if (check.bool(control$normalize)) default$normalize = control$normalize
    if (check.int(control$maxiter)) default$maxiter = floor(control$maxiter)
    if (check.pos(control$eps)) default$eps = control$eps
    if (check.int(control$nafill)) default$nafill = floor(control$nafill)
    if (check.pos(control$tol)) default$tol = control$tol
    if (check.int(control$size)) default$size = floor(control$size)
    if (check.int(control$burn)) default$burn = floor(control$burn)
    if (check.pos(control$rate0)) default$rate0 = control$rate0
    if (check.pos(control$decay)) default$decay = control$decay
    if (check.pos(control$damping)) default$damping = control$damping
    if (check.pos(control$rate)) default$rate = control$rate[1:2]
    if (check.bool(control$parallel)) default$parallel = control$parallel
    if (check.int(control$nthreads)) default$nthreads = floor(control$nthreads)
    if (check.bool(control$verbose)) default$verbose = control$verbose
    if (check.int(control$frequency)) default$frequency = floor(control$frequency)
    if (check.bool(control$progress)) default$progress = control$progress

    # Addional consistency checks
    if (default$nafill > default$maxiter) default$nafill = default$maxiter
    if (default$eps > 1e-01) default$eps = 1e-01
    if (default$rate[1] > 1 - 1e-08) default$rate[1] = 1 - 1e-08
    if (default$rate[2] > 1 - 1e-08) default$rate[2] = 1 - 1e-08
    if (default$burn > default$maxiter) default$burn = default$maxiter
    if (default$frequency > default$maxiter) default$frequency = default$maxiter
  }

  return (default)
}

#' @title Check and set the control parameters for the coordinate-SGD algorithm
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @keywords internal
set.control.csgd = function (control) {
  # Set the default control parameters
  default = list()
  default$normalize = TRUE
  default$maxiter = 1000
  default$eps = 0.01
  default$nafill = 10
  default$tol = 1e-08
  default$size = c(100, 100)
  default$burn = 0.75
  default$rate0 = 0.01
  default$decay = 0.01
  default$damping = 1e-03
  default$rate = c(0.95, 0.99)
  default$parallel = FALSE
  default$nthreads = 1
  default$verbose = TRUE
  default$frequency = 250
  default$progress = FALSE

  if (!is.null(control)) {
    # Standard safety checks
    if (check.bool(control$normalize)) default$normalize = control$normalize
    if (check.int(control$maxiter)) default$maxiter = floor(control$maxiter)
    if (check.pos(control$eps)) default$eps = control$eps
    if (check.int(control$nafill)) default$nafill = floor(control$nafill)
    if (check.pos(control$tol)) default$tol = control$tol
    if (check.int(control$size)) default$size = floor(control$size[1:2])
    if (check.int(control$burn)) default$burn = floor(control$burn)
    if (check.pos(control$rate0)) default$rate0 = control$rate0
    if (check.pos(control$decay)) default$decay = control$decay
    if (check.pos(control$damping)) default$damping = control$damping
    if (check.pos(control$rate)) default$rate = control$rate[1:2]
    if (check.bool(control$parallel)) default$parallel = control$parallel
    if (check.int(control$nthreads)) default$nthreads = floor(control$nthreads)
    if (check.bool(control$verbose)) default$verbose = control$verbose
    if (check.int(control$frequency)) default$frequency = floor(control$frequency)
    if (check.bool(control$progress)) default$progress = control$progress

    # Addional consistency checks
    if (default$nafill > default$maxiter) default$nafill = default$maxiter
    if (default$eps > 1e-01) default$eps = 1e-01
    if (default$rate[1] > 1 - 1e-08) default$rate[1] = 1 - 1e-08
    if (default$rate[2] > 1 - 1e-08) default$rate[2] = 1 - 1e-08
    if (default$burn > default$maxiter) default$burn = default$maxiter
    if (default$frequency > default$maxiter) default$frequency = default$maxiter
  }

  return (default)
}

#' @title Check and set the control parameters for the rowwise-SGD algorithm
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @keywords internal
set.control.rsgd = function (control) {
  # Set the default control parameters
  defualt = list()
  default$normalize = TRUE
  default$maxiter = 100
  default$eps = 1e-04
  default$nafill = 10
  default$tol = 1e-05
  default$size = 100
  default$burn = 90
  default$rate0 = 0.01
  default$decay = 1.0
  default$damping = 1e-04
  default$rate = c(0.05, 0.01)
  default$parallel = FALSE
  default$nthreads = 1
  default$verbose = TRUE
  default$frequency = 10
  default$progress = FALSE

  if (!is.null(control)) {
    # Standard safety checks
    if (check.bool(control$normalize)) default$normalize = control$normalize
    if (check.int(control$maxiter)) default$maxiter = floor(control$maxiter)
    if (check.pos(control$eps)) default$eps = control$eps
    if (check.int(control$nafill)) default$nafill = floor(control$nafill)
    if (check.pos(control$tol)) default$tol = control$tol
    if (check.int(control$size)) default$size = floor(control$size)
    if (check.int(control$burn)) default$burn = floor(control$burn)
    if (check.pos(control$rate0)) default$rate0 = control$rate0
    if (check.pos(control$decay)) default$decay = control$decay
    if (check.pos(control$damping)) default$damping = control$damping
    if (check.pos(control$rate)) default$rate = control$rate[1:2]
    if (check.bool(control$parallel)) default$parallel = control$parallel
    if (check.int(control$nthreads)) default$nthreads = floor(control$nthreads)
    if (check.bool(control$verbose)) default$verbose = control$verbose
    if (check.int(control$frequency)) default$frequency = floor(control$frequency)
    if (check.bool(control$progress)) default$progress = control$progress

    # Addional consistency checks
    if (default$nafill > default$maxiter) default$nafill = default$maxiter
    if (default$eps > 1e-01) default$eps = 1e-01
    if (default$rate[1] > 1 - 1e-08) default$rate[1] = 1 - 1e-08
    if (default$rate[2] > 1 - 1e-08) default$rate[2] = 1 - 1e-08
    if (default$burn > default$maxiter) default$burn = default$maxiter
    if (default$frequency > default$maxiter) default$frequency = default$maxiter
  }

  return (default)
}

#' @title Check and set the control parameters for the blockwise-SGD algorithm
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @keywords internal
set.control.bsgd = function (control) {
  # Set the default control parameters
  default = list()
  default$normalize = TRUE
  default$maxiter = 1000
  default$eps = 0.01
  default$nafill = 10
  default$tol = 1e-08
  default$size = c(100, 100)
  default$burn = 0.75
  default$rate0 = 0.01
  default$decay = 0.01
  default$damping = 1e-03
  default$rate = c(0.95, 0.99)
  default$parallel = FALSE
  default$nthreads = 1
  default$verbose = TRUE
  default$frequency = 250
  default$progress = FALSE

  if (!is.null(control)) {
    # Standard safety checks
    if (check.bool(control$normalize)) default$normalize = control$normalize
    if (check.int(control$maxiter)) default$maxiter = floor(control$maxiter)
    if (check.pos(control$eps)) default$eps = control$eps
    if (check.int(control$nafill)) default$nafill = floor(control$nafill)
    if (check.pos(control$tol)) default$tol = control$tol
    if (check.int(control$size)) default$size = floor(control$size[1:2])
    if (check.int(control$burn)) default$burn = floor(control$burn)
    if (check.pos(control$rate0)) default$rate0 = control$rate0
    if (check.pos(control$decay)) default$decay = control$decay
    if (check.pos(control$damping)) default$damping = control$damping
    if (check.pos(control$rate)) default$rate = control$rate[1:2]
    if (check.bool(control$parallel)) default$parallel = control$parallel
    if (check.int(control$nthreads)) default$nthreads = floor(control$nthreads)
    if (check.bool(control$verbose)) default$verbose = control$verbose
    if (check.int(control$frequency)) default$frequency = floor(control$frequency)
    if (check.bool(control$progress)) default$progress = control$progress

    # Addional consistency checks
    if (default$nafill > default$maxiter) default$nafill = default$maxiter
    if (default$eps > 1e-01) default$eps = 1e-01
    if (default$rate[1] > 1 - 1e-08) default$rate[1] = 1 - 1e-08
    if (default$rate[2] > 1 - 1e-08) default$rate[2] = 1 - 1e-08
    if (default$burn > default$maxiter) default$burn = default$maxiter
    if (default$frequency > default$maxiter) default$frequency = default$maxiter
  }

  return (default)
}

#' @title Check and set the control parameters for the selecte optimization algorithm
#'
#' @description
#' Check if the input control parameters are allowed and set them to default
#' values if they are not. Returns a list of well-defined control parameters.
#'
#' @keywords internal
set.control = function (method, control) {

  if (!(method %in% c("airwls", "newton", "msgd", "csgd", "rsgd", "bsgd"))) {
    stop("The specified optimization method is not implemented.")
  }

  # Set the input values for the control parameters, check whether the
  # they are allowed and, if they aren't, set them to default values
  if (method == "airwls") control = set.control.airwls(control)
  if (method == "newton") control = set.control.newton(control)
  if (method == "msgd") control = set.control.msgd(control)
  if (method == "csgd") control = set.control.csgd(control)
  if (method == "rsgd") control = set.control.rsgd(control)
  if (method == "bsgd") control = set.control.bsgd(control)

  return (control)
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
