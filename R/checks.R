
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
check.neg = function (object) {
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

#' @title Matrix dimension check
#'
#' @description
#' Check whether an object is a matrix with the expected dimensions
#'
#' @keywords internal
check.dim = function (object, n, m) {
  flag = FALSE
  # check if object exists
  if (!is.null(object)) {
    # check if object is a numeric matrix
    if (is.numeric(object) & is.matrix(object)) {
      # check is object is a matrix of appropriate dimensions
      if (nrow(object) == n & ncol(object) == m) {
        flag = TRUE
      }
    }
  }
  return (flag)
}

check.data = function (Y, X = NULL, Z = NULL) {
  # Check Y
  if (!is.numeric(Y)) stop("Y is not numeric.")
  if (!is.matrix(Y)) stop("Y is not a matrix.")
  n = nrow(Y)
  m = ncol(Y)

  # Check X
  if (!is.null(X)) {
    if (!is.numeric(X)) stop("X is not numeric.")
    if (!is.matrix(X)) stop("X is not a matrix.")
    if (nrow(X) != n) stop("The dimensions of X are not compatible with Y.")
    if (anyNA(X)) stop("X contains some NA.")
    if (sum(apply(X, 2, sd) == 0) > 1) stop("X has too many constant columns.")
  }

  # Check Z
  if (!is.null(Z)) {
    if (!is.numeric(Z)) stop("Z is not numeric.")
    if (!is.matrix(Z)) stop("Z is not a matrix.")
    if (nrow(Z) != m) stop("The dimensions of Z are not compatible with Y.")
    if (anyNA(Z)) stop("Z contains some NA.")
    if (sum(apply(Z, 2, sd) == 0) > 1) stop("Z has too many constant columns.")
  }
}
