
#' @title Summarizing a GMF model
#'
#' @description ...
#'
#' @param object an object of class \code{sgdgmf}
#'
#' @rdname summary
#' @export summary
summary = function (object) {
  UseMethod("summary.sgdgmf")
}

#' @title Extract the coefficient of a GMF model
#'
#' @description ...
#'
#' @param object an object of class \code{sgdgmf}
#'
#' @export
coefficients = function (object) {
  UseMethod("coef.sgdgmf")
}

#' @title Extract the residuals of a GMF model
#'
#' @description ...
#'
#' @param object an object of class \code{sgdgmf}
#' @param type the type of residuals which should be returned
#'
#' @export
residuals = function (
    object, type = c("deviance", "pearson", "working", "response")
) {
  UseMethod("residuals.sgdgmf")
}

#' @title Predict method for GMF fits
#'
#' @description ...
#'
#' @param object an object of class \code{sgdgmf}
#' @param newdata optionally, a list containing new values for \code{X} and \code{Z}
#' @param type the type of prediction which should be returned
#'
#' @export
predict = function (
    object, newdata = NULL, type = c("link", "response", "terms")
) {
  UseMethod("predict.sgdgmf")
}


#' @return \code{NULL}
#'
#' @rdname summary
#' @method summary sgdgmf
#' @export
summary.sgdgmf = function (object) {

}

#' @return \code{NULL}
#'
#' @rdname coefficients
#' @method coefficients sgdgmf
#' @export
coefficients.sgdgmf = function (
    object, type = c("all", "colreg", "rowreg", "scores", "loadings")
) {
  type = match.arg(type)
  switch(type,
    "all" = list(colcoef = object$B, rowcoef = object$A,
                 scores = object$U, loadings = object$V),
    "colreg" = object$B,
    "rowreg" = object$A,
    "scores" = object$U,
    "loadings" = object$V)
}


#' @return \code{NULL}
#'
#' @rdname residuals
#' @method residuals sgdgmf
#' @export
residuals.sgdgmf = function (
    object, type = c("deviance", "pearson", "working", "response")
) {

  type = match.arg(type)
  y = object$Y
  mu = object$mu
  var = object$var
  fam = object$family

  switch(type,
    "deviance" = sqrt(y - mu) * sqrt(fam$dev.resids(y, mu, 1)),
    "pearson" = (y - mu) / sqrt(fam$variance(mu)),
    "working" = (y - mu) / fam$mu.eta(mu),
    "response" = y - mu)
}

#' @return \code{NULL}
#'
#' @rdname predict
#' @method predict sgdgmf
#' @export
predict.sgdgmf = function (
    object, newdata = NULL, type = c("link", "response", "terms")
) {

  n = nrow(object$Y)
  m = ncol(object$Y)
  f = object$family
  type = match.arg(type)

  if (is.null(newdata)) {
    pred = switch(type,
      "link" = object$eta,
      "response" = object$mu,
      "terms" = list(
        XB = tcrossprod(object$X, object$B),
        AZ = tcrossprod(object$A, object$Z),
        UV = tcrossprod(object$U, object$V)))
  } else {
    if (is.null(newdata$X)) stop("The new 'X' is missing")
    if (is.null(newdata$Z)) stop("The new 'Z' is missing")
    if (!is.numeric(newdata$X)) stop("The new 'X' is not numeric.")
    if (!is.numeric(newdata$Z)) stop("The new 'Z' is not numeric.")
    if (!is.matrix(newdata$X)) stop("The new 'X' is not a matrix.")
    if (!is.matrix(newdata$Z)) stop("The new 'Z' is not a matrix.")
    if (dim(newdata$X) != c(n, p)) stop("The new 'X' has incompatible dimentions.")
    if (dim(newdata$Z) != c(m, q)) stop("The new 'Z' has incompatible dimentions.")

    pred = switch (object,
      "link" = tcrossprod(
        cbind(newdata$X, object$A, object$U),
        cbind(object$B, newdata$Z, object$V)),
      "response" = f$linkinv(tcrossprod(
        cbind(newdata$X, object$A, object$U),
        cbind(object$B, newdata$Z, object$V))),
      "terms" = list(
        XB = tcrossprod(newdata$X, object$B),
        AZ = tcrossprod(object$A, newdata$Z),
        UV = tcrossprod(object$U, object$V)))
  }

  # Output
  return (pred)
}



