
#' @title Extract the residuals of a GMF model
#'
#' @description ...
#'
#' @param object an object of class \code{sgdgmf}
#' @param type the type of residuals which should be returned
#'
#' @export
residuals.sgdgmf = function (
    object,
    type = c("deviance", "pearson", "working", "response")
) {

  type = match.arg(type)
  y = object$Y
  mu = object$mu
  var = object$var
  fam = object$family

  res = switch(type,
    deviance = sqrt(y - mu) * sqrt(fam$dev.resids(y, mu, 1)),
    pearson = (y - mu) / sqrt(fam$variance(mu)),
    working = (y - mu) / fam$mu.eta(mu),
    response = y - mu)

  return (res)
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
predict.sgdgmf = function (
    object,
    newdata = NULL,
    idx = NULL,
    type = c("link", "response", "terms")
) {

  n = nrow(object$Y)
  m = ncol(object$Y)
  f = object$family
  type = match.arg(type)

  if (is.null(newdata)) {
    pred = switch(type,
      "link" = object$eta,
      "response" = object$mu,
      "terms" = array(NA, dim = c(n, m, 3))
    )
    if (type == "terms") {
      pred[,,1] = tcrossprod(object$X, object$B)
      pred[,,2] = tcrossprod(object$A, object$Z)
      pred[,,3] = tcrossprod(object$U, object$V)
    }
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
      "terms" = array(NA, dim = c(n, m, 3))
    )

    if (type == "terms") {
      pred[,,1] = tcrossprod(newdata$X, object$B)
      pred[,,2] = tcrossprod(object$A, newdata$Z)
      pred[,,3] = tcrossprod(object$U, object$V)
    }
  }

  # Output
  return (pred)
}



