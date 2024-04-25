
#' @title Schwarz information criterion
#'
#' @description
#' Generic function calculating the Schwarz information criterion for one or several
#' fitted model objects for which a log-likelihood, or negative deviance, value can
#' be obtained, according to the formula \eqn{\text{SIC} = - 2 \,\text{log-likelihood} + p \log(n)/n},
#' where \eqn{p} represents the number effective degrees of freedom of the fitted model,
#' i.e. the effective number of parameters, and \eqn{n} is the number of observations.
#'
#' @param object a fitted model object
#' @param ... optional parameters
#'
#' @export
SIC = function (object, ...) {
  UseMethod("SIC", object)
}

