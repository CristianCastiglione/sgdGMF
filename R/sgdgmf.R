

#' @title Factorize a matrix of non-Gaussian observations
#'
#' @description
#' ...
#'
#' @param Y matrix of responses (\eqn{n \times m})
#' @param X matrix of row fixed effects (\eqn{n \times p})
#' @param Z matrix of column fixed effects (\eqn{q \times m})
#' @param family a family as in the \code{\link{glm}} interface
#' @param ncomp number of random effects to estimate (default 2)
#' @param intercept should the model include a row/col intercept? (default FALSE)
#' @param dispersion should the method estimate column-specific disperssion parameters? (default FALSE)
#' @param penalty list of penalty parameters corresponding to u, v, a and b
#' @param init list of initialization options
#' @param control list of optimization options
#'
#' @return
#' ...
#'
#' @details
#' ...
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
sgdgmf = function (
    Y,
    X = NULL,
    Z = NULL,
    family = poisson(),
    ncomp = 2,
    method = "newton",
    intercept = c(FALSE, FALSE),
    dispersion = FALSE,
    penalty = list(),
    init = list(),
    control = list()) {

  # Chose the optimizer
  gmf.fit = NULL
  methods = c("airwls", "newton", "m-sgd", "b-sgd", "c-sgd")
  if (method %in% methods) {
    if (method == "airwls") gmf.fit = gmf.airwls
    if (method == "newton") gmf.fit = gmf.newton
    if (method == "m-sgd") gmf.fit = gmf.memo.sgd
    if (method == "b-sgd") gmf.fit = gmf.block.sgd
    if (method == "c-sgd") gmf.fit = gmf.coord.sgd
  } else {
    stop("Method = `", method,"' is not implemented yet.")
  }

  # Sanity check of the intercept parameter
  if (is.logical(intercept)) {
    if (length(intercept) < 2) intercept = rep(intercept, 2)
    if (length(intercept) > 2) intercept = intercept[1:2]
  } else {
    stop("Intercept must be a a logical vector of two components.")
  }

  # Fit the model using the chosen optimizer
  fit = gmf.fit(
    Y = Y, X = X, Z = Z, family = family, ncomp = ncomp,
    intercept = intercept, dispersion = dispersion,
    penalty = penalty, init = init, control = control)

  # Convert the result in a S3 class
  class(fit) = "sgdgmf"

  # Return the fitted model
  return (fit)
}

