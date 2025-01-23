
#' @title Pointwise deviance of a GMF model
#' @description Compute the pointwise deviance for all the observations in the sample
#' @keywords internal
pointwise.deviance = function (mu, y, family = gaussian()) {
  if (length(mu) == 1) {
    mut = y
    mut[] = mu
    mu = mut
  }
  nona = !is.na(y)
  dev = y
  dev[] = NA
  dev[nona] = family$dev.resids(y[nona], mu[nona], 1)
  return(dev)
}

#' @title Model deviance of a GMF model
#' @description Compute the overall deviance averaging the contributions of all data
#' @keywords internal
matrix.deviance = function (mu, y, family = gaussian()) {
  dev = pointwise.deviance(mu, y, family)
  dev = sum(dev, na.rm = TRUE)
  # dev = mean(dev, na.rm = TRUE)
  return (dev)
}

#' @title Frobenius penalty for the parameters of a GMF model
#' @description Compute the Frobenius penalty for all the parameters in the model
#' @keywords internal
matrix.penalty = function (U, penalty) {
  pen = sum(sweep(U**2, 2, penalty, "*"))
  return (pen)
}

