
#' @title Newton update for the parameters of a GMF model
#' @description Compute the Newton updates for U given V
#' @keywords internal
update.newton = function (U, V, penalty, ddiff, ddratio, stepsize, damping) {

  # The `sweep' command below, i.e.
  #   sweep(U, MARGIN = 2, STATS = penalty, FUN = "*")
  # is equivalent to `U %*% diag(penalty)'.
  # It multiply each column of `U' by the corresponding values of `penalty'
  # such that its output is a matrix with h-th column `penalty[h]*U[,h]'.

  # Compute the score matrix wrt U
  dU = - ddiff %*% V + sweep(U, 2, penalty, "*")

  # Compute the epementwise Hessian matrix wrt U
  ddU = ddratio %*% V**2 + sweep(U**0, 2, penalty, "*") + damping

  # Update the parameter matrix U
  U = matrix(U - stepsize * (dU / ddU), nrow = nrow(U), ncol = ncol(U))

  # Return the updated parameters
  return (U)
}
