
#' @keywords internal
qrrange = function (X, q = c(0.05, 0.95)) {
  range = quantile(c(X),q)
  X[X > range[2]] = range[2]
  X[X < range[1]] = range[1]
  X
}

#' @keywords internal
normx = function (x, p = 2) {
  n = 0
  if (p >= 0) {
    if (p == 0) n = max(x)
    if (p == 1) n = sum(abs(x))
    if (p == 2) n = sqrt(sum(x^2))
    if (!(p %in% c(0, 1, 2))) n = sum(x^p)^(1/p)
  } else {
    stop("p must be non-negative.")
  }
  return (n)
}

#' @keywords internal
max0 = function (x) {
  0.5 * (abs(x) + x)
}

#' @keywords internal
soft.threshold = function (x, a) {
  sign(x) * max0(abs(x) - a)
}


#' @title Rotate the matrices U and V
#'
#' @description
#' Rotate U and V in such a way that the transformed matrices
#' are such that U is orthogonal and V is lower triangular
#'
#' @importFrom whitening whiteningMatrix
#'
#' @keywords internal
correct.uv = function (U, V) {
  S = cov(U)

  if (ncol(U) == 1){
    return(list(u = U / sqrt(c(S)), v = V * sqrt(c(S))))
  }

  # Make cov of U identity
  W = whitening::whiteningMatrix(S)
  U = U %*% W
  V = V %*% t(solve(W))

  # Make V lower triangular
  V.qr = qr(t(V))
  U = U %*% qr.Q(V.qr)
  V = t(qr.R(V.qr))

  # Positive diagonal of V
  D = diag(V)
  V = t(sign(D) * t(V))
  U = t(sign(D) * t(U))

  list(u = U, v = V)
}

#' @title Procrustes distance
#' @description ...
#' @importFrom vegan procrustes
#' @keywords internal
norm.procrustes = function(A, B){
  A = A / norm(A, type = "F")
  B = B / norm(B, type = "F")
  vegan::procrustes(A, B)
}

#' @title Split the data matrix in train and test sets
#' @description ...
#' @keywords internal
partition = function (y, p = 0.3) {
  # Data dimensions
  n = nrow(y)
  m = ncol(y)
  s = floor(p*n*m)

  # Sparsification mask
  mask = cbind(ii = sample.int(n = n, size = s, replace = TRUE),
               jj = sample.int(n = m, size = s, replace = TRUE))

  # Train-test split
  train = test = y
  train[mask] = NA
  test[!is.na(train)] = NA

  # Return the splitted data
  list(train = train, test = test)
}
