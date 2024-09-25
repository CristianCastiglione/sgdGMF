
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
max.zero = function (x) {
  0.5 * (abs(x) + x)
}

#' @keywords internal
soft.threshold = function (x, a) {
  y = abs(x) - a
  0.5 * sign(x) * (abs(y) + y)
}

#' @title Procrustes rotation of two configurations
#' @description Rotates a configuration to maximum similarity with another configuration
#' @param X target matrix
#' @param Y matrix to be rotated
#' @param scale allow scaling of axes of Y
#' @param symmetric if \code{TRUE}, use symmetric Procrustes statistic
#' @keywords internal
procrustes <- function (X, Y, scale = TRUE, symmetric = FALSE) {
  # X <- vegan::scores(X, display = scores, ...)
  # Y <- vegan::scores(Y, display = scores, ...)
  if (nrow(X) != nrow(Y))
    stop(gettextf("matrices have different number of rows: %d and %d", nrow(X), nrow(Y)))
  if (ncol(X) < ncol(Y)) {
    warning("X has fewer axes than Y: X adjusted to comform Y\n")
    addcols <- ncol(Y) - ncol(X)
    for (i in 1:addcols) X <- cbind(X, 0)
  }
  ctrace <- function(mat) sum(mat^2)
  c <- 1
  if (symmetric) {
    X <- scale(X, scale = FALSE)
    Y <- scale(Y, scale = FALSE)
    X <- X / sqrt(ctrace(X))
    Y <- Y / sqrt(ctrace(Y))
  }
  xmean <- colMeans(X, 2, mean)
  ymean <- colMeans(Y, 2, mean)
  if (!symmetric) {
    X <- scale(X, scale = FALSE)
    Y <- scale(Y, scale = FALSE)
  }
  XY <- crossprod(X, Y)
  sol <- svd(XY)
  A <- sol$v %*% t(sol$u)
  if (scale) {
    c <- sum(sol$d) / ctrace(Y)
  }
  Yrot <- c * Y %*% A
  ## Translation (b) needs scale (c) although Mardia et al. do not
  ## have this. Reported by Christian Dudel.
  b <- xmean - c * ymean %*% A
  R2 <- ctrace(X) + c * c * ctrace(Y) - 2 * c * sum(sol$d)
  result <- list(Yrot = Yrot, X = X, ss = R2, rotation = A,
                 translation = b, scale = c, xmean = xmean,
                 symmetric = symmetric, call = match.call())
  result$svd <- sol
  class(reslt) <- "procrustes"
  return(result)
}

#' @title Procrustes distance
#' @description Compute the Procrustes distance between two matrices
#' @param A target matrix
#' @param B matrix to be rotated
#' @keywords internal
norm.procrustes = function(A, B){
  A = A / norm(A, type = "F")
  B = B / norm(B, type = "F")
  procrustes(A, B)
  # vegan::procrustes(A, B)
}

#' @title Fix sign ambiguity of eigen-vectors
#' @description Fix sign ambiguity of eigen-vectors by making U positive diagonal
#' @keywords internal
make.pos.diag = function(U) {
  sweep(U, 2, sign(diag(U)), "*")
}

#' @title Compute the whitening matrix from a given covariance matrix
#' @description Compute the whitening matrix from a given covariance matrix
#' @param sigma covariance matrix.
#' @param method determines the type of whitening transformation.
#' @keywords internal
whitening.matrix = function(sigma, method = c("ZCA", "ZCA-cor", "PCA", "PCA-cor", "Cholesky")) {

  W = NULL
  method = match.arg(method)
  switch(method,
    "ZCA" = {
      eS = eigen(sigma, symmetric = TRUE)
      U = eS$vectors
      lambda = eS$values
      W = U %*% diag(1 / sqrt(lambda)) %*% t(U)
    },
    "ZCA-cor" = {
      R = stats::cov2cor(sigma)
      eR = eigen(R, symmetric = TRUE)
      G = eR$vectors
      theta = eR$values
      W = G %*% diag(1 / sqrt(theta)) %*% t(G) %*% diag(1 / sqrt(v))
    },
    "PCA" = {
      eS = eigen(sigma, symmetric = TRUE)
      U = eS$vectors
      lambda = eS$values
      U = make.pos.diag(U)
      W = diag(1 / sqrt(lambda)) %*% t(U)
    },
    "PCA-cor" = {
      v = diag(sigma)
      R = stats::cov2cor(sigma)
      eR = eigen(R, symmetric = TRUE)
      G = eR$vectors
      theta = eR$values
      G = make.pos.diag(G)
      W = diag(1 / sqrt(theta)) %*% t(G) %*% diag(1 / sqrt(v))
    },
    "Cholesky" = {
      W = solve(t(chol(sigma)))
    })


  # if (method == "ZCA" | method == "PCA") {
  #   eS = eigen(sigma, symmetric = TRUE)
  #   U = eS$vectors
  #   lambda = eS$values
  # }
  # if (method == "ZCA-cor" | method == "PCA-cor") {
  #   R = stats::cov2cor(sigma)
  #   eR = eigen(R, symmetric = TRUE)
  #   G = eR$vectors
  #   theta = eR$values
  # }
  # if (method == "ZCA") {
  #   W = U %*% diag(1 / sqrt(lambda)) %*% t(U)
  # }
  # if (method == "PCA") {
  #   U = make.pos.diag(U)
  #   W = diag(1 / sqrt(lambda)) %*% t(U)
  # }
  # if (method == "Cholesky") {
  #   W = solve(t(chol(sigma)))
  # }
  # if (method == "ZCA-cor"){
  #   W = G %*% diag(1 / sqrt(theta)) %*% t(G) %*% diag(1 / sqrt(v))
  # }
  # if (method == "PCA-cor") {
  #   G = make.pos.diag(G)
  #   W = diag(1 / sqrt(theta)) %*% t(G) %*% diag(1 / sqrt(v))
  # }

  return(W)
}

#' @title Normalize the matrices U and V
#'
#' @description
#' Rotate U and V using either QR or SVD decompositions.
#' The QR methods rotate U and V in such a way to obtain an orthogonal U
#' and a lower triangular V.  The SVD method rotate U and V in such a way
#' to obtain an orthogonal U and a scaled orthogonal V.
#'
#' @keywords internal
normalize.uv = function (U, V, method = c("qr", "svd")) {

  method = match.arg(method)

  if (method == "svd") {
    try({
      ncomp = ncol(U)
      uv = tcrossprod(U, V)
      s = RSpectra::svds(uv, ncomp)
      if (ncomp == 1) {
        U = s$u
        V = s$v *s$d
      } else {
        U = s$u
        V = s$v %*% diag(s$d)
      }
    })
  }
  if (method == "qr") {
    if (is.vector(U)){
      S = sd(c(U))
      U = U / sqrt(S)
      V = V * sqrt(S)
    }
    if (is.matrix(U) & ncol(U) == 1){
      S = sd(c(U))
      U = U / sqrt(S)
      V = V * sqrt(S)
    }
    if (is.matrix(U) & ncol(U) > 1) {
      try({
        # Compute the covariance of U
        S = cov(U)

        # Make the covariance of U identity
        W = whitening.matrix(S)
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
      })
    }
  }

  # Output
  list(U = U, V = V)
}

#' @title Split the data matrix in train and test sets
#'
#' @description
#' Returns a list of two matrices \code{train} and \code{test}.
#' \code{train} corresponds to the input matrix with a fixed persentage of
#' entries masked by NA values. \code{test} is the complement of \code{train}
#' and contains the values of the input matrix in the cells where \code{train}
#' is NA, while all the other entries are filled by NA values.
#'
#' @param y input matrix to be split into train and test sets
#' @param p fraction of observations to be used for the test set
#'
#' @keywords internal
partition = function (y, p = 0.3) {
  # Data dimensions
  n = nrow(y)
  m = ncol(y)
  s = floor(p*n*m)

  # Sparsification mask
  mask = cbind(row = sample.int(n = n, size = s, replace = TRUE),
               col = sample.int(n = m, size = s, replace = TRUE))

  # Train-test split
  train = test = y
  train[mask] = NA
  test[!is.na(train)] = NA

  # Return the splitted data
  list(train = train, test = test)
}
