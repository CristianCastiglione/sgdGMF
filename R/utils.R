



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
#' @param U target matrix
#' @keywords internal
make.pos.diag = function(U) {
  sweep(U, 2, sign(diag(U)), "*")
}

#' @title Compute the whitening matrix from a given covariance matrix
#' @description Compute the whitening matrix from a given covariance matrix
#' @param sigma covariance matrix.
#' @param method determines the type of whitening transformation.
#' @details
#' This function is an internal re-implementation of the function \code{whiteningMatrix}
#' in the \code{whitening} package. See the original documentation to get more details.
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
      v = diag(sigma)
      R = stats::cov2cor(sigma)
      eR = eigen(R, symmetric = TRUE)
      G = eR$vectors
      theta = eR$values
      W = G %*% diag(1 / sqrt(theta)) %*% t(G) %*% diag(1 / sqrt(v))
    },
    "PCA" = {
      eS = eigen(sigma, symmetric = TRUE)
      U = make.pos.diag(eS$vectors)
      lambda = eS$values
      W = diag(1 / sqrt(lambda)) %*% t(U)
    },
    "PCA-cor" = {
      v = diag(sigma)
      R = stats::cov2cor(sigma)
      eR = eigen(R, symmetric = TRUE)
      G = make.pos.diag(eR$vectors)
      theta = eR$values
      W = diag(1 / sqrt(theta)) %*% t(G) %*% diag(1 / sqrt(v))
    },
    "Cholesky" = {
      W = solve(t(chol(sigma)))
    })

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
      S = stats::sd(c(U))
      U = U / sqrt(S)
      V = V * sqrt(S)
    }
    if (is.matrix(U) & ncol(U) == 1){
      S = stats::sd(c(U))
      U = U / sqrt(S)
      V = V * sqrt(S)
    }
    if (is.matrix(U) & ncol(U) > 1) {
      try({
        # Compute the covariance of U
        S = stats::cov(U)

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

#' @title Simulate non-Gaussian data from a GMF model
#'
#' @description
#' Simulate synthetic non-Gaussian data from a generalized matrix factorization (GMF) model.
#'
#' @param n number of observations
#' @param m number of variables
#' @param ncomp rank of the latent matrix factorization
#' @param family a \code{glm} family (see \code{\link{family}} for more details)
#' @param dispersion a positive dispersion parameter
#'
#' @return
#' A list containing the following objects:
#' \itemize{
#'   \item \code{Y}: simulated response matrix
#'   \item \code{U}: simulated factor matrix
#'   \item \code{V}: simulated loading matrix
#'   \item \code{eta}: linear predictor matrix
#'   \item \code{mu}: conditional mean matrix
#'   \item \code{phi}: scalar dispersion parameter
#'   \item \code{family}: model family
#'   \item \code{ncomp}: rank of the latent matrix factorization
#'   \item \code{param}: a list containing time, phase, frequency and amplitude vectors used to generate \code{U}
#' }
#'
#' @details
#' The loadings, \code{V}, are independently sampled from a standard normal distribution.
#' The scores, \code{U}, are simulated according to sinusoidal signals evaluated at different
#' phases, frequencies and amplitudes. These parameters are randomly sampled from independent
#' uniform distributions.
#'
#' @examples
#' library(sgdGMF)
#'
#' # Set the data dimensions
#' n = 100; m = 20; d = 5
#'
#' # Generate data using Poisson, Binomial and Gamma models
#' data_pois = sim.gmf.data(n = n, m = m, ncomp = d, family = poisson())
#' data_bin = sim.gmf.data(n = n, m = m, ncomp = d, family = binomial())
#' data_gam = sim.gmf.data(n = n, m = m, ncomp = d, family = Gamma(link = "log"), dispersion = 0.25)
#'
#' # Compare the results
#' par(mfrow = c(3,3), mar = c(1,1,3,1))
#' image(data_pois$Y, axes = FALSE, main = expression(Y[Pois]))
#' image(data_pois$mu, axes = FALSE, main = expression(mu[Pois]))
#' image(data_pois$U, axes = FALSE, main = expression(U[Pois]))
#' image(data_bin$Y, axes = FALSE, main = expression(Y[Bin]))
#' image(data_bin$mu, axes = FALSE, main = expression(mu[Bin]))
#' image(data_bin$U, axes = FALSE, main = expression(U[Bin]))
#' image(data_gam$Y, axes = FALSE, main = expression(Y[Gam]))
#' image(data_gam$mu, axes = FALSE, main = expression(mu[Gam]))
#' image(data_gam$U, axes = FALSE, main = expression(U[Gam]))
#'
#' @export sim.gmf.data
sim.gmf.data = function (n = 100, m = 20, ncomp = 5, family = gaussian(), dispersion = 1) {

  # Set the time range, phases, frequences and amplitudes of the underlying signals
  time = seq(from = 0, to = 1, length = n)
  phase = sort(stats::runif(ncomp), decreasing = FALSE) * 2
  freq = sort(stats::runif(ncomp), decreasing = FALSE) * 2
  amp = sort(stats::runif(ncomp), decreasing = TRUE)

  # Combine the latent signals using independent Gaussian coefficients
  V = matrix(stats::rnorm(m * ncomp), nrow = m, ncol = ncomp)
  U = do.call("cbind", lapply(1:ncomp, function(h){
    amp[h] * sin(2 * pi * freq[h] * (time + phase[h]))
  }))

  # Compute the linear predictor, the conditional mean and the dispersion parameter
  eta = tcrossprod(U, V)
  mu = family$linkinv(eta)
  phi = 3*rbeta(1, 2, 3)
  phi = ifelse(!is.null(dispersion), abs(dispersion), 3 * stats::rbeta(1, 2, 3))

  # Simulate the data using an dispersion exponential family distribution
  Y = matrix(
    switch(family$family,
      "gaussian" = stats::rnorm(n * m, mean = mu, sd = sqrt(phi)),
      "binomial" = stats::rbinom(n * m, size = 1, prob = mu),
      "poisson" = stats::rpois(n * m, lambda = mu),
      "Gamma" = stats::rgamma(n * m, shape = 1 / phi, scale = phi * mu),
      "negbinom" = MASS::rnegbin(n * m, mu = mu, theta = phi)),
    nrow = n, ncol = m)

  # Return the simulated signals
  list(Y = Y, U = U, V = V,
       eta = eta, mu = mu, phi = phi,
       ncomp = ncomp, family = family)
}

