
#' @title Procrustes rotation of two configurations
#' @description Rotates a configuration to maximum similarity with another configuration
#' @param X target matrix
#' @param Y matrix to be rotated
#' @param scale allow scaling of axes of Y
#' @param symmetric if \code{TRUE}, use symmetric Procrustes statistic
#' @keywords internal
procrustes <- function (X, Y, scale = TRUE, symmetric = FALSE) {

  # Safety checks
  if (nrow(X) != nrow(Y))
    stop(gettextf("matrices have different number of rows: %d and %d", nrow(X), nrow(Y)))
  if (ncol(X) < ncol(Y)) {
    warning("X has fewer axes than Y: X adjusted to comform Y\n")
    addcols <- ncol(Y) - ncol(X)
    for (i in 1:addcols) X <- cbind(X, 0)
  }

  # Squared matrix trace (Frobenius norm)
  ctrace <- function(mat) sum(mat^2)

  # Symmetric standardization
  if (symmetric) {
    X <- scale(X, scale = FALSE)
    Y <- scale(Y, scale = FALSE)
    X <- X / sqrt(ctrace(X))
    Y <- Y / sqrt(ctrace(Y))
  }

  # Get the colum-wise means
  xmean <- colMeans(X, 2, mean)
  ymean <- colMeans(Y, 2, mean)

  # Asymmetric standardization
  if (!symmetric) {
    X <- scale(X, scale = FALSE)
    Y <- scale(Y, scale = FALSE)
  }

  # Get the rotation matrix
  sol <- svd(crossprod(X, Y))
  A <- sol$v %*% t(sol$u)
  c <- 1
  if (scale) {
    c <- sum(sol$d) / ctrace(Y)
  }

  # Rotate and scale the matrix
  Yrot <- c * Y %*% A

  # Transformed mean
  b <- xmean - c * ymean %*% A

  # Residual sum of squares (R2)
  R2 <- ctrace(X) + c * c * ctrace(Y) - 2 * c * sum(sol$d)

  # Output object
  result <- list(Yrot = Yrot, X = X, ss = R2, rotation = A,
                 translation = b, scale = c, xmean = xmean,
                 symmetric = symmetric, svd = sol, call = match.call())
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
  # Check the method
  method = match.arg(method)

  # Orthogonalize U and V
  switch(method,
         "ZCA" = whitening.zca(sigma),
         "PCA" = whitening.pca(sigma),
         "ZCA-cor" = whitening.zca.cor(sigma),
         "PCA-cor" = whitening.pca.cor(sigma),
         "Cholesky" = whitening.chol(sigma))
}


#' @rdname whitening.matrix
#' @keywords internal
whitening.zca = function(sigma) {
  eigS = eigen(sigma, symmetric = TRUE)
  W = eigS$vectors %*% (t(eigS$vectors) / sqrt(eigS$values))
  return(W)
}

#' @rdname whitening.matrix
#' @keywords internal
whitening.zca.cor = function(sigma) {
  v = sqrt(diag(sigma))
  eigR = eigen(stats::cov2cor(sigma), symmetric = TRUE)
  W = eigR$vectors %*% (t(eigR$vectors) / sqrt(eigR$values)) %*% diag(1/v)
  return(W)
}

#' @rdname whitening.matrix
#' @keywords internal
whitening.pca = function(sigma) {
  eigS = eigen(sigma, symmetric = TRUE)
  eigS$vectors = make.pos.diag(eigS$vectors)
  W = t(eigS$vectors) / sqrt(eigS$values)
  return(W)
}

#' @rdname whitening.matrix
#' @keywords internal
whitening.pca.cor = function(sigma) {
  v = sqrt(diag(sigma))
  R = stats::cov2cor(sigma)
  eigR = eigen(R, symmetric = TRUE)
  eigR$vectors = make.pos.diag(eigR$vectors)
  W = (t(eigR$vectors) / sqrt(eigR$values)) %*% diag(1/v)
  return(W)
}

#' @rdname whitening.matrix
#' @keywords internal
whitening.chol = function(sigma) {
  W = solve(t(chol(sigma)))
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
      s = RSpectra::svds(tcrossprod(U, V), ncomp)
      if (ncomp == 1) {
        U = s$u
        V = s$v * s$d
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
        V = t(solve(W, t(V)))

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


#' @title Orthogonalize the matrices U and V with respect to X and Z
#'
#' @description
#' Orthogonalize \code{[A, U]} and \code{V} with respect to \code{X} and \code{Z},
#' respectively, sequentially applying multivariate least squares and residual
#' whitening on U. The result must satisfy the following contraints:
#' \eqn{X^\top A = 0}, \eqn{X^\top U = 0}, \eqn{Z^\top V = 0}, \eqn{U^\top U = 0}.
#'
#' @keywords internal
orthogonalize = function (X, Z, B, A, U, V, method = c("QR", "SVD", "ZCA", "ZCA-cor", "PCA", "PCA-cor", "Cholesky")) {

  # Parameter dimension
  p = ncol(X)
  q = ncol(Z)
  d = ncol(U)

  # Orthogonalize A and U wrt X, and update B
  DU = solve(crossprod(X), crossprod(X, cbind(A, U)))
  A = A - X %*% DU[,  (1:q)]
  U = U - X %*% DU[,q+(1:d)]
  B = B + tcrossprod(cbind(Z, V), DU)
  rm(DU); gc()

  # Orthogonalize V wrt Z, and update A
  DV = solve(crossprod(Z), crossprod(Z, V))
  V = V - Z %*% DV
  A = A + tcrossprod(U, DV)
  rm(DV); gc()

  # Orthogonalize U and V
  QR = orthogonalize.uv(U, V, method=method)
  # qrU = qr(U)
  # U = qr.Q(qrU)
  # V = tcrossprod(V, qr.R(qrU))
  # rm(qrU); gc()

  # Output
  list(B = B, A = A, U = QR$U, V = QR$V)
}

#' @title Normalize the matrices U and V
#'
#' @description
#' Rotate U and V using either QR or SVD decompositions.
#'
#' @details
#' Orthogonalization is implemented using the following methods:
#' \itemize{
#'   \item \code{method = "SVD"}: orthogonal \eqn{U} and scaled orthogonal \eqn{V} based on SVD decomposition;
#'   \item \code{method = "QR"}: orthogonal \eqn{U} and lower triangular \eqn{V} based on QR decomposition;
#'   \item \code{method = "ZCA"}: standardized \eqn{U} and lower triangular \eqn{V} based on ZCA whitening and QR decomposition;
#'   \item \code{method = "ZCA-cor"}: uncorrelated \eqn{U} and lower triangular \eqn{V} based on ZCA whitening and QR decomposition;
#'   \item \code{method = "PCA"}: standardized \eqn{U} and lower triangular \eqn{V} based on PCA whitening and QR decomposition;
#'   \item \code{method = "PCA-cor"}: uncorrelated \eqn{U} and lower triangular \eqn{V} based on PCA whitening and QR decomposition;
#'   \item \code{method = "Cholesky"}: standardized \eqn{U} and lower triangular \eqn{V} based on Cholesky whitening and QR decomposition.
#' }
#'
#' @keywords internal
orthogonalize.uv = function(U, V, method = c("QR", "SVD", "ZCA", "ZCA-cor", "PCA", "PCA-cor", "Cholesky")) {
  # Check the method
  method = match.arg(method)

  # Orthogonalize U and V
  switch(method,
         "SVD" = orthogonalize.svd(U, V),
         "QR" = orthogonalize.qr(U, V),
         "ZCA" = orthogonalize.std(U, V, method),
         "PCA" = orthogonalize.std(U, V, method),
         "ZCA-cor" = orthogonalize.std(U, V, method),
         "PCA-cor" = orthogonalize.std(U, V, method),
         "Cholesky" = orthogonalize.std(U, V, method))
}

#' @rdname orthogonalize.uv
#' @keywords internal
orthogonalize.svd = function(U, V) {
  try({
    # SVD of U*Vt
    ncomp = ncol(U)
    s = RSpectra::svds(tcrossprod(U, V), ncomp)
    if (ncomp == 1) {
      # Orthogonalize U and V
      U = sign(s$v[,1]) * s$u
      V = sign(s$v[,1]) * s$v * s$d
    } else {
      # Orthogonalize U and V
      U = s$u
      V = s$v %*% diag(s$d)

      # Fix the sign of U and V
      D = diag(V)
      V = t(sign(D) * t(V))
      U = t(sign(D) * t(U))
    }
  })
  list(U = U, V = V)
}

#' @rdname orthogonalize.uv
#' @keywords internal
orthogonalize.qr = function(U, V) {
  try({
    if (ncol(U) == 1) {
      # Normalize U
      normU = sqrt(sum(U^2))
      signV = sign(V[,1])
      U = signV * U / normU
      V = signV * V * normU
    } else {
      # Orthogonalize U
      qrU = qr(U)
      U = qr.Q(qrU)
      V = tcrossprod(V, qr.R(qrU))

      # Triangularize V
      qrV = qr(t(V))
      U = U %*% qr.Q(qrV)
      V = t(qr.R(qrV))

      # Fix the sign of U and V
      D = diag(V)
      V = t(sign(D) * t(V))
      U = t(sign(D) * t(U))
    }
  })
  list(U = U, V = V)
}

#' @rdname orthogonalize.uv
#' @keywords internal
orthogonalize.std = function(U, V, method) {
  if (is.vector(U) | (is.matrix(U) & ncol(U) == 1)){
    sdU = stats::sd(c(U))
    U = U / sdU
    V = V * sdU
  }
  if (is.matrix(U) & ncol(U) > 1) {
    try({
      # Compute the covariance of U
      S = stats::cov(U)

      # Make the covariance of U identity
      W = whitening.matrix(S, method=method)
      U = U %*% W
      V = t(solve(W, t(V)))

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
#' oldpar = par(no.readonly = TRUE)
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
#' par(oldpar)
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
  phi = ifelse(!is.null(dispersion), abs(dispersion), 3*stats::rbeta(1, 2, 3))

  # Family name
  fname = family$family
  if (fname == "Gamma") fname = "gamma"
  if (fname == "inverse.gaussian") fname = "invgaussian"
  if (substring(family$family, 1, 17) == "Negative Binomial") fname = "negbinom"

  # Simulate the data using an dispersion exponential family distribution
  Y = matrix(
    switch(fname,
      "gaussian" = stats::rnorm(n * m, mean = mu, sd = sqrt(phi)),
      "binomial" = stats::rbinom(n * m, size = 1, prob = mu),
      "poisson" = stats::rpois(n * m, lambda = mu),
      "gamma" = stats::rgamma(n * m, shape = 1 / phi, scale = phi * mu),
      "invgaussian" = SuppDists::rinvGauss(n * m, nu = mu, lambda = 1 / phi),
      "negbinom" = MASS::rnegbin(n * m, mu = mu, theta = phi)),
    nrow = n, ncol = m)

  # Return the simulated signals
  list(Y = Y, U = U, V = V,
       eta = eta, mu = mu, phi = phi,
       ncomp = ncomp, family = family)
}

