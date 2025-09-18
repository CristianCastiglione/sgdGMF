# file: test-utils.R
# author: Cristian Castiglione
# creation: 05/02/2024
# last change: 25/02/2024

testthat::test_that("Whitening matrix", {
  n = 100; m = 10; d = 5

  U = matrix(rnorm(n*d), n, d)
  V = matrix(rnorm(d*d), d, d)
  X = U %*% V
  S = cov(X)
  W.zca = whitening.matrix(cov(X), method = "ZCA")
  W.pca = whitening.matrix(cov(X), method = "PCA")
  W.zca.cor = whitening.matrix(cov(X), method = "ZCA-cor")
  W.pca.cor = whitening.matrix(cov(X), method = "ZCA-cor")
  W.chol = whitening.matrix(cov(X), method = "Cholesky")

  # Check if the whitening matrix corresponds to the inverse covariance
  testthat::expect_equal(cov(X %*% W.zca), diag(d))
  testthat::expect_equal(crossprod(W.zca), solve(S))
  testthat::expect_equal(crossprod(W.pca), solve(S))
  testthat::expect_equal(crossprod(W.zca.cor), solve(S))
  testthat::expect_equal(crossprod(W.pca.cor), solve(S))
  testthat::expect_equal(crossprod(W.chol), solve(S))
})


testthat::test_that("QR orthogonalization", {
  n = 100; m = 10; d = 5

  old = list(U = matrix(rnorm(n*d), n, d), V = matrix(rnorm(m*d), m, d))
  new = normalize.uv(old$U, old$V, method = "qr")

  # Check if the matrix reconstruction is the same
  testthat::expect_equal(tcrossprod(new$U, new$V), tcrossprod(old$U, old$V))
  # Check if the orthogonalized U has identity variance
  testthat::expect_equal(var(new$U), diag(d))
  # Check if the orthogonalized V is lower triangular
  testthat::expect_equal(new$V[upper.tri(new$V)], rep(0, floor(d * (d-1) / 2)))
})


testthat::test_that("SVD orthogonalization", {
  n = 100; m = 10; d = 5

  old = list(U = matrix(rnorm(n*d), n, d), V = matrix(rnorm(m*d), m, d))
  new = normalize.uv(old$U, old$V, method = "svd")
  svds = RSpectra::svds(tcrossprod(new$U, new$V), d)$d

  # Check if the matrix reconstruction is the same
  testthat::expect_equal(tcrossprod(new$U, new$V), tcrossprod(old$U, old$V))
  # Check if U is a scaled orthogonal matrix
  testthat::expect_equal(crossprod(new$U), diag(d))
  # Check if V is a scaled orthogonal matrix
  testthat::expect_equal(crossprod(new$V), diag(colSums(new$V^2)))
  # Check if U and V have the column norm
  testthat::expect_equal(colSums(new$V^2), RSpectra::svds(tcrossprod(new$U, new$V), d)$d^2)
})

testthat::test_that("Column-space orthogonalization", {
  n = 100; m = 10; p = 3; q = 2; d = 5

  X = matrix(rnorm(n*p), n, p)
  Z = matrix(rnorm(m*q), m, q)
  A = matrix(rnorm(n*q), n, q)
  B = matrix(rnorm(m*p), m, p)
  U = matrix(rnorm(n*d), n, d)
  V = matrix(rnorm(m*d), m, d)

  old = list(B = B, A = A, U = U, V = V)
  new = orthogonalize(X, Z, B, A, U, V)

  Yold = tcrossprod(cbind(X, old$A, old$U), cbind(old$B, Z, old$V))
  Ynew = tcrossprod(cbind(X, new$A, new$U), cbind(new$B, Z, new$V))

  testthat::expect_lt(mean(abs(crossprod(X, new$A))), 1e-10)
  testthat::expect_lt(mean(abs(crossprod(X, new$U))), 1e-10)
  testthat::expect_lt(mean(abs(crossprod(Z, new$V))), 1e-10)
  testthat::expect_lt(mean(abs(crossprod(new$U, new$U) - diag(d))), 1e-10)
  testthat::expect_lt(mean(abs(Yold - Ynew)), 1e-10)
})

testthat::test_that("GMF data simulation", {
  n = 100; m = 10; d = 5

  pois = sim.gmf.data(n = n, m = m, ncomp = d, family = poisson())
  bin = sim.gmf.data(n = n, m = m, ncomp = d, family = binomial())
  gam = sim.gmf.data(n = n, m = m, ncomp = d, family = Gamma(link = "log"))

  # Check if the matrices have the correct dimensions
  testthat::expect_equal(dim(pois$Y), c(n,m))
  testthat::expect_equal(dim(bin$Y), c(n,m))
  testthat::expect_equal(dim(gam$Y), c(n,m))

  testthat::expect_equal(dim(pois$eta), c(n,m))
  testthat::expect_equal(dim(bin$eta), c(n,m))
  testthat::expect_equal(dim(gam$eta), c(n,m))

  testthat::expect_equal(dim(pois$mu), c(n,m))
  testthat::expect_equal(dim(bin$mu), c(n,m))
  testthat::expect_equal(dim(gam$mu), c(n,m))

  testthat::expect_equal(dim(pois$U), c(n,d))
  testthat::expect_equal(dim(bin$U), c(n,d))
  testthat::expect_equal(dim(gam$U), c(n,d))

  testthat::expect_equal(dim(pois$V), c(m,d))
  testthat::expect_equal(dim(bin$V), c(m,d))
  testthat::expect_equal(dim(gam$V), c(m,d))

  # Check if the generated data respect their natural constraints
  testthat::expect_true(all(pois$Y >= 0 & is.integer(pois$Y)))
  testthat::expect_true(all(bin$Y %in% c(0,1) & is.integer(pois$Y)))
  testthat::expect_true(all(gam$Y > 0))

  testthat::expect_true(all(pois$mu >= 0))
  testthat::expect_true(all(bin$mu >= 0 & bin$mu <= 1))
  testthat::expect_true(all(gam$mu > 0))
})
