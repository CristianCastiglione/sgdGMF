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
