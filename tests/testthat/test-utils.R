# file: test-utils.R
# author: Cristian Castiglione
# creation: 05/02/2024
# last change: 25/02/2024

testthat::test_that("UV orthogonalization", {
  n = 100; m = 10; d = 3

  old = list(
    U = matrix(rnorm(n*d), nrow = n, ncol = d),
    V = matrix(rnorm(m*d), nrow = m, ncol = d))
  new = correct.uv(old$U, old$V)

  # Check if the matrix reconstraction is the same
  testthat::expect_equal(tcrossprod(new$U, new$V), tcrossprod(old$U, old$V))
  # Check if the orthogonalized U has identity variance
  testthat::expect_equal(var(new$U), diag(d))
  # Check if the orthogonalized V is lower triangular
  testthat::expect_equal(new$V[upper.tri(new$V)], rep(0, floor(d * (d-1) / 2)))
})

