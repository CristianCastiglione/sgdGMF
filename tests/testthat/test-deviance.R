# file: test-deviance.R
# author: Cristian Castiglione
# creation: 05/02/2024
# last change: 25/02/2024

testthat::test_that("Elementwise Gaussian deviance", {
  n = 100; m = 10

  mu = matrix(rnorm(n*m), nrow = n, ncol = m)
  y = matrix(rnorm(n*m, mean = mu, sd = .1), nrow = n, ncol = m)
  dev = pointwise.deviance(mu, y, gaussian())

  testthat::expect_equal(dim(dev), c(n, m))
  testthat::expect_true(all(dev >= 0))
  testthat::expect_true(all(is.finite(dev)))
  testthat::expect_false(anyNA(dev))
})

testthat::test_that("Elementwise Poisson deviance", {
  n = 100; m = 10

  mu = matrix(exp(rnorm(n*m)), nrow = n, ncol = m)
  y = matrix(rpois(n*m, lambda = mu), nrow = n, ncol = m)
  dev = pointwise.deviance(mu, y, poisson())

  testthat::expect_equal(dim(dev), c(n, m))
  testthat::expect_true(all(dev >= 0))
  testthat::expect_true(all(is.finite(dev)))
  testthat::expect_false(anyNA(dev))
})

testthat::test_that("Elementwise Binomial deviance", {
  n = 100; m = 10

  mu = matrix(plogis(rnorm(n*m)), nrow = n, ncol = m)
  y = matrix(rbinom(n*m, size = 1, prob = mu), nrow = n, ncol = m)
  dev = pointwise.deviance(mu, y, binomial())

  testthat::expect_equal(dim(dev), c(n, m))
  testthat::expect_true(all(dev >= 0))
  testthat::expect_true(all(is.finite(dev)))
  testthat::expect_false(anyNA(dev))
})

testthat::test_that("Elementwise Gamma deviance", {
  n = 100; m = 10

  mu = matrix(exp(rnorm(n*m)), nrow = n, ncol = m)
  y = matrix(rgamma(n*m, shape = 1, rate = mu), nrow = n, ncol = m)
  dev = pointwise.deviance(mu, y, Gamma())

  testthat::expect_equal(dim(dev), c(n, m))
  testthat::expect_true(all(dev >= 0))
  testthat::expect_true(all(is.finite(dev)))
  testthat::expect_false(anyNA(dev))
})

testthat::test_that("Elementwise deviance with missing", {
  n = 100; m = 10; f = floor(.3 * n * m)

  mask = unique(cbind(
    sample(1:n, size = f, replace = TRUE),
    sample(1:m, size = f, replace = TRUE)))

  mu = matrix(exp(rnorm(n*m)), nrow = n, ncol = m)
  y = matrix(rgamma(n*m, shape = 1, rate = mu), nrow = n, ncol = m)
  y[mask] = NA

  dev = pointwise.deviance(mu, y, Gamma())

  testthat::expect_equal(dim(dev), c(n, m))
  testthat::expect_true(all(dev[-mask[,1],-mask[,2]] >= 0))
  testthat::expect_true(all(is.finite(dev[-mask[,1],-mask[,2]])))
  testthat::expect_equal(sum(is.na(dev)), nrow(mask))
})

testthat::test_that("Matrix Gaussian deviance", {
  n = 100; m = 10

  mu = matrix(rnorm(n*m), nrow = n, ncol = m)
  y = matrix(rnorm(n*m, mean = mu, sd = .1), nrow = n, ncol = m)
  dev = matrix.deviance(mu, y, gaussian())

  testthat::expect_true(is.finite(dev))
  testthat::expect_true(dev >= 0)
})

testthat::test_that("Matrix Poisson deviance", {
  n = 100; m = 10

  mu = matrix(exp(rnorm(n*m)), nrow = n, ncol = m)
  y = matrix(rpois(n*m, lambda = mu), nrow = n, ncol = m)
  dev = matrix.deviance(mu, y, poisson())

  testthat::expect_true(is.finite(dev))
  testthat::expect_true(dev >= 0)
})

testthat::test_that("Matrix Binomial deviance", {
  n = 100; m = 10

  mu = matrix(plogis(rnorm(n*m)), nrow = n, ncol = m)
  y = matrix(rbinom(n*m, size = 1, prob = mu), nrow = n, ncol = m)
  dev = matrix.deviance(mu, y, binomial())

  testthat::expect_true(is.finite(dev))
  testthat::expect_true(dev >= 0)
})

testthat::test_that("Matrix Gamma deviance", {
  n = 100; m = 10

  mu = matrix(exp(rnorm(n*m)), nrow = n, ncol = m)
  y = matrix(rgamma(n*m, shape = 1, rate = mu), nrow = n, ncol = m)
  dev = matrix.deviance(mu, y, Gamma())

  testthat::expect_true(is.finite(dev))
  testthat::expect_true(dev >= 0)
})


testthat::test_that("Matrix deviance with missing", {
  n = 100; m = 10; f = floor(.3 * n * m)

  mask = unique(cbind(
    sample(1:n, size = f, replace = TRUE),
    sample(1:m, size = f, replace = TRUE)))

  mu = matrix(exp(rnorm(n*m)), nrow = n, ncol = m)
  y = matrix(rgamma(n*m, shape = 1, rate = mu), nrow = n, ncol = m)
  y[mask] = NA

  dev = matrix.deviance(mu, y, Gamma())

  testthat::expect_true(is.finite(dev))
  testthat::expect_false(is.na(dev))
  testthat::expect_true(dev >= 0)
})

testthat::test_that("Frobenious matrix penalty", {
  n = 100; m = 3

  U = matrix(rnorm(n*m), nrow = n, ncol = m)
  lambda = rexp(m)
  pen = matrix.penalty(U, lambda)

  testthat::expect_equal(pen, sum((U * U) %*% diag(lambda)))
})

