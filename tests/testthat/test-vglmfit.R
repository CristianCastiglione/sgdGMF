# file: test-vglmfit.R
# author: Cristian Castiglione
# creation: 23/03/2024
# last change: 24/03/2024

testthat::test_that("Multivariate OLS fitting", {
  n = 100; m = 10; p = 5; q = p+1

  O = matrix(rexp(n*m, rate = 2.0), nrow = n, ncol = m)
  X = cbind(1, matrix(rnorm(n*p, mean = 0.0, sd = 1.0),  nrow = n, ncol = p))
  B = matrix(rnorm(m*q, mean = 0.1, sd = 0.25), nrow = m, ncol = q)
  E = matrix(rnorm(n*m, mean = 0.0, sd = 0.1), nrow = n, ncol = m)
  Y = O + tcrossprod(X, B) + E

  B.hat = ols.fit.coef(Y, X, offset = O)
  F.hat = O + tcrossprod(X, B.hat)
  E.hat = Y - F.hat

  # Check the dimension and the basic properties of the estimates
  testthat::expect_equal(c(m, q), dim(B.hat))
  testthat::expect_equal(crossprod(X, F.hat), crossprod(X, Y - O))
  testthat::expect_equal(matrix(0, q, m), crossprod(X, E.hat))
  testthat::expect_equal(0, mean(E.hat))
})


testthat::test_that("Binomial VGLM fitting", {
  n = 100; m = 10; p = 5; q = p+1
  family = binomial(link = "probit")

  O = matrix(rexp(n*m, rate = 2.0), nrow = n, ncol = m)
  X = cbind(1, matrix(rnorm(n*p, mean = 0.0, sd = 1.0),  nrow = n, ncol = p))
  B = matrix(rnorm(m*q, mean = 0.1, sd = 0.25), nrow = m, ncol = q)
  eta = O + tcrossprod(X, B)
  mu = family$linkinv(eta)
  Y = matrix(rbinom(n*m, size = 1, prob = mu), nrow = n, ncol = m)

  B.hat = vglm.fit.coef(Y, X, family, offset = O, parallel = FALSE)
  eta.hat = O + tcrossprod(X, B.hat)
  mu.hat = family$linkinv(eta.hat)
  dmu.hat = family$mu.eta(eta.hat)
  var.hat = family$variance(mu.hat)
  res.hat = (Y - mu.hat) * dmu.hat / var.hat

  # Check the dimension and the basic properties of the estimates
  testthat::expect_equal(c(m, q), dim(B.hat))
  testthat::expect_true(mean(crossprod(X, res.hat)) < 1e-04)
})

testthat::test_that("Poisson VGLM fitting", {
  n = 100; m = 10; p = 5; q = p+1
  family = poisson(link = "log")

  O = matrix(rexp(n*m, rate = 2.0), nrow = n, ncol = m)
  X = cbind(1, matrix(rnorm(n*p, mean = 0.0, sd = 1.0),  nrow = n, ncol = p))
  B = matrix(rnorm(m*q, mean = 0.1, sd = 0.25), nrow = m, ncol = q)
  eta = O + tcrossprod(X, B)
  mu = family$linkinv(eta)
  Y = matrix(rpois(n*m, lambda = mu), nrow = n, ncol = m)

  B.hat = vglm.fit.coef(Y, X, family, offset = O, parallel = FALSE)
  eta.hat = O + tcrossprod(X, B.hat)
  mu.hat = family$linkinv(eta.hat)
  dmu.hat = family$mu.eta(eta.hat)
  var.hat = family$variance(mu.hat)
  res.hat = (Y - mu.hat) * dmu.hat / var.hat

  # Check the dimension and the basic properties of the estimates
  testthat::expect_equal(c(m, q), dim(B.hat))
  testthat::expect_true(mean(crossprod(X, res.hat)) < 1e-04)
})

testthat::test_that("Gamma VGLM fitting", {
  n = 100; m = 10; p = 5; q = p+1
  family = Gamma(link = "log")

  O = matrix(rexp(n*m, rate = 2.0), nrow = n, ncol = m)
  X = cbind(1, matrix(rnorm(n*p, mean = 0.0, sd = 1.0),  nrow = n, ncol = p))
  B = matrix(rnorm(m*q, mean = 0.1, sd = 0.25), nrow = m, ncol = q)
  eta = O + tcrossprod(X, B)
  mu = family$linkinv(eta)
  Y = matrix(rgamma(n*m, shape = 1, rate = mu), nrow = n, ncol = m)

  B.hat = vglm.fit.coef(Y, X, family, offset = O, parallel = FALSE)
  eta.hat = O + tcrossprod(X, B.hat)
  mu.hat = family$linkinv(eta.hat)
  dmu.hat = family$mu.eta(eta.hat)
  var.hat = family$variance(mu.hat)
  res.hat = (Y - mu.hat) * dmu.hat / var.hat

  # Check the dimension and the basic properties of the estimates
  testthat::expect_equal(c(m, q), dim(B.hat))
  testthat::expect_true(mean(crossprod(X, res.hat)) < 1e-04)
})
