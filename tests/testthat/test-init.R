# file: test-init.R
# author: Cristian Castiglione
# creation: 05/02/2024
# last change: 04/10/2024

testthat::test_that("OLS initialization", {
  n = 100; m = 10; d = 5; f = Gamma(link = "log")
  dat = sim.gmf.data(n = n, m = m, ncomp = d, family = f, dispersion = 0.5)
  init = sgdgmf.init.ols(dat$Y, ncomp = d, family = f)

  # Output class
  testthat::expect_true(is.list(init))
  # Sub-output classes
  testthat::expect_true(is.matrix(init$U) && is.numeric(init$U))
  testthat::expect_true(is.matrix(init$V) && is.numeric(init$V))
  testthat::expect_true(is.matrix(init$A) && is.numeric(init$A))
  testthat::expect_true(is.matrix(init$B) && is.numeric(init$B))
  # Output dimensions
  testthat::expect_equal(dim(init$U), c(n,d))
  testthat::expect_equal(dim(init$V), c(m,d))
  testthat::expect_equal(dim(init$A), c(n,1))
  testthat::expect_equal(dim(init$B), c(m,1))
})

testthat::test_that("GLM initialization", {
  n = 100; m = 10; d = 5; f = Gamma(link = "log")
  dat = sim.gmf.data(n = n, m = m, ncomp = d, family = f, dispersion = 0.5)
  init = sgdgmf.init.glm(dat$Y, ncomp = d, family = f)

  # Output class
  testthat::expect_true(is.list(init))
  # Sub-output classes
  testthat::expect_true(is.matrix(init$U) && is.numeric(init$U))
  testthat::expect_true(is.matrix(init$V) && is.numeric(init$V))
  testthat::expect_true(is.matrix(init$A) && is.numeric(init$A))
  testthat::expect_true(is.matrix(init$B) && is.numeric(init$B))
  # Output dimensions
  testthat::expect_equal(dim(init$U), c(n,d))
  testthat::expect_equal(dim(init$V), c(m,d))
  testthat::expect_equal(dim(init$A), c(n,1))
  testthat::expect_equal(dim(init$B), c(m,1))
})

testthat::test_that("Random initialization", {
  n = 100; m = 10; d = 5; f = Gamma(link = "log")
  dat = sim.gmf.data(n = n, m = m, ncomp = d, family = f, dispersion = 0.5)
  init = sgdgmf.init.random(dat$Y, ncomp = d, family = f)

  # Output class
  testthat::expect_true(is.list(init))
  # Sub-output classes
  testthat::expect_true(is.matrix(init$U) && is.numeric(init$U))
  testthat::expect_true(is.matrix(init$V) && is.numeric(init$V))
  testthat::expect_true(is.matrix(init$A) && is.numeric(init$A))
  testthat::expect_true(is.matrix(init$B) && is.numeric(init$B))
  # Output dimensions
  testthat::expect_equal(dim(init$U), c(n,d))
  testthat::expect_equal(dim(init$V), c(m,d))
  testthat::expect_equal(dim(init$A), c(n,1))
  testthat::expect_equal(dim(init$B), c(m,1))
})


testthat::test_that("Random initialization", {
  n = 100; m = 10; d = 5; f = Gamma(link = "log")
  dat = sim.gmf.data(n = n, m = m, ncomp = d, family = f, dispersion = 0.5)

  init.ols = sgdgmf.init(dat$Y, ncomp = d, family = f, method = "ols")
  init.glm = sgdgmf.init(dat$Y, ncomp = d, family = f, method = "glm")
  init.rnd = sgdgmf.init(dat$Y, ncomp = d, family = f, method = "random")

  # Output class
  testthat::expect_true(is.list(init.ols))
  testthat::expect_true(is.list(init.glm))
  testthat::expect_true(is.list(init.rnd))
})



