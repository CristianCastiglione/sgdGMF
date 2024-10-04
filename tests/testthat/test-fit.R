# file: test-fit.R
# author: Cristian Castiglione
# creation: 05/02/2024
# last change: 04/10/2024

testthat::test_that("GMF fit", {
  n = 100; m = 20; d = 5

  # Generate data using Poisson, Binomial and Gamma models
  data_pois = sim.gmf.data(n = n, m = m, ncomp = d, family = poisson())
  data_bin = sim.gmf.data(n = n, m = m, ncomp = d, family = binomial())
  data_gam = sim.gmf.data(n = n, m = m, ncomp = d, family = Gamma(link = "log"), dispersion = 0.25)

  # Initialize the GMF parameters assuming 3 latent factors
  gmf_pois = sgdgmf.fit(data_pois$Y, ncomp = 3, family = poisson())
  gmf_bin = sgdgmf.fit(data_bin$Y, ncomp = 3, family = binomial())
  gmf_gam = sgdgmf.fit(data_gam$Y, ncomp = 3, family = Gamma(link = "log"))

  # Output class
  testthat::expect_true(is.list(gmf_pois))
  testthat::expect_true(is.list(gmf_bin))
  testthat::expect_true(is.list(gmf_gam))

  testthat::expect_s3_class(gmf_pois, "sgdgmf")
  testthat::expect_s3_class(gmf_bin, "sgdgmf")
  testthat::expect_s3_class(gmf_gam, "sgdgmf")

  # Sub-output checks
  testthat::expect_true(is.matrix(gmf_pois$U) && is.numeric(gmf_pois$U))
  testthat::expect_true(is.matrix(gmf_pois$V) && is.numeric(gmf_pois$V))
  testthat::expect_true(is.matrix(gmf_pois$A) && is.numeric(gmf_pois$A))
  testthat::expect_true(is.matrix(gmf_pois$B) && is.numeric(gmf_pois$B))
  testthat::expect_true(is.matrix(gmf_pois$eta) && is.numeric(gmf_pois$eta))
  testthat::expect_true(is.matrix(gmf_pois$mu) && is.numeric(gmf_pois$mu))
  testthat::expect_true(all(gmf_pois$mu >= 0))

  testthat::expect_true(is.matrix(gmf_bin$U) && is.numeric(gmf_bin$U))
  testthat::expect_true(is.matrix(gmf_bin$V) && is.numeric(gmf_bin$V))
  testthat::expect_true(is.matrix(gmf_bin$A) && is.numeric(gmf_bin$A))
  testthat::expect_true(is.matrix(gmf_bin$B) && is.numeric(gmf_bin$B))
  testthat::expect_true(is.matrix(gmf_bin$eta) && is.numeric(gmf_bin$eta))
  testthat::expect_true(is.matrix(gmf_bin$mu) && is.numeric(gmf_bin$mu))
  testthat::expect_true(all(gmf_bin$mu >= 0 & gmf_bin$mu <= 1))

  testthat::expect_true(is.matrix(gmf_gam$U) && is.numeric(gmf_gam$U))
  testthat::expect_true(is.matrix(gmf_gam$V) && is.numeric(gmf_gam$V))
  testthat::expect_true(is.matrix(gmf_gam$A) && is.numeric(gmf_gam$A))
  testthat::expect_true(is.matrix(gmf_gam$B) && is.numeric(gmf_gam$B))
  testthat::expect_true(is.matrix(gmf_gam$eta) && is.numeric(gmf_gam$eta))
  testthat::expect_true(is.matrix(gmf_gam$mu) && is.numeric(gmf_gam$mu))
  testthat::expect_true(all(gmf_gam$mu >= 0))
})
