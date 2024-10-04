# file: test-eigengap.R
# author: Cristian Castiglione
# creation: 23/03/2024
# last change: 04/10/2024

testthat::test_that("Rank selecion", {
  n = 100; m = 20; d = 5

  # Generate data using Poisson, Binomial and Gamma models
  data_pois = sim.gmf.data(n = n, m = m, ncomp = d, family = poisson())
  data_bin = sim.gmf.data(n = n, m = m, ncomp = d, family = binomial())
  data_gam = sim.gmf.data(n = n, m = m, ncomp = d, family = Gamma(link = "log"), dispersion = 0.25)

  # Initialize the GMF parameters assuming 3 latent factors
  ncomp_pois = sgdgmf.rank(data_pois$Y, family = poisson(), normalize = TRUE)
  ncomp_bin = sgdgmf.rank(data_bin$Y, family = binomial(), normalize = TRUE)
  ncomp_gam = sgdgmf.rank(data_gam$Y, family = Gamma(link = "log"), normalize = TRUE)

  # Output class
  testthat::expect_true(is.numeric(ncomp_pois$ncomp))
  testthat::expect_true(is.numeric(ncomp_bin$ncomp))
  testthat::expect_true(is.numeric(ncomp_gam$ncomp))

  # Output bounds
  testthat::expect_true(ncomp_pois$ncomp > 0 & ncomp_pois$ncomp <= m)
  testthat::expect_true(ncomp_bin$ncomp > 0 & ncomp_bin$ncomp <= m)
  testthat::expect_true(ncomp_gam$ncomp > 0 & ncomp_gam$ncomp <= m)
})
