# file: test-control.R
# author: Cristian Castiglione
# creation: 05/02/2024
# last change: 04/10/2024

testthat::test_that("Set AIRWLS control parameters", {
  # Empty call
  testthat::expect_true(is.list(set.control.airwls()))
  testthat::expect_true(is.list(set.control.airwls(maxiter = 200, stepsize = 0.5)))
  # Wrongly parametrized call I: right parameter, but wrong value
  testthat::expect_warning(set.control.airwls(stepsize = -1))
  testthat::expect_warning(set.control.airwls(stepsize = TRUE))
  # Wrongly parametrixed call II: inexistent parameter
  testthat::expect_error(set.control.airwls(foo = TRUE))
})

testthat::test_that("Set Newton control parameters", {
  # Empty call
  testthat::expect_true(is.list(set.control.newton()))
  testthat::expect_true(is.list(set.control.newton(maxiter = 200, stepsize = 0.5)))
  # Wrongly parametrized call I: right parameter, but wrong value
  testthat::expect_warning(set.control.newton(stepsize = -1))
  testthat::expect_warning(set.control.newton(stepsize = TRUE))
  # Wrongly parametrixed call II: inexistent parameter
  testthat::expect_error(set.control.newton(foo = TRUE))
})

testthat::test_that("Set C-SGD control parameters", {
  # Empty call
  testthat::expect_true(is.list(set.control.coord.sgd()))
  testthat::expect_true(is.list(set.control.coord.sgd(maxiter = 500, rate0 = 0.5)))
  # Wrongly parametrized call I: right parameter, but wrong value
  testthat::expect_warning(set.control.coord.sgd(rate0 = -1))
  testthat::expect_warning(set.control.coord.sgd(rate0 = TRUE))
  # Wrongly parametrixed call II: inexistent parameter
  testthat::expect_error(set.control.coord.sgd(foo = TRUE))
})

testthat::test_that("Set B-SGD control parameters", {
  # Empty call
  testthat::expect_true(is.list(set.control.block.sgd()))
  testthat::expect_true(is.list(set.control.block.sgd(maxiter = 500, rate0 = 0.5)))
  # Wrongly parametrized call I: right parameter, but wrong value
  testthat::expect_warning(set.control.block.sgd(rate0 = -1))
  testthat::expect_warning(set.control.block.sgd(rate0 = TRUE))
  # Wrongly parametrixed call II: inexistent parameter
  testthat::expect_error(set.control.block.sgd(foo = TRUE))
})

testthat::test_that("Set generic control parameters", {
  ctr.airwls = set.control.alg(method = "airwls", control = list())
  ctr.newton = set.control.alg(method = "newton", control = list())
  ctr.csgd = set.control.alg(method = "sgd", sampling = "coord", control = list())
  ctr.bsgd = set.control.alg(method = "sgd", sampling = "block", control = list())

  testthat::expect_true(is.list(ctr.airwls))
  testthat::expect_true(is.list(ctr.newton))
  testthat::expect_true(is.list(ctr.csgd))
  testthat::expect_true(is.list(ctr.bsgd))
})
