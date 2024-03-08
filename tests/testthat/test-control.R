# file: test-control.R
# author: Cristian Castiglione
# creation: 05/02/2024
# last change: 05/02/2024

testthat::test_that("Set AIRWLS control parameters", {
  control = set.control("airwls", list())
})

testthat::test_that("Set Newton control parameters", {
  control = set.control("newton", list())
})

testthat::test_that("Set M-SGD control parameters", {
  control = set.control("msgd", list())
})

testthat::test_that("Set C-SGD control parameters", {
  control = set.control("csgd", list())
})

testthat::test_that("Set R-SGD control parameters", {
  control = set.control("rsgd", list())
})

testthat::test_that("Set B-SGD control parameters", {
  control = set.control("bsgd", list())
})
