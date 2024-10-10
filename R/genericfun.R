
#' @title Object eigenvalues
#'
#' @description
#' Generic function to compute and return the eigenvalues of a matrix
#'
#' @param object an object to extract the eivenalues from
#' @param ... additional arguments passed to or from other methods
#'
#' @export
eigenval = function (object, ...) UseMethod("eigenval")

## #' @title Simulate new data
## #'
## #' @description
## #' Generic function to simulate new data from a statistical model
## #'
## #' @param object an object from which simulate new data
## #' @param ... additional arguments passed to or from other methods
## #'
## #' @export
## simulate = function (object, ...) UseMethod("simulate")
