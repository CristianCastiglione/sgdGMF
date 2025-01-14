
#' @export
generics::refit

#' @title Simulate new data
#'
#' @description
#' Generic function to simulate new data from a statistical model
#'
#' @param object an object from which simulate new data
#' @param ... additional arguments passed to or from other methods
#'
#' @export
simulate = function (object, ...) UseMethod("simulate")
