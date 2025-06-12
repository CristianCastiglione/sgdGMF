
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
#' @return An array containing the simulated data.
#'
#' @export
simulate = function (object, ...) UseMethod("simulate")


#' @title Store data into an object
#'
#' @description
#' Generic function to store data into an object, typically a statistical model
#'
#' @param object an object from which simulate new data
#' @param ... additional arguments passed to or from other methods
#'
#' @return An object of the same class as the input containing new data
#'
#' @export
storedata = function (object, ...) UseMethod("storedata")
