#' @keywords internal
"_PACKAGE"

#' @useDynLib sgdGMF, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @import Rcpp
#' @import RcppArmadillo
#' @importFrom stats glm.fit
#' @importFrom svd propack.svd
#' @importFrom RSpectra eigs svds
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %do% %dopar% foreach
NULL
