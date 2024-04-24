#' @keywords internal
"_PACKAGE"

#' @useDynLib sgdGMF, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @import Rcpp
#' @import RcppArmadillo
#' @importFrom stats glm.fit family
#' @importFrom stats gaussian quasi
#' @importFrom stats binomial quasibinomial
#' @importFrom stats poisson quasipoisson
#' @importFrom stats Gamma inverse.gaussian
#' @importFrom stats biplot screeplot heatmap
#' @importFrom MASS neg.bin negative.binomial
#' @importFrom generics refit
#' @importFrom stats heatmap
#' @importFrom RSpectra svds eigs eigs_sym
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %do% %dopar% foreach
NULL
