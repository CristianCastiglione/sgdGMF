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
#' @importFrom stats Gamma invese.gaussian
#' @importFrom MASS neg.bin negative.binomial
#' @importFrom svd propack.svd
#' @importFrom RSpectra eigs svds
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %do% %dopar% foreach
NULL
