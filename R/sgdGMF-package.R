#' @keywords internal
"_PACKAGE"

#' @useDynLib sgdGMF, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @import Rcpp
#' @import RcppArmadillo
#' @importFrom stats glm.fit
#' @importFrom stats family
#' @importFrom stats gaussian
#' @importFrom stats binomial
#' @importFrom stats poisson
#' @importFrom stats Gamma
#' @importFrom stats inverse.gaussian
#' @importFrom stats quasi
#' @importFrom stats quasibinomial
#' @importFrom stats quasipoisson
#' @importFrom MASS neg.bin
#' @importFrom MASS negative.binomial
#' @importFrom RSpectra svds
#' @importFrom RSpectra eigs
#' @importFrom RSpectra eigs_sym
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom stats fitted
#' @importFrom stats predict
#' @importFrom stats coef coefficients
#' @importFrom stats resid residuals
#' @importFrom stats deviance
#' @importFrom stats BIC
#' @importFrom stats deviance
#' @importFrom stats biplot
#' @importFrom stats screeplot
#' @importFrom graphics image
#' @importFrom generics refit
NULL
