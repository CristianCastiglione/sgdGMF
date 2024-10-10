// test-utils.cpp
// author: Cristian Castiglione
// creation: 29/09/2023
// last change: 29/09/2023

#include "utils.h"

//' @keywords internal
// [[Rcpp::export("cpp.utils.dabsmax")]]
double cpp_dabsmax (const double & u, const double & v) {return utils::absmax(u, v);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.vabsmax")]]
double cpp_vabsmax (const arma::vec & u, const arma::vec & v) {return utils::absmax(u, v);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.trim")]]
arma::vec cpp_trim (const arma::vec & x, double a, double b) {
    arma::vec y = x;
    utils::trim(y, a, b);
    return y;
}

//' @keywords internal
// [[Rcpp::export("cpp.utils.xlogx")]]
arma::vec cpp_xlogx (const arma::vec & x) {return utils::xlogx(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.log1pexp")]]
arma::vec cpp_log1pexp (const arma::vec & x) {return utils::log1pexp(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.log1mexp")]]
arma::vec cpp_log1mexp (const arma::vec & x) {return utils::log1mexp(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.logit")]]
arma::vec cpp_logit (const arma::vec & x) {return utils::logit(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.expit")]]
arma::vec cpp_expit (const arma::vec & x) {return utils::expit(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.expit2")]]
arma::vec cpp_expit2 (const arma::vec & x) {return utils::expit2(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.expitn")]]
arma::vec cpp_expitn (const arma::vec & x, double n = 1) {return utils::expitn(x, n);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.cloglog")]]
arma::vec cpp_cloglog (const arma::vec & x) {return utils::cloglog(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.cexpexp")]]
arma::vec cpp_cexpexp (const arma::vec & x) {return utils::cexpexp(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.loglog")]]
arma::vec cpp_loglog (const arma::vec & x) {return utils::loglog(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.expexp")]]
arma::vec cpp_expexp (const arma::vec & x) {return utils::expexp(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.pdfn")]]
arma::vec cpp_pdfn (const arma::vec & x) {return utils::pdfn(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.cdfn")]]
arma::vec cpp_cdfn (const arma::vec & x) {return utils::cdfn(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.logpdfn")]]
arma::vec cpp_logpdfn (const arma::vec & x) {return utils::logpdfn(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.logcdfn")]]
arma::vec cpp_logcdfn (const arma::vec & x) {return utils::logcdfn(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.gamma")]]
arma::vec cpp_gamma (const arma::vec & x) {return utils::gamma(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.loggamma")]]
arma::vec cpp_loggamma (const arma::vec & x) {return utils::loggamma(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.digamma")]]
arma::vec cpp_digamma (const arma::vec & x) {return utils::digamma(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.trigamma")]]
arma::vec cpp_trigamma (const arma::vec & x) {return utils::trigamma(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.beta")]]
arma::vec cpp_beta (const arma::vec & x, const arma::vec & y) {return utils::beta(x, y);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.logbeta")]]
arma::vec cpp_logbeta (const arma::vec & x, const arma::vec & y) {return utils::logbeta(x, y);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.dibeta")]]
arma::vec cpp_dibeta (const arma::vec & x, const arma::vec & y) {return utils::dibeta(x, y);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.tribeta")]]
arma::vec cpp_tribeta (const arma::vec & x, const arma::vec & y) {return utils::tribeta(x, y);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.hinge")]]
arma::vec cpp_hinge (const arma::vec & x) {return utils::hinge(x);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.dirac")]]
arma::vec cpp_dirac (const arma::vec & x, double a = 0) {return utils::dirac(x, a);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.step")]]
arma::vec cpp_step (const arma::vec & x, double a = 0, bool lower = true) {return utils::step(x, a, lower);}

//' @keywords internal
// [[Rcpp::export("cpp.utils.vech")]]
arma::vec cpp_vech(const arma::mat & A) {return utils::vech(A);}
