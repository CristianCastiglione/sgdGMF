// test-utils.cpp
// author: Cristian Castiglione
// creation: 29/09/2023
// last change: 29/09/2023

#include "utils.h"

// [[Rcpp::export("cpp.dabsmax")]]
double cpp_dabsmax (const double & u, const double & v) {return utils::absmax(u, v);}

// [[Rcpp::export("cpp.vabsmax")]]
double cpp_vabsmax (const arma::vec & u, const arma::vec & v) {return utils::absmax(u, v);}

// [[Rcpp::export("cpp.trim")]]
arma::vec cpp_trim (const arma::vec & x, double a, double b) {
    arma::vec y = x;
    utils::trim(y, a, b);
    return y;
}

// [[Rcpp::export("cpp.xlogx")]]
arma::vec cpp_xlogx (const arma::vec & x) {return utils::xlogx(x);}

// [[Rcpp::export("cpp.log1pexp")]]
arma::vec cpp_log1pexp (const arma::vec & x) {return utils::log1pexp(x);}

// [[Rcpp::export("cpp.log1mexp")]]
arma::vec cpp_log1mexp (const arma::vec & x) {return utils::log1mexp(x);}

// [[Rcpp::export("cpp.logit")]]
arma::vec cpp_logit (const arma::vec & x) {return utils::logit(x);}

// [[Rcpp::export("cpp.expit")]]
arma::vec cpp_expit (const arma::vec & x) {return utils::expit(x);}

// [[Rcpp::export("cpp.expit2")]]
arma::vec cpp_expit2 (const arma::vec & x) {return utils::expit2(x);}

// [[Rcpp::export("cpp.expitn")]]
arma::vec cpp_expitn (const arma::vec & x, double n = 1) {return utils::expitn(x, n);}

// [[Rcpp::export("cpp.cloglog")]]
arma::vec cpp_cloglog (const arma::vec & x) {return utils::cloglog(x);}

// [[Rcpp::export("cpp.cexpexp")]]
arma::vec cpp_cexpexp (const arma::vec & x) {return utils::cexpexp(x);}

// [[Rcpp::export("cpp.loglog")]]
arma::vec cpp_loglog (const arma::vec & x) {return utils::loglog(x);}

// [[Rcpp::export("cpp.expexp")]]
arma::vec cpp_expexp (const arma::vec & x) {return utils::expexp(x);}

// [[Rcpp::export("cpp.pdfn")]]
arma::vec cpp_pdfn (const arma::vec & x) {return utils::pdfn(x);}

// [[Rcpp::export("cpp.cdfn")]]
arma::vec cpp_cdfn (const arma::vec & x) {return utils::cdfn(x);}

// [[Rcpp::export("cpp.logpdfn")]]
arma::vec cpp_logpdfn (const arma::vec & x) {return utils::logpdfn(x);}

// [[Rcpp::export("cpp.logcdfn")]]
arma::vec cpp_logcdfn (const arma::vec & x) {return utils::logcdfn(x);}

// [[Rcpp::export("cpp.gamma")]]
arma::vec cpp_gamma (const arma::vec & x) {return utils::gamma(x);}

// [[Rcpp::export("cpp.loggamma")]]
arma::vec cpp_loggamma (const arma::vec & x) {return utils::loggamma(x);}

// [[Rcpp::export("cpp.digamma")]]
arma::vec cpp_digamma (const arma::vec & x) {return utils::digamma(x);}

// [[Rcpp::export("cpp.trigamma")]]
arma::vec cpp_trigamma (const arma::vec & x) {return utils::trigamma(x);}

// [[Rcpp::export("cpp.beta")]]
arma::vec cpp_beta (const arma::vec & x, const arma::vec & y) {return utils::beta(x, y);}

// [[Rcpp::export("cpp.logbeta")]]
arma::vec cpp_logbeta (const arma::vec & x, const arma::vec & y) {return utils::logbeta(x, y);}

// [[Rcpp::export("cpp.dibeta")]]
arma::vec cpp_dibeta (const arma::vec & x, const arma::vec & y) {return utils::dibeta(x, y);}

// [[Rcpp::export("cpp.tribeta")]]
arma::vec cpp_tribeta (const arma::vec & x, const arma::vec & y) {return utils::tribeta(x, y);}

// [[Rcpp::export("cpp.hinge")]]
arma::vec cpp_hinge (const arma::vec & x) {return utils::hinge(x);}

// [[Rcpp::export("cpp.dirac")]]
arma::vec cpp_dirac (const arma::vec & x, double a = 0) {return utils::dirac(x, a);}

// [[Rcpp::export("cpp.step")]]
arma::vec cpp_step (const arma::vec & x, double a = 0, bool lower = true) {return utils::step(x, a, lower);}

// [[Rcpp::export("cpp.vech")]]
arma::vec cpp_vech(const arma::mat & A) {return utils::vech(A);}
