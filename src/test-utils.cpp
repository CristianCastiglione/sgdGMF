// test-utils.cpp
// author: Cristian Castiglione
// creation: 29/09/2023
// last change: 29/09/2023

#include "utils.h"

// [[Rcpp::export]]
double c_dabsmax (const double & u, const double & v) {return utils::absmax(u, v);}

// [[Rcpp::export]]
double c_vabsmax (const arma::vec & u, const arma::vec & v) {return utils::absmax(u, v);}

// [[Rcpp::export]]
arma::vec c_trim (const arma::vec & x, double a, double b) {
    arma::vec y = x;
    utils::trim(y, a, b);
    return y;
}

// [[Rcpp::export]]
arma::vec c_xlogx (const arma::vec & x) {return utils::xlogx(x);}

// [[Rcpp::export]]
arma::vec c_log1pexp (const arma::vec & x) {return utils::log1pexp(x);}

// [[Rcpp::export]]
arma::vec c_log1mexp (const arma::vec & x) {return utils::log1mexp(x);}

// [[Rcpp::export]]
arma::vec c_logit (const arma::vec & x) {return utils::logit(x);}

// [[Rcpp::export]]
arma::vec c_expit (const arma::vec & x) {return utils::expit(x);}

// [[Rcpp::export]]
arma::vec c_expit2 (const arma::vec & x) {return utils::expit2(x);}

// [[Rcpp::export]]
arma::vec c_expitn (const arma::vec & x, double n = 1) {return utils::expitn(x, n);}

// [[Rcpp::export]]
arma::vec c_cloglog (const arma::vec & x) {return utils::cloglog(x);}

// [[Rcpp::export]]
arma::vec c_cexpexp (const arma::vec & x) {return utils::cexpexp(x);}

// [[Rcpp::export]]
arma::vec c_loglog (const arma::vec & x) {return utils::loglog(x);}

// [[Rcpp::export]]
arma::vec c_expexp (const arma::vec & x) {return utils::expexp(x);}

// [[Rcpp::export]]
arma::vec c_pdfn (const arma::vec & x) {return utils::pdfn(x);}

// [[Rcpp::export]]
arma::vec c_cdfn (const arma::vec & x) {return utils::cdfn(x);}

// [[Rcpp::export]]
arma::vec c_logpdfn (const arma::vec & x) {return utils::logpdfn(x);}

// [[Rcpp::export]]
arma::vec c_logcdfn (const arma::vec & x) {return utils::logcdfn(x);}

// [[Rcpp::export]]
arma::vec c_gamma (const arma::vec & x) {return utils::gamma(x);}

// [[Rcpp::export]]
arma::vec c_loggamma (const arma::vec & x) {return utils::loggamma(x);}

// [[Rcpp::export]]
arma::vec c_digamma (const arma::vec & x) {return utils::digamma(x);}

// [[Rcpp::export]]
arma::vec c_trigamma (const arma::vec & x) {return utils::trigamma(x);}

// [[Rcpp::export]]
arma::vec c_beta (const arma::vec & x, const arma::vec & y) {return utils::beta(x, y);}

// [[Rcpp::export]]
arma::vec c_logbeta (const arma::vec & x, const arma::vec & y) {return utils::logbeta(x, y);}

// [[Rcpp::export]]
arma::vec c_dibeta (const arma::vec & x, const arma::vec & y) {return utils::dibeta(x, y);}

// [[Rcpp::export]]
arma::vec c_tribeta (const arma::vec & x, const arma::vec & y) {return utils::tribeta(x, y);}

// [[Rcpp::export]]
arma::vec c_hinge (const arma::vec & x) {return utils::hinge(x);}

// [[Rcpp::export]]
arma::vec c_dirac (const arma::vec & x, double a = 0) {return utils::dirac(x, a);}

// [[Rcpp::export]]
arma::vec c_step (const arma::vec & x, double a = 0, bool lower = true) {return utils::step(x, a, lower);}

// [[Rcpp::export]]
arma::vec c_vech(const arma::mat & A) {return utils::vech(A);}
