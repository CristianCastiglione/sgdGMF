// test-family.cpp
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 01/10/2023

#include "family.h"

// [[Rcpp::export]]
arma::vec c_gaussian_variance (const arma::vec & mu) {
    Family::Gaussian f; return f.variance(mu);
}

// [[Rcpp::export]]
arma::vec c_gaussian_initialize (const arma::vec & y) {
    Family::Gaussian f; return f.initialize(y);
}

// [[Rcpp::export]]
arma::vec c_gaussian_devresid (const arma::vec & y, const arma::vec & mu) {
    Family::Gaussian f; return f.devresid(y, mu);
}

// [[Rcpp::export]]
arma::vec c_binomial_variance (const arma::vec & mu) {
    Family::Binomial f; return f.variance(mu);
}

// [[Rcpp::export]]
arma::vec c_binomial_initialize (const arma::vec & y) {
    Family::Binomial f; return f.initialize(y);
}

// [[Rcpp::export]]
arma::vec c_binomial_devresid (const arma::vec & y, const arma::vec & mu) {
    Family::Binomial f; return f.devresid(y, mu);
}

// [[Rcpp::export]]
arma::vec c_poisson_variance (const arma::vec & mu) {
    Family::Poisson f; return f.variance(mu);
}

// [[Rcpp::export]]
arma::vec c_poisson_initialize (const arma::vec & y) {
    Family::Poisson f; return f.initialize(y);
}

// [[Rcpp::export]]
arma::vec c_poisson_devresid (const arma::vec & y, const arma::vec & mu) {
    Family::Poisson f; return f.devresid(y, mu);
}

// [[Rcpp::export]]
arma::vec c_gamma_variance (const arma::vec & mu) {
    Family::Gamma f; return f.variance(mu);
}

// [[Rcpp::export]]
arma::vec c_gamma_initialize (const arma::vec & y) {
    Family::Gamma f; return f.initialize(y);
}

// [[Rcpp::export]]
arma::vec c_gamma_devresid (const arma::vec & y, const arma::vec & mu) {
    Family::Gamma f; return f.devresid(y, mu);
}