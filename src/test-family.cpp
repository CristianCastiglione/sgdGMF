// test-family.cpp
// author: Cristian Castiglione
// creation: 29/09/2023
// last change: 29/09/2023

#include "family.h"

// [[Rcpp::export]]
arma::vec c_binomial_logit_loglik (const arma::vec & y, const arma::vec & eta) {
    Link::Logit l;
    Family::Binomial f(l);
    arma::vec mu = f.linkinv(eta);
    arma::vec dev = 0.5 * f.devresid(y, mu);
    return dev;
}


