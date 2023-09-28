// deviance.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 28/09/2023

#ifndef DEVIANCE_H
#define DEVIANCE_H

#include <RcppArmadillo.h>

// Pointwise deviance
template<class F> void deviance (const arma::mat & mu, const amra::mat & y, F family);
template<class F> arma::mat deviance (const arma::mat & mu, const amra::mat & y, F family);

// Penalty function
void penalty (const arma::mat & u, const arma::vec & p);
double penalty (const arma::mat & u, const arma::vec & p);

#endif