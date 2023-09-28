// deviance.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 28/09/2023

#ifndef DEVIANCE_H
#define DEVIANCE_H

#include <RcppArmadillo.h>

// Pointwise deviance
template<class F> void deviance (arma::mat & dev, const arma::mat & y, const arma::mat & mu, F family);
template<class F> arma::mat deviance (const arma::mat & y, const arma::mat & mu, F family);

// Penalty function
void penalty (double & pen, const arma::mat & u, const arma::vec & p);
double penalty (const arma::mat & u, const arma::vec & p);

#endif