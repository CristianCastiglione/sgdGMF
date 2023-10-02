// deviance.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 28/09/2023

#ifndef DEVIANCE_H
#define DEVIANCE_H

#include <RcppArmadillo.h>
#include <memory>
#include "family.h"

// Pointwise deviance
void deviance (
    arma::mat & dev, const arma::mat & y, const arma::mat & mu, 
    const std::unique_ptr<Family::Family> & family);
arma::mat deviance (
    const arma::mat & y, const arma::mat & mu, 
    const std::unique_ptr<Family::Family> & family);

// Penalty function
void penalty (double & pen, const arma::mat & u, const arma::vec & p);
double penalty (const arma::mat & u, const arma::vec & p);

#endif