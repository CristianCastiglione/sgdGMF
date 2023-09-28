// deviance.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 28/09/2023

#ifndef NEWTON_H
#define NEWTON_H

#include <RcppArmadillo.h>

namespace Newton {

void update (
    arma::mat & U, const arma::mat & V, const arma::vec & pen, 
    const arma::mat & ddiff, const arma::mat & ddratio, 
    const double & stepsize, const double & damping);

arma::mat update (
    const arma::mat & U, const arma::mat & V, const arma::vec & pen, 
    const arma::mat & ddiff, const arma::mat & ddratio, 
    const double & stepsize, const double & damping);

}

#endif