// deviance.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 28/09/2023

#include "newton.h"

void Newton::update (
    arma::mat & U, const arma::mat & V, const arma::vec & pen, 
    const arma::mat & ddiff, const arma::mat & ddratio, 
    const double & stepsize, const double & damping) {

    unsigned int n, m;
    arma::mat dU, ddU;
    n = U.n_rows; m = U.n_cols;
    dU = - ddiff * V + U * arma::diagmat(pen);
    ddU = ddratio * (V % V) + arma::ones(n, m) * arma::diagmat(pen) + damping;
    U = U - stepsize * (dU / ddU);
}

arma::mat Newton::update (
    const arma::mat & U, const arma::mat & V, const arma::vec & pen, 
    const arma::mat & ddiff, const arma::mat & ddratio, 
    const double & stepsize, const double & damping) {
    
    arma::mat Ut = U;
    Newton::update(Ut, V, pen, ddiff, ddratio, stepsize, damping);
    return Ut;
}