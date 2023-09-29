// deviance.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 28/09/2023

#ifndef NEWTON_H
#define NEWTON_H

#include <RcppArmadillo.h>
#include <time.h>
#include "deviance.h"
#include "utils.h"

namespace Newton {

template<class F>
Rcpp::List fit (
    const arma::mat & Y, 
    const arma::mat & X, const arma::mat & Z,
    const arma::mat & B, const arma::mat & A,
    const arma::mat & U, const arma::mat & V,
    const F & family, const int & ncomp, const arma::vec & pen,
    const int & maxiter, const double & stepsize, const double & eps,
    const int & nafill, const double & tol, const double & damping, 
    const bool & verbose, const int & frequency);

void update (
    arma::mat & U, const arma::mat & V, 
    const arma::vec & pen, const arma::uvec & idx,
    const arma::mat & deta, const arma::mat & ddeta, 
    const double & stepsize, const double & damping);


}

#endif