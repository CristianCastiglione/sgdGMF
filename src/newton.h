// newton.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 01/10/2023

#ifndef NEWTON_H
#define NEWTON_H

#include <RcppArmadillo.h>
#include <time.h>
#include "link.h"
#include "family.h"
#include "deviance.h"
#include "utils.h"
#include "misc.h"

namespace Newton {

void check (
    const int & maxiter, const int & stepsize, 
    const double & eps, const int & nafill, 
    const double & tol, const double & damping,
    const bool & verbose, const int & frequency);

void update (
    arma::mat & u, const arma::mat & v, 
    const arma::vec & pen, const arma::uvec & idx,
    const arma::mat & deta, const arma::mat & ddeta, 
    const double & stepsize, const double & damping);

template<class F, class L>
Rcpp::List fit (
    arma::mat & Y, // to fill NA values we need Y to be non-const
    const arma::mat & X, const arma::mat & B, 
    const arma::mat & A, const arma::mat & Z,
    const arma::mat & U, const arma::mat & V,
    const F & family, const L & link,
    const int & ncomp, const arma::vec & lambda,
    const int & maxiter = 500, const int & stepsize = 0.01, 
    const double & eps = 1e-02, const int & nafill = 1, 
    const double & tol = 1e-05, const double & damping = 1e-03,
    const bool & verbose = true, const int & frequency = 10);

}

#endif