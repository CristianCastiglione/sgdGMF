// airwls.h
// author: Cristian Castiglione
// creation: 30/09/2023
// last change: 30/09/2023

#ifndef AIRWLS_H
#define AIRWLS_H

#include <RcppArmadillo.h>
#include <time.h>
#include "link.h"
#include "family.h"
#include "deviance.h"
#include "utils.h"
#include "misc.h"

class AIRWLS {
    public:
        unsigned int maxiter = 500;
        double stepsize = 0.1;
        unsigned int steps = 1;
        double eps = 0.01;
        unsigned int nafill = 1;
        double tol = 1e-05;
        double damping = 1e-03;
        bool verbose = true;
        unsigned int frequency = 10;

        void update (
            arma::mat & u, const arma::mat & v, 
            const arma::vec & pen, const arma::uvec & idx,
            const arma::mat & deta, const arma::mat & ddeta, 
            const double & stepsize, const double & damping);

        template<class F>
        Rcpp::List fit (
            arma::mat & Y, // to fill NA values we need Y to be non-const
            const arma::mat & X, const arma::mat & Z,
            const arma::mat & B, const arma::mat & A,
            const arma::mat & U, const arma::mat & V,
            const F & family, const int & ncomp, const arma::vec & pen);
        
        AIRWLS ();
        AIRWLS (
            const int & maxiter, const int & stepsize, const double & eps,
            const int & nafill, const double & tol, const double & damping,
            const bool & verbose, const int & frequency);
};


#endif