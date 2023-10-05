// newton.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 02/10/2023

#ifndef NEWTON_H
#define NEWTON_H

#include <RcppArmadillo.h>
#include <memory>
#include <string>
#include <time.h>
#include "link.h"
#include "family.h"
#include "deviance.h"
#include "utils.h"
#include "misc.h"


class Newton {
    public:
        int maxiter;
        double stepsize;
        double eps;
        int nafill;
        double tol;
        double damping;
        bool verbose;
        int frequency;

        // Print the class attributs
        void summary ();

        // Quasi-Newton block update of the parameters
        void update (
            arma::mat & u, const arma::mat & v, 
            const arma::vec & pen, const arma::uvec & idx,
            const arma::mat & deta, const arma::mat & ddeta);

        // Model fitting via quasi-Newton
        Rcpp::List fit (
            arma::mat & Y, // to fill NA values we need Y to be non-const
            const arma::mat & X, const arma::mat & B, 
            const arma::mat & A, const arma::mat & Z,
            const arma::mat & U, const arma::mat & V,
            const std::unique_ptr<Family::Family> & family, 
            const int & ncomp, const arma::vec & lambda);
        
        // Class constructor
        Newton (
            const int & maxiter, const double & stepsize, 
            const double & eps, const int & nafill, 
            const double & tol, const double & damping,
            const bool & verbose, const int & frequency
        ) {
            if (maxiter > 0) {this->maxiter = maxiter;} else {this->maxiter = 500;}
            if (stepsize > 0) {this->stepsize = stepsize;} else {this->stepsize = 0.01;}
            if (eps >= 0 && eps < 0.5) {this->eps = eps;} else {this->eps = 1e-08;}
            if (nafill > 0) {this->nafill = nafill;} else {this->nafill = 1;}
            if (tol > 0) {this->tol = tol;} else {this->tol = 1e-05;}
            if (damping >= 0) {this->damping = damping;} else {this->damping = 1e-03;}
            if (frequency > 1) {this->frequency = frequency;} else {this->frequency = 50;}
            this->verbose = verbose;
        }
};

#endif