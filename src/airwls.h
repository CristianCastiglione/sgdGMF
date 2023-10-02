// airwls.h
// author: Cristian Castiglione
// creation: 02/10/2023
// last change: 02/10/2023

#ifndef AIRWLS_H
#define AIRWLS_H

#include <RcppArmadillo.h>
#include <memory>
#include <time.h>
#include "link.h"
#include "family.h"
#include "deviance.h"
#include "utils.h"
#include "misc.h"


class AIRWLS {
    public:
        int maxiter = 250;
        int nsteps = 1;
        double stepsize = 0.1;
        double eps = 1e-08;
        int nafill = 1;
        double tol = 1e-05;
        double damping = 1e-03;
        bool verbose = true;
        int frequency = 25;

        // Basic weighted least-squares solver for GLM steps
        void wlsfit (
            arma::vec & beta, const arma::vec & y, const arma::mat & X,
            const std::unique_ptr<Family::Family> & family, 
            const arma::vec & offset, const double & penalty);

        // PIRLS algorithm for GLM fitting
        void glmfit (
            arma::vec & beta, const arma::vec & y, const arma::mat & X,
            const std::unique_ptr<Family::Family> & family, 
            const arma::vec & offset, const double & penalty);

        // Sliced updates for penalized V-GLM  
        void update (
            arma::mat & beta, const arma::mat & Y, const arma::mat & X,
            const std::unique_ptr<Family::Family> & family,
            const int & nslices, const arma::mat & offset, 
            const double & penalty);

        // Model fitting via AIRWLS
        Rcpp::List fit (
            arma::mat & Y, // to fill NA values we need Y to be non-const
            const arma::mat & X, const arma::mat & B, 
            const arma::mat & A, const arma::mat & Z,
            const arma::mat & U, const arma::mat & V,
            const std::unique_ptr<Family::Family> & family, 
            const int & ncomp, const arma::vec & lambda);
        
        // Class constructor
        AIRWLS (
            const int & maxiter = 250, const int & nsteps = 1, 
            const int & stepsize = 0.01, const double & eps = 1e-08, 
            const int & nafill = 1, const double & tol = 1e-05, 
            const double & damping = 1e-03, const bool & verbose = true, 
            const int & frequency = 25
        ) {
            if (maxiter > 0) {this->maxiter = maxiter;} else {this->maxiter = 250;}
            if (nsteps > 0) {this->nsteps = nsteps;} else {this->nsteps = 1;}
            if (stepsize > 0) {this->stepsize = stepsize;} else {this->stepsize = 0.01;}
            if (eps >= 0 && eps < 0.5) {this->eps = eps;} else {this->eps = 1e-08;}
            if (nafill > 0) {this->nafill = nafill;} else {this->nafill = 1;}
            if (tol > 0) {this->tol = tol;} else {this->tol = 1e-05;}
            if (damping >= 0) {this->damping = damping;} else {this->damping = 1e-03;}
            if (frequency > 1) {this->frequency = frequency;} else {this->frequency = 25;}
            this->verbose = verbose;
        }
};


#endif