// airwls.h
// author: Cristian Castiglione
// creation: 02/10/2023
// last change: 05/10/2023

#ifndef AIRWLS_H
#define AIRWLS_H

#include <RcppArmadillo.h>

#include <memory>   // for dynamic pointers management with std::unique_ptr and std::make_unique 
#include <thread>   // for checking the number of cores with std::thread::hardware_concurrency()
#include <time.h>   // for checking the CPU clocks and the execution time
#include <omp.h>    // for parrallelizing the code via openMP

#include "link.h"
#include "family.h"
#include "deviance.h"
#include "utils.h"
#include "misc.h"


class AIRWLS {
    public:
        int maxiter;
        int nsteps;
        double stepsize;
        double eps;
        int nafill;
        double tol;
        double damping;
        bool verbose;
        int frequency;
        bool parallel;
        
        // Print the class attributes
        void summary ();

        // Basic weighted least-squares solver for GLM steps
        void glmstep (
            arma::vec & beta, const arma::vec & y, const arma::mat & X,
            const std::unique_ptr<Family::Family> & family, 
            const arma::vec & offset, const arma::vec & penalty);

        // PIRLS algorithm for GLM fitting
        void glmfit (
            arma::vec & beta, const arma::vec & y, const arma::mat & X,
            const std::unique_ptr<Family::Family> & family, 
            const arma::vec & offset, const arma::vec & penalty);

        // Sliced updates for penalized V-GLM  
        void update (
            arma::mat & beta, const arma::mat & Y, const arma::mat & X,
            const std::unique_ptr<Family::Family> & family,
            const arma::uvec & idx, const arma::mat & offset, 
            const arma::vec & penalty, const bool & transp);

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
            const int & maxiter, const int & nsteps, 
            const double & stepsize, const double & eps, 
            const int & nafill, const double & tol, 
            const double & damping, const bool & verbose, 
            const int & frequency, const bool & parallel
        ) {
            if (maxiter > 0) {this->maxiter = maxiter;} else {this->maxiter = 250;}
            if (nsteps > 0) {this->nsteps = nsteps;} else {this->nsteps = 1;}
            if (stepsize > 0) {this->stepsize = stepsize;} else {this->stepsize = 0.1;}
            if (eps >= 0 && eps < 0.5) {this->eps = eps;} else {this->eps = 1e-08;}
            if (nafill > 0) {this->nafill = nafill;} else {this->nafill = 1;}
            if (tol > 0) {this->tol = tol;} else {this->tol = 1e-05;}
            if (damping >= 0) {this->damping = damping;} else {this->damping = 1e-03;}
            if (frequency > 1) {this->frequency = frequency;} else {this->frequency = 25;}
            this->verbose = verbose;
            this->parallel = parallel;
        }
};


#endif