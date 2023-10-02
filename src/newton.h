// newton.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 02/10/2023

#ifndef NEWTON_H
#define NEWTON_H

#include <RcppArmadillo.h>
#include <memory>
#include <time.h>
#include "link.h"
#include "family.h"
#include "deviance.h"
#include "utils.h"
#include "misc.h"


class Newton {
    public:
        int maxiter = 500;
        double stepsize = 0.1;
        double eps = 0.01;
        int nafill = 1;
        double tol = 1e-05;
        double damping = 1e-03;
        bool verbose = true;
        int frequency = 10;

        void update (
            arma::mat & u, const arma::mat & v, 
            const arma::vec & pen, const arma::uvec & idx,
            const arma::mat & deta, const arma::mat & ddeta, 
            const double & stepsize, const double & damping);

        Rcpp::List fit (
            arma::mat & Y, // to fill NA values we need Y to be non-const
            const arma::mat & X, const arma::mat & B, 
            const arma::mat & A, const arma::mat & Z,
            const arma::mat & U, const arma::mat & V,
            const std::unique_ptr<Family::Family> & family, 
            const int & ncomp, const arma::vec & lambda);
        
        // Newton () {}
        Newton (
            const int & maxiter = 500, const int & stepsize = 0.01, 
            const double & eps = 1e-02, const int & nafill = 1, 
            const double & tol = 1e-05, const double & damping = 1e-03,
            const bool & verbose = true, const int & frequency = 10
        ) {
            if (maxiter > 0) {this->maxiter = maxiter;} else {this->maxiter = 500;}
            if (stepsize > 0) {this->stepsize = stepsize;} else {this->stepsize = 0.01;}
            if (eps >= 0 && eps < 0.5) {this->eps = eps;} else {this->eps = 1e-02;}
            if (nafill > 0) {this->nafill = nafill;} else {this->nafill = 1;}
            if (tol > 0) {this->tol = tol;} else {this->tol = 1e-05;}
            if (damping >= 0) {this->damping = damping;} else {this->damping = 1e-03;}
            if (frequency > 1) {this->frequency = frequency;} else {this->frequency = 10;}
            this->verbose = verbose;
        }
};

/*
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
*/

#endif