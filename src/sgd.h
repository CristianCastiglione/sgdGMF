// sgd.h
// author: Cristian Castiglione
// creation: 05/10/2023
// last change: 10/10/2023

#ifndef SGD_H
#define SGD_H

#include <RcppArmadillo.h>

#include <string>
#include <memory>   // for dynamic pointers management with std::unique_ptr and std::make_unique 
#include <thread>   // for checking the number of cores with std::thread::hardware_concurrency()
#include <time.h>   // for checking the CPU clocks and the execution time
#include <omp.h>    // for parrallelizing the code via openMP

#include "link.h"
#include "family.h"
#include "deviance.h"
#include "utils.h"
#include "misc.h"
#include "minibatch.h"

// Log-likelihood differentials with respect to the linear predictor 
struct dEta {
    arma::mat deta;
    arma::mat ddeta;
    dEta (const int & n, const int & m)
        : deta(n, m), ddeta(n, m) {}
};

// Log-likelihood differentials with respect to the parameters
struct dPar {
    arma::mat dpar;
    arma::mat ddpar;
    dPar (const int & n, const int & m)
        : dpar(n, m), ddpar(n, m) {}
};

// Log-likelihood differentials with respect to the parameters
struct dCube {
    arma::cube dpar;
    arma::cube ddpar;
    dCube (const int & n, const int & m, const int & k)
        : dpar(n, m, k), ddpar(n, m, k) {}
};

// Blockwise-SGD
class BSGD {
    public:
        int maxiter;
        double eps;
        int nafill;
        double tol;
        int size1;
        int size2;
        double burn;
        double rate0;
        double decay;
        double damping;
        double rate1;
        double rate2;
        bool parallel;
        bool verbose;
        int frequency;
        bool progress;

        // Print the class attributes
        void summary ();

        // Update the learning rate at iteration t
        void update_rate (double & rate, const int & iter);
        
        // Update the log-likelihood differentials wrt eta
        void update_deta (
            dEta & deta, const arma::uvec & idx, const arma::uvec & idy, 
            const arma::mat & Y, const arma::mat & eta, const arma::mat & mu, 
            const std::unique_ptr<Family::Family> & family);

        // Update the deviance differentials wrt the parameters
        void update_dpar (
            dPar & dpar, const dEta & deta, 
            const arma::uvec & idx, const arma::uvec & idy, 
            const arma::mat & u, const arma::mat & v, const arma::vec & penalty, 
            const double & scale, const bool & transp);

        // Update the parameter estimates (chunk-wise)
        void update_par (
            arma::mat & par, const dPar & dpar, const double & rate,
            const arma::uvec & idx, const arma::uvec & idy);

        // Smooth the parameter estimates by averaging over iterations
        void smooth_par (
            arma::mat & u, const arma::mat & ut, const int & iter,
            const arma::uvec & idx, const arma::uvec & idy);

        // Model fitting via SGD (1) - SCATTERED UPDATES
        Rcpp::List fit (
            arma::mat & Y, // to fill NA values we need Y to be non-const
            const arma::mat & X, const arma::mat & B, 
            const arma::mat & A, const arma::mat & Z,
            const arma::mat & U, const arma::mat & V,
            const std::unique_ptr<Family::Family> & family, 
            const int & ncomp, const arma::vec & lambda);

        // Model fitting via SGD (2) - SLICEWISE UPDATES
        Rcpp::List fit2 (
            arma::mat & Y, // to fill NA values we need Y to be non-const
            const arma::mat & X, const arma::mat & B, 
            const arma::mat & A, const arma::mat & Z,
            const arma::mat & U, const arma::mat & V,
            const std::unique_ptr<Family::Family> & family, 
            const int & ncomp, const arma::vec & lambda);

        // Class constructor
        BSGD (
            const int & maxiter, const double & eps,
            const int & nafill, const double & tol, const int & size1, 
            const int & size2, const double & burn, const double & rate0, 
            const double & decay, const double & damping, const double & rate1, 
            const double & rate2, const bool & parallel, const bool & verbose, 
            const int & frequency, const bool & progress
        ) {
            if (maxiter > 0) {this->maxiter = maxiter;} else {this->maxiter = 100;}
            if (eps > 0 && eps < 0.5) {this->eps = eps;} else {this->eps = 1e-08;}
            if (nafill > 0) {this->nafill = nafill;} else {this->nafill = 10;}
            if (tol >= 0) {this->tol = tol;} else {this->tol = 1e-05;}
            if (size1 > 0) {this->size1 = size1;} else {this->size1 = 100;}
            if (size2 > 0) {this->size2 = size2;} else {this->size2 = 100;}
            if (burn > 0 && burn <= 1) {this->burn = burn;} else {this->burn = 0.5;}
            if (rate0 > 0) {this->rate0 = rate0;} else {this->rate0 = 0.01;}
            if (decay > 0) {this->decay = decay;} else {this->decay = 1.0;}
            if (damping >= 0) {this->damping = damping;} else {this->damping = 1e-03;}
            if (rate1 > 0 && rate1 < 1) {this->rate1 = rate1;} else {this->rate1 = 0.05;}
            if (rate2 > 0 && rate2 < 1) {this->rate2 = rate2;} else {this->rate2 = 0.01;}
            if (frequency > 0) {this->frequency = frequency;} else {this->frequency = 10;}
            this->parallel = parallel;
            this->verbose = verbose;
            this->progress = progress;
        }
};

// Coordinate-SGD
class CSGD {
    public:
        int maxiter;
        double eps;
        int nafill;
        double tol;
        int size1;
        int size2;
        double burn;
        double rate0;
        double decay;
        double damping;
        double rate1;
        double rate2;
        bool parallel;
        bool verbose;
        int frequency;
        bool progress;

        // Print the class attributes
        void summary ();

        // Update the learning rate at iteration t
        void update_rate (double & rate, const int & iter);
        
        // Update the log-likelihood differentials wrt eta
        void update_deta (
            dEta & deta, const arma::uvec & idx, 
            const arma::mat & Y, const arma::mat & eta, const arma::mat & mu, 
            const std::unique_ptr<Family::Family> & family, const bool & transp);

        // Update the deviance differentials wrt the parameters
        void update_dpar (
            dPar & dpar, const dEta & deta, const arma::uvec & idx,
            const arma::mat & u, const arma::mat & v, const arma::vec & penalty, 
            const double & scale, const bool & transp);

        // Update the parameter estimates (chunk-wise)
        void update_par (
            arma::mat & par, const dPar & dpar, 
            const double & rate, const arma::uvec & idx);

        // Smooth the parameter estimates by averaging over iterations
        void smooth_par (
            arma::mat & u, const arma::mat & ut, 
            const int & iter, const arma::uvec & idx);

        // Model fitting via SGD
        Rcpp::List fit (
            arma::mat & Y, // to fill NA values we need Y to be non-const
            const arma::mat & X, const arma::mat & B, 
            const arma::mat & A, const arma::mat & Z,
            const arma::mat & U, const arma::mat & V,
            const std::unique_ptr<Family::Family> & family, 
            const int & ncomp, const arma::vec & lambda);

        // Class constructor
        CSGD (
            const int & maxiter, const double & eps,
            const int & nafill, const double & tol, const int & size1, 
            const int & size2, const double & burn, const double & rate0, 
            const double & decay, const double & damping, const double & rate1, 
            const double & rate2, const bool & parallel, const bool & verbose, 
            const int & frequency, const bool & progress
        ) {
            if (maxiter > 0) {this->maxiter = maxiter;} else {this->maxiter = 100;}
            if (eps > 0 && eps < 0.5) {this->eps = eps;} else {this->eps = 1e-08;}
            if (nafill > 0) {this->nafill = nafill;} else {this->nafill = 10;}
            if (tol >= 0) {this->tol = tol;} else {this->tol = 1e-05;}
            if (size1 > 0) {this->size1 = size1;} else {this->size1 = 100;}
            if (size2 > 0) {this->size2 = size2;} else {this->size2 = 100;}
            if (burn > 0 && burn <= 1) {this->burn = burn;} else {this->burn = 0.5;}
            if (rate0 > 0) {this->rate0 = rate0;} else {this->rate0 = 0.01;}
            if (decay > 0) {this->decay = decay;} else {this->decay = 1.0;}
            if (damping >= 0) {this->damping = damping;} else {this->damping = 1e-03;}
            if (rate1 > 0 && rate1 < 1) {this->rate1 = rate1;} else {this->rate1 = 0.05;}
            if (rate2 > 0 && rate2 < 1) {this->rate2 = rate2;} else {this->rate2 = 0.01;}
            if (frequency > 0) {this->frequency = frequency;} else {this->frequency = 10;}
            this->parallel = parallel;
            this->verbose = verbose;
            this->progress = progress;
        }
};

// Blockwise-SGD
class MSGD {
    public:
        int maxiter;
        double eps;
        int nafill;
        double tol;
        int size;
        double burn;
        double rate0;
        double decay;
        double damping;
        double rate1;
        double rate2;
        bool parallel;
        bool verbose;
        int frequency;
        bool progress;

        // Print the class attributes
        void summary ();

        // Update the learning rate at iteration t
        void update_rate (double & rate, const int & iter);
        
        // Update the log-likelihood differentials wrt eta
        void update_deta (
            dEta & deta, const arma::uvec & idx,
            const arma::mat & Y, const arma::mat & eta, const arma::mat & mu, 
            const std::unique_ptr<Family::Family> & family);

        // Update the deviance differentials wrt the parameters
        void update_dpar (
            dPar & dpar, const dEta & deta, const arma::uvec & idx, 
            const arma::mat & u, const arma::mat & v, const arma::vec & penalty, 
            const double & scale, const bool & transp);

        // Update the parameter estimates (block-wise)
        void update_par (
            arma::mat & par, const dPar & dpar, 
            const double & rate, const arma::uvec & idy);
        
        // Update the parameter estimates (chunk-wise)
        void update_par (
            arma::mat & par, const dPar & dpar, const double & rate,
            const arma::uvec & idx, const arma::uvec & idy);

        // Smooth the parameter estimates by averaging over iterations
        void smooth_par (
            arma::mat & u, const arma::mat & ut, 
            const int & iter, const arma::uvec & idy);

        // Smooth the parameter estimates by averaging over iterations
        void smooth_par (
            arma::mat & u, const arma::mat & ut, const int & iter,
            const arma::uvec & idx, const arma::uvec & idy);

        // Model fitting via SGD (1) - SLICEWISE UPDATES
        Rcpp::List fit (
            arma::mat & Y, // to fill NA values we need Y to be non-const
            const arma::mat & X, const arma::mat & B, 
            const arma::mat & A, const arma::mat & Z,
            const arma::mat & U, const arma::mat & V,
            const std::unique_ptr<Family::Family> & family, 
            const int & ncomp, const arma::vec & lambda);

        // Model fitting via SGD (2) - PARALLEL UPDATES
        Rcpp::List fit2 (
            arma::mat & Y, // to fill NA values we need Y to be non-const
            const arma::mat & X, const arma::mat & B, 
            const arma::mat & A, const arma::mat & Z,
            const arma::mat & U, const arma::mat & V,
            const std::unique_ptr<Family::Family> & family, 
            const int & ncomp, const arma::vec & lambda);

        // Class constructor
        MSGD (
            const int & maxiter, const double & eps, const int & nafill, 
            const double & tol, const int & size, const double & burn,
            const double & rate0, const double & decay, const double & damping, 
            const double & rate1, const double & rate2, const bool & parallel, 
            const bool & verbose, const int & frequency, const bool & progress
        ) {
            if (maxiter > 0) {this->maxiter = maxiter;} else {this->maxiter = 100;}
            if (eps > 0 && eps < 0.5) {this->eps = eps;} else {this->eps = 1e-08;}
            if (nafill > 0) {this->nafill = nafill;} else {this->nafill = 10;}
            if (tol >= 0) {this->tol = tol;} else {this->tol = 1e-05;}
            if (size > 0) {this->size = size;} else {this->size = 100;}
            if (burn > 0 && burn <= 1) {this->burn = burn;} else {this->burn = 0.5;}
            if (rate0 > 0) {this->rate0 = rate0;} else {this->rate0 = 0.01;}
            if (decay > 0) {this->decay = decay;} else {this->decay = 1.0;}
            if (damping >= 0) {this->damping = damping;} else {this->damping = 1e-03;}
            if (rate1 > 0 && rate1 < 1) {this->rate1 = rate1;} else {this->rate1 = 0.05;}
            if (rate2 > 0 && rate2 < 1) {this->rate2 = rate2;} else {this->rate2 = 0.01;}
            if (frequency > 0) {this->frequency = frequency;} else {this->frequency = 10;}
            this->parallel = parallel;
            this->verbose = verbose;
            this->progress = progress;
        }
};

#endif