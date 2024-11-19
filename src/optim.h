// optim.h
// author: Cristian Castiglione
// creation: 05/10/2023
// last change: 19/11/2023

#ifndef OPTIM_H
#define OPTIM_H

#define ARMA_WARN_LEVEL 1

#include <RcppArmadillo.h>

#include <string>
#include <memory>   // for dynamic pointers management with std::unique_ptr and std::make_unique 
#include <time.h>   // for checking the CPU clocks and the execution time
#include <thread>   // for checking the number of cores with std::thread::hardware_concurrency()

#ifdef _OPENMP      // for parrallelizing the code via openMP
    #include <omp.h>
    const bool OMP_CHECK = true;
#else               // optionally define a dummy macro or handle non-parallel versions
    #define omp_get_num_threads() 1
    #define omp_get_thread_num() 0
    #define omp_set_num_threads(x)
    const bool OMP_CHECK = false;
#endif

#include "link.h"
#include "variance.h"
#include "family.h"
#include "deviance.h"
#include "utils.h"
#include "misc.h"
#include "minibatch.h"

using namespace glm;

// Log-likelihood differentials with respect to the linear predictor 
struct dEta {
    arma::mat deta;
    arma::mat ddeta;
    dEta (const int & n, const int & m)
        : deta(n, m), ddeta(n, m) {}
};

// Sufficient statistics for the calculation of the differentials
struct dStat {
    arma::mat eta;
    arma::mat mu;
    arma::mat var;
    arma::mat mueta;
    arma::mat dev;
    dStat (const int & n, const int & m)
        : eta(n, m), mu(n, m), var(n, m), mueta(n, m), dev(n, m) {}
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

// AIRWLS (i.e. row-wise Fisher scoring) optimizer
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
        int nthreads;
        
        // Print the class attributes
        void summary ();

        // Dispersion parameter initialization
        void init_phi (
            double & phi, const int & df, 
            const arma::mat & Y, const arma::mat & weights,
            const arma::mat & mu, const arma::mat & var, 
            const std::unique_ptr<Family> & family);

        // Dispersion parameter update
        void update_phi (
            double & phi, const int & df, 
            const arma::mat & Y, const arma::mat & weights,
            const arma::mat & mu, const arma::mat & var, 
            const std::unique_ptr<Family> & family);

        // Basic weighted least-squares solver for GLM steps
        void glmstep (
            arma::vec & beta, const arma::vec & y, const arma::mat & X,
            const std::unique_ptr<Family> & family, 
            const arma::vec & offset, const arma::vec & weights, 
            const arma::vec & penalty);

        // PIRLS algorithm for GLM fitting
        void glmfit (
            arma::vec & beta, const arma::vec & y, const arma::mat & X,
            const std::unique_ptr<Family> & family, 
            const arma::vec & offset, const arma::vec & weights, 
            const arma::vec & penalty);

        // Sliced updates for penalized V-GLM (sequential implementation)
        void sequential_update (
            arma::mat & beta, const arma::mat & Y, const arma::mat & X,
            const std::unique_ptr<Family> & family, const arma::uvec & idx, 
            const arma::mat & offset, const arma::mat & weights, 
            const arma::vec & penalty, const bool & transp);

        // Sliced updates for penalized V-GLM (parallel implementation)
        void parallel_update (
            arma::mat & beta, const arma::mat & Y, const arma::mat & X,
            const std::unique_ptr<Family> & family, const arma::uvec & idx, 
            const arma::mat & offset, const arma::mat & weights, 
            const arma::vec & penalty, const bool & transp);

        // Sliced updates for penalized V-GLM  
        void update (
            arma::mat & beta, const arma::mat & Y, const arma::mat & X,
            const std::unique_ptr<Family> & family, const arma::uvec & idx, 
            const arma::mat & offset, const arma::mat & weights, 
            const arma::vec & penalty, const bool & transp);

        // Model fitting via AIRWLS
        Rcpp::List fit (
            arma::mat & Y, // to fill NA values we need Y to be non-const
            const arma::mat & X, const arma::mat & B, 
            const arma::mat & A, const arma::mat & Z,
            const arma::mat & U, const arma::mat & V,
            const arma::mat & O, const arma::mat & W,
            const std::unique_ptr<Family> & family, 
            const int & ncomp, const arma::vec & lambda);
        
        // Class constructor
        AIRWLS (
            const int & maxiter, const int & nsteps, 
            const double & stepsize, const double & eps, 
            const int & nafill, const double & tol, 
            const double & damping, const bool & verbose, 
            const int & frequency, const bool & parallel,
            const int & nthreads
        ) {
            if (maxiter > 0) {this->maxiter = maxiter;} else {this->maxiter = 250;}
            if (nsteps > 0) {this->nsteps = nsteps;} else {this->nsteps = 1;}
            if (stepsize > 0) {this->stepsize = stepsize;} else {this->stepsize = 0.1;}
            if (eps >= 0 && eps < 0.5) {this->eps = eps;} else {this->eps = 1e-08;}
            if (nafill > 0) {this->nafill = nafill;} else {this->nafill = 1;}
            if (tol > 0) {this->tol = tol;} else {this->tol = 1e-05;}
            if (damping >= 0) {this->damping = damping;} else {this->damping = 1e-03;}
            if (frequency > 0) {this->frequency = frequency;} else {this->frequency = 25;}
            if (nthreads > 0) {this->nthreads = nthreads;} else {this->nthreads = 1;}
            this->verbose = verbose;
            
            // If OpenMP is not available, switch off parallel by default
            #ifdef _OPENMP
                this->parallel = parallel;
            #else
                this->parallel = false;
                this->nthreads = 1;
            #endif
        }
};

// Quasi-Newton optimizer
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
        bool parallel;
        int nthreads;

        // Print the class attributes
        void summary ();

        // Dispersion parameter initialization
        void init_phi (
            double & phi, const int & df, 
            const arma::mat & Y, const arma::mat & weights,
            const arma::mat & mu, const arma::mat & var, 
            const std::unique_ptr<Family> & family);

        // Dispersion parameter update
        void update_phi (
            double & phi, const int & df, 
            const arma::mat & Y, const arma::mat & weights,
            const arma::mat & mu, const arma::mat & var, 
            const std::unique_ptr<Family> & family);

        void update_dstat (
            dStat & dstat, 
            const arma::mat & Y, const arma::mat & offset,
            const arma::mat & u, const arma::mat & v, 
            const double & lo, const double & up, 
            const std::unique_ptr<Family> & family);
        
        void update_deta (
            dEta & deta, const dStat & dstat, 
            const arma::mat & Y, const arma::mat & weights,
            const std::unique_ptr<Family> & family);

        // Quasi-Newton block update of the parameters (block implementation)
        void blocked_update (
            arma::mat & u, const arma::mat & v, 
            const arma::vec & pen, const arma::uvec & idx,
            const arma::mat & deta, const arma::mat & ddeta);

        // Quasi-Newton block update of the parameters (parallel implementation)
        void parallel_update (
            arma::mat & u, const arma::mat & v, 
            const arma::vec & pen, const arma::uvec & idx,
            const arma::mat & deta, const arma::mat & ddeta);

        // Quasi-Newton block update of the parameters
        void update_par (
            arma::mat & u, const arma::mat & v, 
            const arma::vec & pen, const arma::uvec & idx,
            const arma::mat & deta, const arma::mat & ddeta);

        // Model fitting via quasi-Newton
        Rcpp::List fit (
            arma::mat & Y, // to fill NA values we need Y to be non-const
            const arma::mat & X, const arma::mat & B, 
            const arma::mat & A, const arma::mat & Z,
            const arma::mat & U, const arma::mat & V,
            const arma::mat & O, const arma::mat & W,
            const std::unique_ptr<Family> & family, 
            const int & ncomp, const arma::vec & lambda);
        
        // Class constructor
        Newton (
            const int & maxiter, const double & stepsize, 
            const double & eps, const int & nafill, 
            const double & tol, const double & damping,
            const bool & verbose, const int & frequency,
            const bool & parallel, const int & nthreads
        ) {
            if (maxiter > 0) {this->maxiter = maxiter;} else {this->maxiter = 500;}
            if (stepsize > 0) {this->stepsize = stepsize;} else {this->stepsize = 0.01;}
            if (eps >= 0 && eps < 0.5) {this->eps = eps;} else {this->eps = 1e-08;}
            if (nafill > 0) {this->nafill = nafill;} else {this->nafill = 1;}
            if (tol > 0) {this->tol = tol;} else {this->tol = 1e-05;}
            if (damping >= 0) {this->damping = damping;} else {this->damping = 1e-03;}
            if (frequency > 0) {this->frequency = frequency;} else {this->frequency = 50;}
            if (nthreads > 0) {this->nthreads = nthreads;} else {this->nthreads = 1;}
            this->verbose = verbose;

            // If OpenMP is not available, switch off parallel by default
            #ifdef _OPENMP
                this->parallel = parallel;
            #else
                this->parallel = false;
                this->nthreads = 1;
            #endif
        }
};

// Blockwise-SGD optimizer
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
        int nthreads;
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
            const arma::mat & Y, const arma::mat & weights, 
            const arma::mat & eta, const arma::mat & mu, 
            const std::unique_ptr<Family> & family);

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
        
        // Initialize the dispersion parameter estimate
        void init_phi (
            double & phi, const int & df, const arma::mat & Y, 
            const arma::mat & weights, const arma::mat & mu, 
            const std::unique_ptr<Family> & family);

        // Update and smooth the dispersion parameter estimate
        void update_phi (
            double & phi, const double & rate, 
            const int & nm, const int & df, const arma::mat & Y, 
            const arma::mat & weights, const arma::mat & mu, 
            const arma::uvec & idx, const arma::uvec & idy, 
            const std::unique_ptr<Family> & family);

        // Model fitting via SGD (1) - SCATTERED UPDATES
        Rcpp::List fit (
            arma::mat & Y, // to fill NA values we need Y to be non-const
            const arma::mat & X, const arma::mat & B, 
            const arma::mat & A, const arma::mat & Z,
            const arma::mat & U, const arma::mat & V,
            const arma::mat & O, const arma::mat & W,
            const std::unique_ptr<Family> & family, 
            const int & ncomp, const arma::vec & lambda);

        // Model fitting via SGD (2) - SLICEWISE UPDATES
        Rcpp::List fit2 (
            arma::mat & Y, // to fill NA values we need Y to be non-const
            const arma::mat & X, const arma::mat & B, 
            const arma::mat & A, const arma::mat & Z,
            const arma::mat & U, const arma::mat & V,
            const arma::mat & O, const arma::mat & W,
            const std::unique_ptr<Family> & family, 
            const int & ncomp, const arma::vec & lambda);

        // Class constructor
        BSGD (
            const int & maxiter, const double & eps,
            const int & nafill, const double & tol, const int & size1, 
            const int & size2, const double & burn, const double & rate0, 
            const double & decay, const double & damping, const double & rate1, 
            const double & rate2, const bool & parallel, const int & nthreads, 
            const bool & verbose, const int & frequency, const bool & progress
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
            if (nthreads > 0) {this->nthreads = nthreads;} else {this->nthreads = 1;}
            this->parallel = parallel;
            this->verbose = verbose;
            this->progress = progress;
        }
};

// Coordinate-SGD optimizer
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
        int nthreads;
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
            const arma::mat & Y, const arma::mat & weights, 
            const arma::mat & eta, const arma::mat & mu, 
            const std::unique_ptr<Family> & family, const bool & transp);

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
        
        // Initialize the dispersion parameter estimate
        void init_phi (
            double & phi, const int & df, const arma::mat & Y, 
            const arma::mat & weights, const arma::mat & mu, 
            const std::unique_ptr<Family> & family);

        // Update and smooth the dispersion parameter estimate
        void update_phi (
            double & phi, const double & rate, 
            const int & nm, const int & df, const arma::mat & Y, 
            const arma::mat & weights, const arma::mat & mu, 
            const arma::uvec & idx, const arma::uvec & idy, 
            const std::unique_ptr<Family> & family);

        // Model fitting via SGD
        Rcpp::List fit (
            arma::mat & Y, // to fill NA values we need Y to be non-const
            const arma::mat & X, const arma::mat & B, 
            const arma::mat & A, const arma::mat & Z,
            const arma::mat & U, const arma::mat & V,
            const arma::mat & O, const arma::mat & W,
            const std::unique_ptr<Family> & family, 
            const int & ncomp, const arma::vec & lambda);

        // Class constructor
        CSGD (
            const int & maxiter, const double & eps,
            const int & nafill, const double & tol, const int & size1, 
            const int & size2, const double & burn, const double & rate0, 
            const double & decay, const double & damping, const double & rate1, 
            const double & rate2, const bool & parallel, const int & nthreads, 
            const bool & verbose, const int & frequency, const bool & progress
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
            if (nthreads > 0) {this->nthreads = nthreads;} else {this->nthreads = 1;}
            this->parallel = parallel;
            this->verbose = verbose;
            this->progress = progress;
        }
};

#endif