// test-optim-msgd.cpp
// author: Cristian Castiglione
// creation: 08/10/2023
// last change: 10/10/2023

#include "optim.h"

// [[Rcpp::export]]
Rcpp::List c_fit_msgd (
    const arma::mat & Y, 
    const arma::mat & X, const arma::mat & B, 
    const arma::mat & A, const arma::mat & Z,
    const arma::mat & U, const arma::mat & V,
    const std::string & familyname,
    const std::string & linkname, 
    const int & ncomp, 
    const arma::vec & lambda,
    const int & maxiter = 1000,
    const double & eps = 0.01,
    const int & nafill = 10,
    const double & tol = 1e-08,
    const int & size = 100,
    const double & burn = 0.75,
    const double & rate0 = 0.01,
    const double & decay = 0.01,
    const double & damping = 1e-03,
    const double & rate1 = 0.95,
    const double & rate2 = 0.99,
    const bool & parallel = false,
    const int & nthreads = 1,
    const bool & verbose = true,
    const int & frequency = 250,
    const bool & progress = false
) {
    arma::mat y = Y;

    // Instantiate the parametrized family object
    std::unique_ptr<Family::Family> family = make_family(familyname, linkname);
    
    // Instantiate the Newton optimizer
    MSGD sgd(
        maxiter, eps, nafill, tol, size, burn, rate0, decay, 
        damping, rate1, rate2, parallel, nthreads, verbose, frequency, progress);

    // Perform the optimization via Newton algorithm
    Rcpp::List output = sgd.fit(y, X, B, A, Z, U, V, family, ncomp, lambda);

    // Return the estimated model
    return output;
}
