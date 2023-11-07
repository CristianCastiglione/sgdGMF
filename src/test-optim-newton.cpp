// test-newton.cpp
// author: Cristian Castiglione
// creation: 30/09/2023
// last change: 10/10/2023

#include "optim.h"

// [[Rcpp::export]]
Rcpp::List c_fit_newton (
    const arma::mat & Y, 
    const arma::mat & X, const arma::mat & B, 
    const arma::mat & A, const arma::mat & Z,
    const arma::mat & U, const arma::mat & V,
    const std::string & familyname,
    const std::string & linkname, 
    const int & ncomp, 
    const arma::vec & lambda,
    const int & maxiter = 500,
    const double & stepsize = 0.1,
    const double & eps = 1e-08,
    const int & nafill = 1,
    const double & tol = 1e-05,
    const double & damping = 1e-03,
    const bool & verbose = true,
    const int & frequency = 10,
    const bool & parallel = false,
    const int & nthreads = 1
) {
    arma::mat y = Y;

    // Instantiate the parametrized family object
    std::unique_ptr<Family::Family> family = make_family(familyname, linkname);
    
    // Instantiate the Newton optimizer
    Newton newton(maxiter, stepsize, eps, nafill, tol, damping, verbose, frequency, parallel, nthreads);

    // Perform the optimization via Newton algorithm
    Rcpp::List output = newton.fit(y, X, B, A, Z, U, V, family, ncomp, lambda);

    // Return the estimated model
    return output;
}