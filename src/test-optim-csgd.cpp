// test-optim-csgd.cpp
// author: Cristian Castiglione
// creation: 08/10/2023
// last change: 10/10/2023

#include "optim.h"

// [[Rcpp::export]]
Rcpp::List c_fit_csgd (
    const arma::mat & Y, 
    const arma::mat & X, const arma::mat & B, 
    const arma::mat & A, const arma::mat & Z,
    const arma::mat & U, const arma::mat & V,
    const std::string & familyname,
    const std::string & linkname, 
    const int & ncomp, 
    const arma::vec & lambda,
    int maxiter = 1000,
    double eps = 0.01,
    int nafill = 10,
    double tol = 1e-08,
    int size1 = 100,
    int size2 = 100,
    double burn = 0.75,
    double rate0 = 0.01,
    double decay = 0.01,
    double damping = 1e-03,
    double rate1 = 0.95,
    double rate2 = 0.99,
    bool parallel = false,
    bool verbose = true,
    int frequency = 250,
    bool progress = false
) {
    arma::mat y = Y;

    // Instantiate the parametrized family object
    std::unique_ptr<Family::Family> family = make_family(familyname, linkname);
    
    // Instantiate the Newton optimizer
    CSGD sgd(
        maxiter, eps, nafill, tol, size1, size2, burn, rate0, decay, 
        damping, rate1, rate2, parallel, verbose, frequency, progress);

    // Perform the optimization via Newton algorithm
    Rcpp::List output = sgd.fit(y, X, B, A, Z, U, V, family, ncomp, lambda);

    // Return the estimated model
    return output;
}
