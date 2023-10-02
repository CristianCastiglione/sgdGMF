// test-newton.cpp
// author: Cristian Castiglione
// creation: 30/09/2023
// last change: 01/10/2023

#include "newton.h"

// [[Rcpp::export]]
Rcpp::List c_fit_newton (
    arma::mat & Y, 
    const arma::mat & X, const arma::mat & B, 
    const arma::mat & A, const arma::mat & Z,
    const arma::mat & U, const arma::mat & V,
    const std::string & familyname,
    const std::string & linkname, 
    const int & ncomp, 
    const arma::vec & lambda,
    const int & maxiter = 500,
    const double & stepsize = 0.1,
    const double & eps = 0.01,
    const int & nafill = 1,
    const double & tol = 1e-05,
    const double & damping = 1e-03,
    const bool & verbose = true,
    const int & frequency = 10
) {
    // Instantiate the output object
    Rcpp::List output;
    std::unique_ptr<Family::Family> family = make_family(familyname, linkname);
    // std::unique_ptr<Link::Link> link = std::make_unique<Link::Identity>();
    // std::unique_ptr<Family::Family> family = std::make_unique<Family::Gausian>(link); 
    
    // output = Newton::fit(
    //     Y, X, B, A, Z, U, V, family, link, ncomp, 
    //     lambda, maxiter, stepsize, eps, nafill, 
    //     tol, damping, verbose, frequency);
 
    // Return the estimated model
    return output;
}