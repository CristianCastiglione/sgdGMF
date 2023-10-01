// test-newton.cpp
// author: Cristian Castiglione
// creation: 30/09/2023
// last change: 30/09/2023

#include "newton.h"
#include "rcppglm.h"

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
    Rcpp::List output;

    // Instantiate the Newton object and perform the optimization
    Newton newton(maxiter, stepsize, eps, nafill, tol, damping, verbose, frequency);
    
    /*
    // Get the correct Family::Family object matching familyname and linkname
    if (familyname == "gaussian") {
        Rcpp::XPtr<Family::Gaussian> family = make_gaussian(linkname);
        output = newton.fit(Y, X, B, A, Z, U, V, *family, ncomp, lambda);
    } else if (familyname == "binomial") {
        Rcpp::XPtr<Family::Binomial> family = make_gaussian(linkname);
        output = newton.fit(Y, X, B, A, Z, U, V, *family, ncomp, lambda);
    } else if (familyname == "poisson") {
        Rcpp::XPtr<Family::Poisson> family = make_gaussian(linkname);
        output = newton.fit(Y, X, B, A, Z, U, V, *family, ncomp, lambda);
    } else if (familyname == "gamma") {
        Rcpp::XPtr<Family::Gamma> family = make_gaussian(linkname);
        output = newton.fit(Y, X, B, A, Z, U, V, *family, ncomp, lambda);
    } else {
        Rcpp::stop("Model family not available.");
    }
    */

    // Get the correct Family::Family object matching familyname and linkname
    if (familyname == "gaussian") {
        Family::Gaussian family = make_gaussian(linkname);
        output = newton.fit(Y, X, B, A, Z, U, V, *family, ncomp, lambda);
    } else if (familyname == "binomial") {
        Family::Binomial family = make_gaussian(linkname);
        output = newton.fit(Y, X, B, A, Z, U, V, *family, ncomp, lambda);
    } else if (familyname == "poisson") {
        Family::Poisson family = make_gaussian(linkname);
        output = newton.fit(Y, X, B, A, Z, U, V, *family, ncomp, lambda);
    } else if (familyname == "gamma") {
        Family::Gamma family = make_gaussian(linkname);
        output = newton.fit(Y, X, B, A, Z, U, V, *family, ncomp, lambda);
    } else {
        Rcpp::stop("Model family not available.");
    }

    // F family
    // switch (familyname) {
    //     case "gaussian": family = make_gaussian(linkname); break;
    //     case "binomial": family = make_binomial(linkname); break;
    //     case "poisson": family = make_poisson(linkname); break;
    //     case "gamma": family = make_gamma(linkname); break;
    //     default: break;
    // }

    // Return the estimated model
    return output;
}