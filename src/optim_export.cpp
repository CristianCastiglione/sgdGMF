// optim_export.cpp
// author: Cristian Castiglione
// creation: 03/02/2024
// last change: 03/02/2024

#include "optim.h"

using namespace glm;

// [[Rcpp::export("cpp.fit.airwls")]]
Rcpp::List cpp_fit_airwls (
    const arma::mat & Y, 
    const arma::mat & X, const arma::mat & B, 
    const arma::mat & A, const arma::mat & Z,
    const arma::mat & U, const arma::mat & V,
    const std::string & familyname,
    const std::string & linkname, 
    const int & ncomp, 
    const arma::vec & lambda,
    const int & maxiter = 500,
    const int & nsteps = 1,
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
    std::unique_ptr<Family> family = make_family(familyname, linkname);
    
    // Instantiate the AIRWLS optimizer
    AIRWLS airwls(maxiter, nsteps, stepsize, eps, nafill, tol, damping, verbose, frequency, parallel, nthreads);

    // Perform the optimization via Newton algorithm
    Rcpp::List output = airwls.fit(y, X, B, A, Z, U, V, family, ncomp, lambda);

    // Return the estimated model
    return output;
}

// [[Rcpp::export("cpp.fit.newton")]]
Rcpp::List cpp_fit_newton (
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
    std::unique_ptr<Family> family = make_family(familyname, linkname);
    
    // Instantiate the Newton optimizer
    Newton newton(maxiter, stepsize, eps, nafill, tol, damping, verbose, frequency, parallel, nthreads);

    // Perform the optimization via Newton algorithm
    Rcpp::List output = newton.fit(y, X, B, A, Z, U, V, family, ncomp, lambda);

    // Return the estimated model
    return output;
}

// [[Rcpp::export("cpp.fit.csgd")]]
Rcpp::List cpp_fit_csgd (
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
    const int & size1 = 100,
    const int & size2 = 100,
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
    std::unique_ptr<Family> family = make_family(familyname, linkname);
    
    // Instantiate the Newton optimizer
    CSGD sgd(
        maxiter, eps, nafill, tol, size1, size2, burn, rate0, decay, 
        damping, rate1, rate2, parallel, nthreads, verbose, frequency, progress);

    // Perform the optimization via Newton algorithm
    Rcpp::List output = sgd.fit(y, X, B, A, Z, U, V, family, ncomp, lambda);

    // Return the estimated model
    return output;
}

// [[Rcpp::export("cpp.fit.rsgd")]]
Rcpp::List cpp_fit_rsgd (
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
    const int & size1 = 100,
    const int & size2 = 100,
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
    std::unique_ptr<Family> family = make_family(familyname, linkname);
    
    // Instantiate the Newton optimizer
    RSGD sgd(
        maxiter, eps, nafill, tol, size1, size2, burn, rate0, decay, 
        damping, rate1, rate2, parallel, nthreads, verbose, frequency, progress);

    // Perform the optimization via Newton algorithm
    Rcpp::List output = sgd.fit(y, X, B, A, Z, U, V, family, ncomp, lambda);

    // Return the estimated model
    return output;
}


// [[Rcpp::export("cpp.fit.bsgd")]]
Rcpp::List cpp_fit_bsgd (
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
    const int & size1 = 100,
    const int & size2 = 100,
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
    std::unique_ptr<Family> family = make_family(familyname, linkname);
    
    // Instantiate the Newton optimizer
    BSGD sgd(
        maxiter, eps, nafill, tol, size1, size2, burn, rate0, decay, 
        damping, rate1, rate2, parallel, nthreads, verbose, frequency, progress);

    // Perform the optimization via Newton algorithm
    Rcpp::List output = sgd.fit2(y, X, B, A, Z, U, V, family, ncomp, lambda);

    // Return the estimated model
    return output;
}









