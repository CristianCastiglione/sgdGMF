// optim_export.cpp
// author: Cristian Castiglione
// creation: 03/02/2024
// last change: 03/02/2024

#include "optim.h"

using namespace glm;

// [[Rcpp::export("cpp.airwls.glmstep")]]
arma::vec cpp_airwls_glmstep (
    const arma::vec & beta, const arma::vec & y, const arma::mat & X,
    const std::string & familyname, const std::string & linkname, 
    const arma::vec & offset, const arma::vec & penalty
) {
    // Instantiate the parametrized family object
    std::unique_ptr<Family> family = make_family(familyname, linkname);

    // Instantiate the AIRWLS optimizer
    bool verbose = false, parallel = false;
    int maxiter = 100, nsteps = 10, nafill = 10, frequency = 25, nthreads = 1;
    double stepsize = 0.1, eps = 1e-08, tol = 1e-05, damping = 1e-03;
    AIRWLS airwls(maxiter, nsteps, stepsize, eps, nafill, tol, damping, verbose, frequency, parallel, nthreads);

    // Update the current parameter estimate via WLS
    arma::vec coef = beta;
    airwls.glmstep(coef, y, X, family, offset, penalty);

    return coef;
}

// [[Rcpp::export("cpp.airwls.glmfit")]]
arma::vec cpp_airwls_glmfit (
    const arma::vec & beta, const arma::vec & y, const arma::mat & X,
    const std::string & familyname, const std::string & linkname, 
    const arma::vec & offset, const arma::vec & penalty,
    const int & nsteps = 100, const double & stepsize = 0.1,
    const bool & print = false
) {
    // Instantiate the parametrized family object
    std::unique_ptr<Family> family = make_family(familyname, linkname);

    // Instantiate the AIRWLS optimizer
    bool verbose = false, parallel = false;
    int maxiter = 100, nafill = 10, frequency = 25, nthreads = 1;
    double eps = 1e-08, tol = 1e-05, damping = 1e-04;
    AIRWLS airwls(maxiter, nsteps, stepsize, eps, nafill, tol, damping, verbose, frequency, parallel, nthreads);
    if (print) {airwls.summary();}

    // GLM fit via PERLS
    // arma::vec coef = beta;
    const int p = X.n_cols;
    arma::vec coef = arma::solve(X.t() * X + 0.1 * arma::eye(p,p), X.t() * family->initialize(y));
    airwls.glmfit(coef, y, X, family, offset, penalty);

    return coef;    
}

// [[Rcpp::export("cpp.airwls.update")]]
arma::mat cpp_airwls_update (
    const arma::mat & beta, const arma::mat & Y, const arma::mat & X,
    const std::string & familyname, const std::string & linkname, 
    const arma::uvec & idx, const arma::mat & offset, 
    const arma::vec & penalty, const bool & transp = false, 
    const int & nsteps = 100, const double & stepsize = 0.1, 
    const bool & print = false, const bool & parallel = false,
    const int & nthreads = 1
) {
    // Instantiate the parametrized family object
    std::unique_ptr<Family> family = make_family(familyname, linkname);

    // Instantiate the AIRWLS optimizer
    bool verbose = false;
    int maxiter = 100, nafill = 10, frequency = 25;
    double eps = 1e-08, tol = 1e-05, damping = 1e-04;
    AIRWLS airwls(maxiter, nsteps, stepsize, eps, nafill, tol, damping, verbose, frequency, parallel, nthreads);
    if (print) {airwls.summary();}

    // GLM fit via PERLS
    arma::mat xtx; 
    arma::mat xty;
    arma::mat coef;
    if (transp) {
        xtx = X.t() * X;
        xty = X.t() * family->initialize(Y).t();
        coef = arma::solve(xtx, xty).t();
    } else {
        xtx = X.t() * X;
        xty = X.t() * family->initialize(Y);
        coef = arma::solve(xtx, xty).t();
    }
    // Rcpp::Rcout << "\n" << arma::size(coef);
    airwls.update(coef, Y, X, family, idx, offset, penalty, transp);

    return coef;    
}

// [[Rcpp::export("cpp.fit.airwls")]]
Rcpp::List cpp_fit_airwls (
    const arma::mat & Y, 
    const arma::mat & X, 
    const arma::mat & B, 
    const arma::mat & A, 
    const arma::mat & Z,
    const arma::mat & U, 
    const arma::mat & V,
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
    const arma::mat & X, 
    const arma::mat & B, 
    const arma::mat & A, 
    const arma::mat & Z,
    const arma::mat & U, 
    const arma::mat & V,
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
    const arma::mat & X, 
    const arma::mat & B, 
    const arma::mat & A, 
    const arma::mat & Z,
    const arma::mat & U, 
    const arma::mat & V,
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


// [[Rcpp::export("cpp.fit.bsgd")]]
Rcpp::List cpp_fit_bsgd (
    const arma::mat & Y, 
    const arma::mat & X, 
    const arma::mat & B, 
    const arma::mat & A, 
    const arma::mat & Z,
    const arma::mat & U, 
    const arma::mat & V,
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
