// optim_export.cpp
// author: Cristian Castiglione
// creation: 03/02/2024
// last change: 03/02/2024

#include "optim.h"

using namespace glm;

//' @title Compute one Fisher scoring step for GLMs
//' 
//' @description
//' Internal function to compute one Fisher scoring step for GLMs.
//' It constitutes the building block of the AIRWLS algorithm for the
//' estimation of GMF models.
//' 
//' @param beta current value of the regression coefficients to be updated
//' @param y response vector
//' @param X design matrix
//' @param familyname model family name
//' @param linkname link function name
//' @param offset vector of constants to be added to the linear predictor
//' @param penalty penalty parameter of a ridge-type penalty
//' 
//' @keywords internal
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

//' @title Fisher scoring algorithm for GLMs
//' 
//' @description
//' Internal function implementing the Fisher scoring algorithms for the
//' estimation of GLMs. It is used in the AIRWLS algorithm for the 
//' estimation of GMF models.
//' 
//' @param beta initial value of the regression coefficients to be estimated
//' @param y response vector
//' @param X design matrix
//' @param familyname model family name
//' @param linkname link function name
//' @param offset vector of constants to be added to the linear predictor
//' @param penalty penalty parameter of a ridge-type penalty
//' @param nsteps number of iterations
//' @param stepsize stepsize parameter of the Fisher scoring algorithm
//' @param print if \code{TRUE}, print the algorithm history
//' 
//' @keywords internal
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

    // GLM fit via PIRLS
    // arma::vec coef = beta;
    const int p = X.n_cols;
    arma::vec coef = arma::solve(X.t() * X + 0.1 * arma::eye(p,p), X.t() * family->initialize(y));
    airwls.glmfit(coef, y, X, family, offset, penalty);

    return coef;    
}

//' @title AIRWLS update for GMF models
//' 
//' @description
//' Internal function implementing one step of AIRWLS for the
//' estimation of GMF models. 
//' 
//' @param beta initial value of the regression coefficients to be estimated
//' @param Y response vector
//' @param X design matrix
//' @param familyname model family name
//' @param linkname link function name
//' @param idx index identifying the parameters to be updated in \code{beta}
//' @param offset vector of constants to be added to the linear predictor
//' @param penalty penalty parameter of a ridge-type penalty
//' @param transp if \code{TRUE}, transpose the data
//' @param nsteps number of iterations
//' @param stepsize stepsize parameter of the Fisher scoring algorithm
//' @param print if \code{TRUE}, print the algorithm history
//' @param parallel if \code{TRUE}, run the updates in parallel using \code{openMP}
//' @param nthreads number of threads to be run in parallel (only if \code{parallel=TRUE})
//' 
//' @keywords internal
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

//' @title Fit a GMF model using the AIRWLS algorithm
//'
//' @description Fit a GMF model using the AIRWLS algorithm
//'
//' @param Y matrix of responses (\eqn{n \times m})
//' @param X matrix of row fixed effects (\eqn{n \times p})
//' @param B initial row-effect matrix (\eqn{n \times p})
//' @param A initial column-effect matrix (\eqn{n \times q})
//' @param Z matrix of column fixed effects (\eqn{m \times q})
//' @param U initial factor matrix (\eqn{n \times d})
//' @param V initial loading matrix (\eqn{m \times d})
//' @param familyname a \code{glm} model family name
//' @param linkname a \code{glm} link function name
//' @param ncomp rank of the latent matrix factorization
//' @param lambda penalization parameters
//' @param maxiter maximum number of iterations
//' @param nsteps number of inner Fisher scoring iterations
//' @param stepsize stepsize of the inner Fisher scoring algorithm
//' @param eps shrinkage factor for extreme predictions
//' @param nafill how often the missing values are updated
//' @param tol tolerance threshold for the stopping criterion
//' @param damping diagonal dumping factor for the Hessian matrix
//' @param verbose if \code{TRUE}, print the optimization status
//' @param frequency how often the optimization status is printed
//' @param parallel if \code{TRUE}, allows for parallel computing
//' @param nthreads number of cores to be used in parallel
//'
//' @keywords internal
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

//' @title Fit a GMF model using the diagonal quasi-Newton algorithm
//'
//' @description Fit a GMF model using the diagonal quasi-Newton algorithm
//'
//' @param Y matrix of responses (\eqn{n \times m})
//' @param X matrix of row fixed effects (\eqn{n \times p})
//' @param B initial row-effect matrix (\eqn{n \times p})
//' @param A initial column-effect matrix (\eqn{n \times q})
//' @param Z matrix of column fixed effects (\eqn{m \times q})
//' @param U initial factor matrix (\eqn{n \times d})
//' @param V initial loading matrix (\eqn{m \times d})
//' @param familyname a \code{glm} model family name
//' @param linkname a \code{glm} link function name
//' @param ncomp rank of the latent matrix factorization
//' @param lambda penalization parameters
//' @param maxiter maximum number of iterations
//' @param stepsize stepsize of the quasi-Newton update
//' @param eps shrinkage factor for extreme predictions
//' @param nafill how often the missing values are updated
//' @param tol tolerance threshold for the stopping criterion
//' @param damping diagonal dumping factor for the Hessian matrix
//' @param verbose if \code{TRUE}, print the optimization status
//' @param frequency how often the optimization status is printed
//' @param parallel if \code{TRUE}, allows for parallel computing
//' @param nthreads number of cores to be used in parallel
//'
//' @keywords internal
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

//' @title Fit a GMF model using the adaptive SGD with coordinate-wise minibatch subsampling algorithm
//'
//' @description Fit a GMF model using the adaptive SGD with coordinate-wise minibatch subsampling algorithm
//'
//' @param Y matrix of responses (\eqn{n \times m})
//' @param X matrix of row fixed effects (\eqn{n \times p})
//' @param B initial row-effect matrix (\eqn{n \times p})
//' @param A initial column-effect matrix (\eqn{n \times q})
//' @param Z matrix of column fixed effects (\eqn{m \times q})
//' @param U initial factor matrix (\eqn{n \times d})
//' @param V initial loading matrix (\eqn{m \times d})
//' @param familyname a \code{glm} model family name
//' @param linkname a \code{glm} link function name
//' @param ncomp rank of the latent matrix factorization
//' @param lambda penalization parameters
//' @param maxiter maximum number of iterations
//' @param eps shrinkage factor for extreme predictions
//' @param nafill how often the missing values are updated
//' @param tol tolerance threshold for the stopping criterion
//' @param size1 row-minibatch dimension
//' @param size2 column-minibatch dimension
//' @param burn burn-in period in which the learning late is not decreased
//' @param rate0 initial learning rate
//' @param decay decay rate of the learning rate
//' @param damping diagonal dumping factor for the Hessian matrix
//' @param rate1 decay rate of the first moment estimate of the gradient
//' @param rate2 decay rate of the second moment estimate of the gradient
//' @param parallel if \code{TRUE}, allows for parallel computing
//' @param nthreads number of cores to be used in parallel
//' @param verbose if \code{TRUE}, print the optimization status
//' @param frequency how often the optimization status is printed
//' @param progress if \code{TRUE}, print an progress bar
//' 
//' @keywords internal
// [[Rcpp::export("cpp.fit.coord.sgd")]]
Rcpp::List cpp_fit_coord_sgd (
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

//' @title Fit a GMF model using the adaptive SGD with block-wise minibatch subsampling
//'
//' @description Fit a GMF model using the adaptive SGD with block-wise minibatch subsampling
//'
//' @param Y matrix of responses (\eqn{n \times m})
//' @param X matrix of row fixed effects (\eqn{n \times p})
//' @param B initial row-effect matrix (\eqn{n \times p})
//' @param A initial column-effect matrix (\eqn{n \times q})
//' @param Z matrix of column fixed effects (\eqn{m \times q})
//' @param U initial factor matrix (\eqn{n \times d})
//' @param V initial loading matrix (\eqn{m \times d})
//' @param familyname a \code{glm} model family name
//' @param linkname a \code{glm} link function name
//' @param ncomp rank of the latent matrix factorization
//' @param lambda penalization parameters
//' @param maxiter maximum number of iterations
//' @param eps shrinkage factor for extreme predictions
//' @param nafill how often the missing values are updated
//' @param tol tolerance threshold for the stopping criterion
//' @param size1 row-minibatch dimension
//' @param size2 column-minibatch dimension
//' @param burn burn-in period in which the learning late is not decreased
//' @param rate0 initial learning rate
//' @param decay decay rate of the learning rate
//' @param damping diagonal dumping factor for the Hessian matrix
//' @param rate1 decay rate of the first moment estimate of the gradient
//' @param rate2 decay rate of the second moment estimate of the gradient
//' @param parallel if \code{TRUE}, allows for parallel computing
//' @param nthreads number of cores to be used in parallel
//' @param verbose if \code{TRUE}, print the optimization status
//' @param frequency how often the optimization status is printed
//' @param progress if \code{TRUE}, print an progress bar
//' 
//' @keywords internal
// [[Rcpp::export("cpp.fit.block.sgd")]]
Rcpp::List cpp_fit_block_sgd (
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
