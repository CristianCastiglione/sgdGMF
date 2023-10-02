// newton.cpp
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 30/09/2023

#include "newton.h"


void Newton::update (
    arma::mat & u, const arma::mat & v, 
    const arma::vec & pen, const arma::uvec & idx,
    const arma::mat & deta, const arma::mat & ddeta
) {

    const int n = u.n_rows;
    const int m = idx.n_rows;
    arma::mat du(n,m), ddu(n,m);
    du = - deta * v.cols(idx) + u.cols(idx) * arma::diagmat(pen(idx));
    ddu = ddeta * arma::square(v.cols(idx)) + arma::ones(n, m) * arma::diagmat(pen(idx)) + this->damping;
    u.cols(idx) = u.cols(idx) - this->stepsize * (du / ddu);
}

Rcpp::List Newton::fit (
    arma::mat & Y, 
    const arma::mat & X, const arma::mat & B, 
    const arma::mat & A, const arma::mat & Z,
    const arma::mat & U, const arma::mat & V,
    const std::unique_ptr<Family::Family> & family,
    const int & ncomp, const arma::vec & lambda
) {

    // Get the initial CPU time
    clock_t start = clock();

    // Get the data dimensions
    const int n = Y.n_rows; 
    const int m = Y.n_cols; 
    const int d = ncomp; 
    const int p = X.n_cols; 
    const int q = Z.n_cols;
    const double nm = std::sqrt(n * m);

    // Get the range of the data, and the lower and upper bounds
    double mulo, muup, etalo, etaup;
    set_data_bounds(mulo, muup, etalo, etaup, this->eps, Y.min(), Y.max(), family);

    // Build the left and right decomposition matrices
    arma::mat u, v, ut, vt;
    u = arma::join_rows(X, A, U);
    v = arma::join_rows(B, Z, V);
    ut = u;
    vt = v;

    // Build the penalization vectors for the U and V columns
    arma::vec penu, penv;
    set_uv_penalty(penu, penv, lambda, p, q, d);

    // Get the column indices of (A,U) and (B,V)
    arma::uvec idu, idv;
    set_uv_indices(idu, idv, p, q, d);

    // Save the optimization history
    arma::vec state(6);
    arma::mat trace(maxiter, 6);

    // Check if and where there are some NA values
    bool anyna = !Y.is_finite();
    arma::uvec isna = arma::find_nonfinite(Y);

    // Get the linear predictor, the mean and the variance matrices
    arma::mat eta(n,m), mu(n,m);
    eta = u * v.t();
    mu = family->linkinv(eta);
    
    // Truncate all the extreme values
    utils::trim(mu, mulo, muup);
    utils::trim(eta, etalo, etaup);

    // Get the variances and the mu-differentials
    arma::mat var(n,m), mueta(n,m);
    var = family->variance(mu);
    mueta = family->mueta(eta);

    // Fill the missing values with the initial predictions
    Y.elem(isna) = mu.elem(isna);

    // Get the initial deviance, penalty and objective function
    double dev, pen, obj, objt, change;
    dev = arma::accu(deviance(Y, mu, family));
    pen = penalty(u, penu) + penalty(v, penv);
    obj = dev + 0.5 * pen; objt = obj;
    change = INFINITY;

    // Initialize the differentials
    arma::mat deta(n,m), ddeta(n,m);

    // Get the current execution time
    clock_t end = clock();
    double time = exetime(start, end);

    // Save the objective status before starting the optimization loop
    state = arma::vec{.0, dev, pen, obj, change, time};
    trace.row(0) = state.t();

    // Print the optimization state
    if (verbose) {
        std::printf("-------------------------------------------\n");
        std::printf(" Iteration    Deviance   Change   Exe-Time \n");
        print_state(0, dev / nm, 1., time);
    }

    // Optimization loop
    int iter; double diter;
    for (iter = 1; iter < this->maxiter; iter++) {

        // Fill the missing values with the current predictions
        if (iter > 0 && anyna && iter % this->nafill == 0) {
            Y.elem(isna) = mu.elem(isna);
        }

        // Set up helper matrices for computing differentials
        deta = (Y - mu) % mueta / var;
        ddeta = (mueta % mueta) / var;

        // Update U and V elementwise via quasi-Newton
        this->update(ut, v, penu, idu, deta, ddeta);
        this->update(vt, u, penv, idv, deta.t(), ddeta.t());
        u = ut;
        v = vt;

        // Update the linear predictor and the mean matrix
        eta = u * v.t();
        mu = family->linkinv(eta);

        // Truncate all the extreme values
        utils::trim(mu, mulo, muup);
        utils::trim(eta, etalo, etaup);
        
        // Update the variances and the mu-derivatives
        var = family->variance(mu);
        mueta = family->mueta(eta);

        // Update the initial deviance, penalty and objective function
        dev = arma::accu(deviance(Y, mu, family));
        pen = penalty(u, penu) + penalty(v, penv);
        objt = obj; obj = dev + 0.5 * pen;
        change = std::abs(obj - objt) / (std::abs(objt) + 1e-04);

        // Get the current execution time
        end = clock();
        time = exetime(start, end);
        
        // Store the optimization state at the current iteration
        diter = iter;
        state = arma::vec{diter, dev, pen, obj, change, time};
        trace.row(iter) = state.t();

        if (this->verbose && iter % frequency == 0) {
            print_state(iter, dev / nm, change, time);
        }
        
        // Check for convergence
        if (change < this->tol) {break;}
    }

    // Get the final execution time
    end = clock();
    time = exetime(start, end);

    if (this->verbose) {
        print_state(iter, dev / nm, change, time);
        std::printf("-------------------------------------------\n");
    }
    
    // Get the final output
    Rcpp::List output;
    output["method"] = std::string("Newton");
    output["family"] = family->family;
    output["link"] = family->link;
    output["idu"] = idu;
    output["idv"] = idv;
    output["U"] = u;
    output["V"] = v;
    output["eta"] = eta;
    output["mu"] = mu;
    output["var"] = var;
    output["penalty"] = pen;
    output["deviance"] = dev;
    output["objective"] = obj;
    output["exe.time"] = time;
    output["trace"] = trace.rows(0, iter-1);

    // Return the estimated model
    return output;
}
