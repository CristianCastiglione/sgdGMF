// optim_newton.cpp
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 10/10/2023

#include "optim.h"

using namespace glm;

void Newton::summary () {
    std::printf("------------------\n");
    std::printf(" maxiter = %i \n", this->maxiter);
    std::printf(" stepsize = %.4f \n", this->stepsize);
    std::printf(" eps = %.4f \n", this->eps);
    std::printf(" nafill = %i \n", this->nafill);
    std::printf(" tol = %.5f \n", this->tol);
    std::printf(" damping = %.5f \n", this->damping);
    std::printf(" parallel = %s \n", this->parallel ? "true" : "false");
    std::printf(" verbose = %s \n", this->verbose ? "true" : "false");
    std::printf(" frequency = %i \n", this->frequency);
    std::printf("------------------\n");
}

void Newton::init_phi (
    double & phi, const int & df, const arma::mat & Y, 
    const arma::mat & mu, const arma::mat & var, 
    const std::unique_ptr<Family> & family
) {
    double ssq;
    if (family->estdisp()) {
        if (family->getfamily() == "NegativeBinomial") {
            ssq = arma::accu(arma::square(Y - mu) - arma::accu(mu));
            phi = std::max(1e-08, ssq / arma::accu(mu % mu));
            family->setdisp(1 / phi);
        } else {
            ssq = arma::accu(arma::square(Y - mu) / var);
            phi = std::max(1e-08, ssq / df);
            family->setdisp(phi);
        }
    }
}

void Newton::update_phi (
    double & phi, const int & df, const arma::mat & Y, 
    const arma::mat & mu, const arma::mat & var, 
    const std::unique_ptr<Family> & family
) {
    double ssq, phit, lambda = 0.5;
    if (family->estdisp()) {
        if (family->getfamily() == "NegativeBinomial") {
            ssq = arma::accu(arma::square(Y - mu));
            phit = (ssq - arma::accu(mu)) / arma::accu(mu % mu);
            phit = std::max(1e-08, phit);
            phi = (1 - lambda) * phi + lambda * phit;
            family->setdisp(1 / phi);
        } else {
            ssq = arma::accu(arma::square(Y - mu) / var);
            phit = std::max(1e-08, ssq / df);
            phi = (1 - lambda) * phi + lambda * phit;
            family->setdisp(phi);
        }
    }
}

void Newton::blocked_update (
    arma::mat & u, const arma::mat & v, 
    const arma::vec & pen, const arma::uvec & idx,
    const arma::mat & deta, const arma::mat & ddeta
) {
    // Block update via elementwise matrix operations
    const unsigned int n = u.n_rows;
    const unsigned int m = idx.n_rows;
    arma::mat du(n,m), ddu(n,m);
    du = - deta * v.cols(idx) + u.cols(idx) * arma::diagmat(pen(idx));
    ddu = ddeta * arma::square(v.cols(idx)) + arma::ones(n, m) * arma::diagmat(pen(idx)) + this->damping;
    u.cols(idx) = u.cols(idx) - this->stepsize * (du / ddu);
}

void Newton::parallel_update (
    arma::mat & u, const arma::mat & v, 
    const arma::vec & pen, const arma::uvec & idx,
    const arma::mat & deta, const arma::mat & ddeta
) {
    // Block update via parallel operations over matrix slices
    const unsigned int n = u.n_rows;
    const unsigned int m = idx.n_rows;
    // A convenient parallelization strategy is taken depending on the row and column dimensions
    if (n >= m) {
        // If n >= m, we perform the computations row-wise
        #pragma omp parallel for
        for (unsigned int i = 0; i < n; i++) {
            arma::uvec ii = {i};
            arma::rowvec dui = - deta.row(i) * v.cols(idx) + u(ii,idx) % pen(idx).t();
            arma::rowvec ddui = ddeta.row(i) * arma::square(v.cols(idx)) + pen(idx).t() + this->damping;
            u(ii,idx) = u(ii,idx) - this->stepsize * (dui / ddui);
        }
    } else {
        // If n < m, we perform the computations column-wise
        #pragma omp parallel for
        for (const unsigned int & j : idx) {
            arma::vec duj = - deta * v.col(j) + u.col(j) * pen(j);
            arma::vec dduj = ddeta * arma::square(v.col(j)) + pen(j) + this->damping;
            u.col(j) = u.col(j) - this->stepsize * (duj / dduj);
        }
    }
}

void Newton::update (
    arma::mat & u, const arma::mat & v, 
    const arma::vec & pen, const arma::uvec & idx,
    const arma::mat & deta, const arma::mat & ddeta
) {
    if (this->parallel) {
        // Block update via parallel operations over matrix slices     
        this->parallel_update(u, v, pen, idx, deta, ddeta);
    } else {
        // Block update via elementwise matrix operations
        this->blocked_update(u, v, pen, idx, deta, ddeta);
    }
}

Rcpp::List Newton::fit (
    arma::mat & Y, 
    const arma::mat & X, const arma::mat & B, 
    const arma::mat & A, const arma::mat & Z,
    const arma::mat & U, const arma::mat & V,
    const std::unique_ptr<Family> & family,
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
    const double nm = n * m;
    const double df = n * m - n * (q + d) - m * (p + d);

    // Get the number of cores and threads
    const unsigned int mcores = this->nthreads;
    const unsigned int ncores = std::thread::hardware_concurrency();
    const unsigned int threads = parallel ? std::min(mcores, ncores-1) : 1;

    // Set the number of threads for openMP
    omp_set_num_threads(threads);

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

    // Get the predicted values, the variances and the mu-differentials
    arma::mat eta(n,m), mu(n,m), var(n,m), mueta(n,m);
    eta = get_eta(u, v, etalo, etaup);
    mu = family->linkinv(eta);
    var = family->variance(mu);
    mueta = family->mueta(eta);

    // Fill the missing values with the initial predictions
    Y.elem(isna) = mu.elem(isna);

    // Get the initial dispersion parameter
    double phi = 1;
    this->init_phi(phi, df, Y, mu, var, family);

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
        std::printf("--------------------------------------------\n");
        std::printf(" Iteration    Deviance    Change   Exe-Time \n");
        print_state(0, dev / nm, 1., time);
    }

    // Optimization loop
    int iter = 0;
    for (iter = 1; iter < this->maxiter; iter++) {

        // Fill the missing values with the current predictions
        if (anyna && iter % this->nafill == 0) {
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

        // Update the predictions, the variances and the mu-differentials
        eta = get_eta(u, v, etalo, etaup);
        mu = family->linkinv(eta);
        var = family->variance(mu);
        mueta = family->mueta(eta);

        // Update the dispersion parameter
        if (iter % 10 == 0){
            this->update_phi(phi, df, Y, mu, var, family);
        }
        
        // Update the initial deviance, penalty and objective function
        dev = arma::accu(deviance(Y, mu, family));
        pen = penalty(u, penu) + penalty(v, penv);
        objt = obj; obj = dev + 0.5 * pen;
        change = std::abs(obj - objt) / (std::abs(objt) + 1e-04);

        // Get the current execution time
        end = clock();
        time = exetime(start, end);
        
        // Store the optimization state at the current iteration
        state = arma::vec{double(iter), dev, pen, obj, change, time};
        trace.row(iter) = state.t();

        if (this->verbose && iter % frequency == 0) {
            print_state(iter, dev / nm, change, time);
        }
        
        // Check for convergence
        if (change < this->tol) {break;}
    }

    // Update the dispersion parameter
    this->update_phi(phi, df, Y, mu, var, family);

    // Get the final execution time
    end = clock();
    time = exetime(start, end);

    if (this->verbose) {
        print_state(iter, dev / nm, change, time);
        std::printf("--------------------------------------------\n");
    }
    
    // Get the final output
    Rcpp::List output;
    output["method"] = std::string("Newton");
    output["family"] = family->getfamily();
    output["link"] = family->getlink();
    output["idu"] = idu;
    output["idv"] = idv;
    output["U"] = u;
    output["V"] = v;
    output["eta"] = eta;
    output["mu"] = mu;
    output["var"] = var;
    output["phi"] = phi;
    output["penalty"] = pen;
    output["deviance"] = dev;
    output["objective"] = obj;
    output["exe.time"] = time;
    output["trace"] = trace.rows(0, iter-1);

    // Return the estimated model
    return output;
}
