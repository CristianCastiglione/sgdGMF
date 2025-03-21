// optim_newton.cpp
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 19/10/2023

#include "optim.h"

using namespace glm;

void Newton::summary () {
    Rprintf("------------------\n");
    Rprintf(" maxiter = %i \n", this->maxiter);
    Rprintf(" stepsize = %.4f \n", this->stepsize);
    Rprintf(" eps = %.4f \n", this->eps);
    Rprintf(" nafill = %i \n", this->nafill);
    Rprintf(" tol = %.5f \n", this->tol);
    Rprintf(" damping = %.5f \n", this->damping);
    Rprintf(" parallel = %s \n", this->parallel ? "true" : "false");
    Rprintf(" verbose = %s \n", this->verbose ? "true" : "false");
    Rprintf(" frequency = %i \n", this->frequency);
    Rprintf("------------------\n");
}

void Newton::init_phi (
    double & phi, const int & df, 
    const arma::mat & Y, const arma::mat & weights, 
    const arma::mat & mu, const arma::mat & var, 
    const std::unique_ptr<Family> & family
) {
    double ssq, sm, ssm;
    if (family->estdisp()) {
        if (family->getfamily() == "NegativeBinomial") {
            // ssq = arma::accu(arma::square(Y - mu) - arma::accu(mu));
            ssq = arma::accu(weights % arma::square(Y - mu));
            sm = arma::accu(weights % mu);
            ssm = arma::accu(weights % mu % mu);
            phi = std::max(1e-08, (ssq - sm) / ssm);
            family->setdisp(1 / phi);
        } else {
            ssq = arma::accu(arma::square(Y - mu) % weights / var);
            phi = std::max(1e-08, ssq / df);
            family->setdisp(phi);
        }
    }
}

void Newton::update_phi (
    double & phi, const int & df, 
    const arma::mat & Y, const arma::mat & weights,
    const arma::mat & mu, const arma::mat & var, 
    const std::unique_ptr<Family> & family
) {
    double ssq, sm, ssm, phit, lambda = 0.5;
    if (family->estdisp()) {
        if (family->getfamily() == "NegativeBinomial") {
            ssq = arma::accu(weights % arma::square(Y - mu));
            sm = arma::accu(weights % mu);
            ssm = arma::accu(weights % mu % mu);
            phit = (ssq - sm) / ssm;
            phit = std::max(1e-08, phit);
            phi = (1 - lambda) * phi + lambda * phit;
            family->setdisp(1 / phi);
        } else {
            ssq = arma::accu(arma::square(Y - mu) % weights / var);
            phit = std::max(1e-08, ssq / df);
            phi = (1 - lambda) * phi + lambda * phit;
            family->setdisp(phi);
        }
    }
}

void Newton::update_dstat (
    dStat & dstat, 
    const arma::mat & Y, const arma::mat & offset, 
    const arma::mat & u, const arma::mat & v, 
    const double & lo, const double & up,
    const std::unique_ptr<Family> & family
) {
    if (this->parallel) {
        const unsigned int n = dstat.eta.n_rows;
        const unsigned int m = dstat.eta.n_cols;
        if (n > m) {
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (unsigned int i = 0; i < n; i++) {
                dstat.eta.row(i) = get_eta(offset.row(i), u.row(i), v, lo, up);
                dstat.mu.row(i) = family->linkinv(dstat.eta.row(i));
                dstat.var.row(i) = family->variance(dstat.mu.row(i));
                dstat.mueta.row(i) = family->mueta(dstat.eta.row(i));
                dstat.dev.row(i) = family->devresid(Y.row(i), dstat.mu.row(i));
            }
        } else {
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (unsigned int j = 0; j < m; j++) {
                dstat.eta.col(j) = get_eta(offset.col(j), u, v.row(j), lo, up);
                dstat.mu.col(j) = family->linkinv(dstat.eta.col(j));
                dstat.var.col(j) = family->variance(dstat.mu.col(j));
                dstat.mueta.col(j) = family->mueta(dstat.eta.col(j));
                dstat.dev.col(j) = family->devresid(Y.col(j), dstat.mu.col(j));
            }
        }
    } else {
        dstat.eta = get_eta(offset, u, v, lo, up);
        dstat.mu = family->linkinv(dstat.eta);
        dstat.var = family->variance(dstat.mu);
        dstat.mueta = family->mueta(dstat.eta);
        dstat.dev = family->devresid(Y, dstat.mu);
    }
}

void Newton::update_deta (
    dEta & deta, const dStat & dstat, 
    const arma::mat & Y, const arma::mat & weights,
    const std::unique_ptr<Family> & family
) {
    if (this->parallel) {
        const unsigned int n = deta.deta.n_rows;
        const unsigned int m = deta.deta.n_cols;
        if (n > m) {
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (unsigned int i = 0; i < n; i++) {
                deta.deta.row(i) = weights.row(i) % (Y.row(i) - dstat.mu.row(i)) % dstat.mueta.row(i) / dstat.var.row(i);
                deta.ddeta.row(i) = weights.row(i) % arma::square(dstat.mueta.row(i)) / dstat.var.row(i);
            }
        } else {
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (unsigned int j = 0; j < m; j++) {
                deta.deta.col(j) = weights.col(j) % (Y.col(j) - dstat.mu.col(j)) % dstat.mueta.col(j) / dstat.var.col(j);
                deta.ddeta.col(j) = weights.col(j) % arma::square(dstat.mueta.col(j)) / dstat.var.col(j);
            }
        }
    } else {
        deta.deta = weights % (Y - dstat.mu) % dstat.mueta / dstat.var;
        deta.ddeta = weights % (dstat.mueta % dstat.mueta) / dstat.var;
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
    #ifdef _OPENMP
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
    #else 
    this->blocked_update(u, v, pen, idx, deta, ddeta);
    #endif
}

void Newton::update_par (
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
    const arma::mat & O, const arma::mat & W,
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
    dStat dstat(n, m);
    dEta deta(n, m);
    
    this->update_dstat(dstat, Y, O, u, v, etalo, etaup, family);

    // Fill the missing values with the initial predictions
    Y.elem(isna) = dstat.mu.elem(isna);

    // Get the initial dispersion parameter
    double phi = 1;
    this->init_phi(phi, df, Y, W, dstat.mu, dstat.var, family);

    // Get the initial deviance, penalty and objective function
    double dev, pen, obj, objt, change;
    dev = arma::accu(W % family->devresid(Y, dstat.mu));
    pen = penalty(u, penu) + penalty(v, penv);
    obj = dev + 0.5 * pen; objt = obj;
    change = INFINITY;

    // Get the current execution time
    clock_t end = clock();
    double time = exetime(start, end);

    // Save the objective status before starting the optimization loop
    state = arma::vec{.0, dev, pen, obj, change, time};
    trace.row(0) = state.t();

    // Print the optimization state
    if (verbose) {
        Rprintf("--------------------------------------------\n");
        Rprintf(" Iteration    Deviance    Change   Exe-Time \n");
        print_state(0, dev / nm, 1., time);
    }

    // Optimization loop
    int iter = 0;
    for (iter = 1; iter < this->maxiter; iter++) {

        // Fill the missing values with the current predictions
        if (anyna && iter % this->nafill == 0) {
            Y.elem(isna) = dstat.mu.elem(isna);
        }

        // Set up helper matrices for computing differentials
        this->update_deta(deta, dstat, Y, W, family);

        // Update U and V elementwise via quasi-Newton
        this->update_par(ut, v, penu, idu, deta.deta, deta.ddeta);
        this->update_par(vt, u, penv, idv, deta.deta.t(), deta.ddeta.t());
        u = ut;
        v = vt;

        // Update the predictions, the variances and the mu-differentials
        this->update_dstat(dstat, Y, O, u, v, etalo, etaup, family);

        // Update the dispersion parameter
        if (iter % 10 == 0){
            this->update_phi(phi, df, Y, W, dstat.mu, dstat.var, family);
        }
        
        // Update the initial deviance, penalty and objective function
        dev = arma::accu(W % dstat.dev);
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
    this->update_phi(phi, df, Y, W, dstat.mu, dstat.var, family);

    // Get the final execution time
    end = clock();
    time = exetime(start, end);

    if (this->verbose) {
        print_state(iter, dev / nm, change, time);
        Rprintf("--------------------------------------------\n");
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
    output["eta"] = dstat.eta;
    output["mu"] = dstat.mu;
    output["var"] = dstat.var;
    output["phi"] = phi;
    output["penalty"] = pen;
    output["deviance"] = dev;
    output["objective"] = obj;
    output["exe.time"] = time;
    output["trace"] = trace.rows(0, iter-1);

    // Return the estimated model
    return output;
}
