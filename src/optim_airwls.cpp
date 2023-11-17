// optim_airwls.cpp
// author: Cristian Castiglione
// creation: 02/10/2023
// last change: 10/10/2023

#include "optim.h"

using namespace glm;

void AIRWLS::summary () {
    std::printf("------------------\n");
    std::printf(" maxiter = %i \n", this->maxiter);
    std::printf(" nsteps = %i \n", this->nsteps);
    std::printf(" stepsize = %.4f \n", this->stepsize);
    std::printf(" eps = %.4f \n", this->eps);
    std::printf(" nafill = %i \n", this->nafill);
    std::printf(" tol = %.5f \n", this->tol);
    std::printf(" damping = %.5f \n", this->damping);
    std::printf(" verbose = %s \n", this->verbose ? "true" : "false");
    std::printf(" frequency = %i \n", this->frequency);
    std::printf(" parallel = %s \n", this->parallel ? "true" : "false");
    std::printf("------------------\n");
}

void AIRWLS::glmstep (
    arma::vec & beta, const arma::vec & y, const arma::mat & X,
    const std::unique_ptr<Family> & family, 
    const arma::vec & offset, const arma::vec & penalty
) {
    // Set the tollerance threshold and the parameter dimension
    const double thr = 1e+06;
    // const int p = beta.n_rows;

    // Compute the approximate suffucient statistics
    arma::vec eta = offset + X * beta;
    arma::vec mu = family->linkinv(eta);
    arma::vec mueta = family->mueta(eta);
    arma::vec varmu = family->variance(mu);
    arma::vec w = (mueta % mueta) / varmu;
    arma::vec z = (eta - offset) + (y - mu) / mueta;

    // Truncate extremely low and high weight values
    utils::trim(w, 1/thr, thr);

    // Compute the new estimate via penalized WLS
    arma::mat pen = arma::diagmat(penalty + this->damping);
    arma::mat xtwx = X.t() * arma::diagmat(w) * X;
    arma::vec xtwz = X.t() * (w % z);
    arma::vec betat = beta;
    bool status = arma::solve(betat, xtwx + pen, xtwz);
    
    if (!status) {
        betat = beta;
    }
    

    // try {
    //     // This system might fail in the case where xtwx is not positive definite
    //     betat = arma::solve(xtwx, xtwz);
    // } catch (...) {
    //     // Here we handle ill-conditioned systems by approximating the matrix inverse
    //     // using only the first positive (enough) eigenvalues/vectors
    //     arma::vec eigval(p);
    //     arma::mat eigvec(p, p);
    //     arma::eig_sym(eigval, eigvec, xtwx);
    //     arma::uvec idx = arma::find(eigval > 1/thr);
    //     betat = eigvec.cols(idx) * arma::diagmat(1 / eigval(idx)) * eigvec.cols(idx).t() * xtwz;
    // }

    // Smooth the updated coefficient vector with the previous guess 
    beta = (1 - this->stepsize) * beta + this->stepsize * betat;
}

void AIRWLS::glmfit (
    arma::vec & beta, const arma::vec & y, const arma::mat & X,
    const std::unique_ptr<Family> & family, 
    const arma::vec & offset, const arma::vec & penalty
) {
    // Set the convergence tolerance 
    const double tol = 1e-05;
    const int p = beta.n_rows;

    // Optimization cycle with Fisher scroring steps 
    arma::vec betaold(p);
    for (int iter = 0; iter < this->nsteps; iter++) {
        betaold = beta;
        this->glmstep(beta, y, X, family, offset, penalty);
        if (utils::absmax(beta, betaold) < tol) {break;}
    }
}

void AIRWLS::sequential_update (
    arma::mat & beta, const arma::mat & Y, const arma::mat & X,
    const std::unique_ptr<Family> & family,
    const arma::uvec & idx, const arma::mat & offset, 
    const arma::vec & penalty, const bool & transp
) {
    // Set the number of slices and the number of coefficients
    unsigned int nslices = transp ? Y.n_rows : Y.n_cols;
    unsigned int ncoefs = idx.n_elem;

    arma::uvec ids(1);
    arma::vec coef(ncoefs);

    // Check whether we need to optimize by row or by column
    if (transp) {
        // We need to transpose all the matrices to update U
        for (unsigned int slice = 0; slice < nslices; slice++) {
            ids = {slice};
            coef = beta(ids, idx).t();
            this->glmfit(coef, Y.row(slice).t(), X.cols(idx), family, offset.row(slice).t(), penalty(idx));
            beta(ids, idx) = coef.t();
        }
    } else {
        // We don't need to transpose anything to update V
        for (unsigned int slice = 0; slice < nslices; slice++) {
            ids = {slice};
            coef = beta(ids, idx).t();
            this->glmfit(coef, Y.col(slice), X.cols(idx), family, offset.col(slice), penalty(idx));
            beta(ids, idx) = coef.t();
        }
    }
}

void AIRWLS::parallel_update (
    arma::mat & beta, const arma::mat & Y, const arma::mat & X,
    const std::unique_ptr<Family> & family,
    const arma::uvec & idx, const arma::mat & offset, 
    const arma::vec & penalty, const bool & transp
) {
    // Set the number of slices and the number of coefficients
    unsigned int nslices = transp ? Y.n_rows : Y.n_cols;
    unsigned int ncoefs = idx.n_elem;

    // Get the number of cores and threads
    const unsigned int mcores = this->nthreads;
    unsigned int ncores = std::thread::hardware_concurrency();
    unsigned int threads = this->parallel ? std::min(mcores, ncores-1) : 1;

    // Set the number of threads for openMP
    omp_set_num_threads(threads);
            
    // Check whether we need to optimize by row or by column
    if (transp) {
        // We need to transpose all the matrices to update U
        #pragma omp parallel 
        {
            // We need to instantiate ids and coefs here in order to make 
            // openMP aware that they are private thread-specific variables
            arma::uvec ids(1);
            arma::vec coef(ncoefs);
            // This loop can be parallelized with openMP since any iteration 
            // is independent of each other
            #pragma omp for
            for (unsigned int slice = 0; slice < nslices; slice++) {
                ids = {slice};
                coef = beta(ids, idx).t();
                this->glmfit(coef, Y.row(slice).t(), X.cols(idx), family, offset.row(slice).t(), penalty(idx));
                beta(ids, idx) = coef.t();
            }
        }
    } else {
        // We don't need to transpose anything to update V
        #pragma omp parallel 
        {
            // We need to instantiate ids and coefs here in order to make 
            // openMP aware that they are private thread-specific variables
            arma::uvec ids(1);
            arma::vec coef(ncoefs);
            // This loop can be parallelized with openMP since any iteration 
            // is independent of each other
            #pragma omp for
            for (unsigned int slice = 0; slice < nslices; slice++) {
                ids = {slice};
                coef = beta(ids, idx).t();
                this->glmfit(coef, Y.col(slice), X.cols(idx), family, offset.col(slice), penalty(idx));
                beta(ids, idx) = coef.t();
            }
        } 
    }
}

void AIRWLS::update (
    arma::mat & beta, const arma::mat & Y, const arma::mat & X,
    const std::unique_ptr<Family> & family,
    const arma::uvec & idx, const arma::mat & offset, 
    const arma::vec & penalty, const bool & transp
) {
    if (this->parallel) {
        parallel_update(beta, Y, X, family, idx, offset, penalty, transp);
    } else {
        sequential_update(beta, Y, X, family, idx, offset, penalty, transp);
    }
}

void AIRWLS::init_phi (
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

void AIRWLS::update_phi (
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


Rcpp::List AIRWLS::fit (
    arma::mat & Y, // to fill NA values we need Y to be non-const
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

    // Get the range of the data, and the lower and upper bounds
    double mulo, muup, etalo, etaup;
    set_data_bounds(mulo, muup, etalo, etaup, this->eps, Y.min(), Y.max(), family);

    // Build the left and right decomposition matrices
    arma::mat u, v, ut, vt;
    u = arma::join_rows(X, A, U);
    v = arma::join_rows(B, Z, V);

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

    // Get the linear predictor and the mean matrices
    arma::mat eta(n,m), mu(n,m), var(n,m), offset(n,m);
    eta = get_eta(u, v, etalo, etaup);
    mu = family->linkinv(eta);
    var = family->variance(mu);

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

        // Update V slicewise via PIRLS
        offset = u.cols(p, p+q-1) * v.cols(p, p+q-1).t(); // = A*Zt
        this->update(v, Y, u, family, idv, offset, penv, false);

        // Update U slicewise via PIRLS
        offset = u.cols(0, p-1) * v.cols(0, p-1).t(); // = X*Bt
        this->update(u, Y, v, family, idu, offset, penu, true);

        // Update the linear predictor and the mean matrix
        eta = get_eta(u, v, etalo, etaup);
        mu = family->linkinv(eta);

        // Update the dispersion parameter
        if (iter % 10 == 0){
            var = family->variance(mu);
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
    var = family->variance(mu);
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
    output["method"] = std::string("AIRWLS");
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
