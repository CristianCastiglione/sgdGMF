// airwls.cpp
// author: Cristian Castiglione
// creation: 02/10/2023
// last change: 05/10/2023

#include "airwls.h"

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
    const std::unique_ptr<Family::Family> & family, 
    const arma::vec & offset, const arma::vec & penalty
) {
    // Set the tollerance threshold for the weight vector
    const double thr = 1e+06;

    // Compute the approximate suffucient statistics
    arma::vec eta = offset + X * beta;
    arma::vec mu = family->linkinv(eta);
    arma::vec mueta = family->mueta(eta);
    arma::vec varmu = family->variance(mu);
    arma::vec w = (mueta % mueta) / varmu;
    arma::vec z = (eta - offset) + (y - mu) / mueta;

    // Truncate extremely low and high weight values
    utils::trim(w, 1/thr, thr);

    // Keep only the obsevations passing some safety controls on varmu and etamu
    // arma::uvec keep = arma::find((varmu > 0) && (mueta != 0) && (w > 0));

    // Compute the new estimate via penalized WLS
    arma::mat xtwx = X.t() * arma::diagmat(w) * X + arma::diagmat(penalty);
    arma::vec xtwz = X.t() * arma::diagmat(w) * z;
    arma::vec betat = arma::solve(xtwx, xtwz);

    // Smooth the updated coefficient vector with the previous guess 
    beta = (1 - this->stepsize) * beta + this->stepsize * betat;
}

void AIRWLS::glmfit (
    arma::vec & beta, const arma::vec & y, const arma::mat & X,
    const std::unique_ptr<Family::Family> & family, 
    const arma::vec & offset, const arma::vec & penalty
) {
    // Set the convergence tolerance 
    const double tol = 1e-05;
    const int p = beta.n_rows; 
    arma::vec betaold(p);

    // Optimization cycle with Fisher scroring steps
    for (int iter = 0; iter < this->nsteps; iter++) {
        betaold = beta;
        this->glmstep(beta, y, X, family, offset, penalty);
        if (utils::absmax(beta, betaold) < tol) {break;}
    }
}

void AIRWLS::update (
    arma::mat & beta, const arma::mat & Y, const arma::mat & X,
    const std::unique_ptr<Family::Family> & family,
    const arma::uvec & idx, const arma::mat & offset, 
    const arma::vec & penalty, const bool & transp
) {
    // Set the number of slices and the number of coefficients
    unsigned int nslices = transp ? Y.n_rows : Y.n_cols;
    unsigned int ncoefs = idx.n_elem;
    arma::vec coef(ncoefs);
    arma::uvec ids(1);

    unsigned int ncores = std::thread::hardware_concurrency();
    unsigned int nthreads = parallel ? ncores-1 : 1;
    omp_set_num_threads(nthreads);
        
    // Check whether we need to optimize row- or column-wise
    if (transp) {
        // We need to transpose all the matrices to update U
        // This loop can be parrallelized with openMP
        #pragma omp parallel for
        for (unsigned int slice = 0; slice < nslices; slice++) {
            ids = {slice};
            coef = beta(ids, idx).t();
            this->glmfit(coef, Y.row(slice).t(), X.cols(idx), family, offset.row(slice).t(), penalty(idx));
            beta(ids, idx) = coef.t();
        }
    } else {
        // We don't need to transpose anything to update V
        // This loop can be parrallelized with openMP
        #pragma omp parallel for
        for (unsigned int slice = 0; slice < nslices; slice++) {
            ids = {slice};
            coef = beta(ids, idx).t();
            this->glmfit(coef, Y.col(slice), X.cols(idx), family, offset.col(slice), penalty(idx));
            beta(ids, idx) = coef.t();
        }
    }
}


Rcpp::List AIRWLS::fit (
    arma::mat & Y, // to fill NA values we need Y to be non-const
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
    arma::mat eta(n,m), mu(n,m), offset(n,m);
    eta = u * v.t();
    mu = family->linkinv(eta);
    
    // Truncate all the extreme values
    utils::trim(mu, mulo, muup);
    utils::trim(eta, etalo, etaup);

    // Fill the missing values with the initial predictions
    Y.elem(isna) = mu.elem(isna);

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

        // Update V slicewise via PIRLS
        offset = u.cols(p, p+q-1) * v.cols(p, p+q-1).t(); // = A*Zt
        this->update(v, Y, u, family, idv, offset, penv, false);

        // Update U slicewise via PIRLS
        offset = u.cols(0, p-1) * v.cols(0, p-1).t(); // = X*Bt
        this->update(u, Y, v, family, idu, offset, penu, true);

        // Update the linear predictor and the mean matrix
        eta = u * v.t();
        mu = family->linkinv(eta);

        // Truncate all the extreme values
        utils::trim(mu, mulo, muup);
        utils::trim(eta, etalo, etaup);

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
    output["method"] = std::string("AIRWLS");
    output["family"] = family->family;
    output["link"] = family->link;
    output["idu"] = idu;
    output["idv"] = idv;
    output["U"] = u;
    output["V"] = v;
    output["eta"] = eta;
    output["mu"] = mu;
    output["var"] = family->variance(mu);
    output["penalty"] = pen;
    output["deviance"] = dev;
    output["objective"] = obj;
    output["exe.time"] = time;
    output["trace"] = trace.rows(0, iter-1);

    // Return the estimated model
    return output;
}