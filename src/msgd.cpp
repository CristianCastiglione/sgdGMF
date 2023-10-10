// msgd.cpp
// author: Cristian Castiglione
// creation: 09/10/2023
// last change: 10/10/2023

#include "sgd.h"

void MSGD::summary () {
    std::printf("------------------\n");
    std::printf(" maxiter = %i \n", this->maxiter);
    std::printf(" eps = %.5f \n", this->eps);
    std::printf(" nafill = %i \n", this->nafill);
    std::printf(" tol = %.5f \n", this->tol);
    std::printf(" size1 = %i \n", this->size);
    std::printf(" burn = %.5f \n", this->burn);
    std::printf(" rate0 = %.5f \n", this->rate0);
    std::printf(" decay = %.5f \n", this->decay);
    std::printf(" damping = %.5f \n", this->damping);
    std::printf(" rate1 = %.5f \n", this->rate1);
    std::printf(" rate2 = %.5f \n", this->rate2);
    std::printf(" parallel = %s \n", this->parallel ? "true" : "false");
    std::printf(" verbose = %s \n", this->verbose ? "true" : "false");
    std::printf(" frequency = %i \n", this->frequency);
    std::printf(" progress = %s \n", this->progress ? "true" : "false");
    std::printf("------------------\n");
}

void MSGD::update_rate (double & rate, const int & iter) {
    rate = this->rate0 / std::pow(1 + this->decay * this->rate0 * iter, .75);
}

void MSGD::update_deta (
    dEta & deta, const arma::uvec & idx,
    const arma::mat & Y,  const arma::mat & eta, const arma::mat & mu, 
    const std::unique_ptr<Family::Family> & family
) {
    arma::mat var = family->variance(mu);
    arma::mat mueta = family->mueta(eta);
    deta.deta.rows(idx) = (Y - mu) % mueta / var;
    deta.ddeta.rows(idx) = (mueta % mueta) / var;
}

void MSGD::update_dpar (
    dPar & dpar, const dEta & deta, const arma::uvec & idx, 
    const arma::mat & u, const arma::mat & v, 
    const arma::vec & penalty, const double & scale,
    const bool & transp
) {
    if (transp) {
        // This update is for the loading matrix
        const int m = u.n_rows;
        const int d = u.n_cols;
        arma::mat grad(m,d), hess(m,d), pen(d,d);
        pen = arma::diagmat(penalty);
        grad = - scale * (deta.deta.rows(idx).t() * v.rows(idx)) + u * pen;
        hess = scale * (deta.ddeta.rows(idx).t() * arma::square(v.rows(idx))) + arma::ones(m, d) * pen + this->damping;
        dpar.dpar = (1 - this->rate1) * dpar.dpar + this->rate1 * grad;
        dpar.ddpar = (1 - this->rate2) * dpar.ddpar + this->rate2 * hess;
    } else {
        // This update is for the factor scores
        const int n = idx.n_elem;
        const int d = u.n_cols;
        arma::mat grad(n,d), hess(n,d), pen(d,d);
        pen = arma::diagmat(penalty);
        grad = - scale * (deta.deta.rows(idx) * v) + u.rows(idx) * pen;
        hess = scale * (deta.ddeta.rows(idx) * arma::square(v)) + arma::ones(n, d) * pen + this->damping;
        dpar.dpar.rows(idx) = (1 - this->rate1) * dpar.dpar.rows(idx) + this->rate1 * grad;
        dpar.ddpar.rows(idx) = (1 - this->rate2) * dpar.ddpar.rows(idx) + this->rate2 * hess;
    }
}

void MSGD::update_par (
    arma::mat & par, const dPar & dpar, 
    const double & rate, const arma::uvec & idy
) {
    par.cols(idy) = par.cols(idy) - rate * (dpar.dpar / dpar.ddpar);
}

void MSGD::update_par (
    arma::mat & par, const dPar & dpar, const double & rate,
    const arma::uvec & idx, const arma::uvec & idy
) {
    par(idx, idy) = par(idx, idy) - rate * (dpar.dpar.rows(idx) / dpar.ddpar.rows(idx));
}

void MSGD::smooth_par (
    arma::mat & u, const arma::mat & ut, 
    const int & iter, const arma::uvec & idy
) {
    int thr = floor(double(this->maxiter) * this->burn);
    if (iter > thr) {
        double rate = 1 / (iter - thr);
        u.cols(idy) = (1 - rate) * u.cols(idy) + rate * ut.cols(idy);
    } else {
        u.cols(idy) = ut.cols(idy);
    }
}

void MSGD::smooth_par (
    arma::mat & u, const arma::mat & ut, const int & iter,
    const arma::uvec & idx, const arma::uvec & idy
) {
    int thr = floor(double(this->maxiter) * this->burn);
    if (iter > thr) {
        double rate = 1 / (iter - thr);
        u(idx, idy) = (1 - rate) * u(idx, idy) + rate * ut(idx, idy);
    } else {
        u(idx, idy) = ut(idx, idy);
    }
}

Rcpp::List MSGD::fit (
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

    // Get the row and column minibatch partition of the data
    this->size = std::min(n, this->size);
    Chunks chunks(n, this->size, true);
    int nepochs = chunks.nchunks;

    // Set the chunk dimensions, scale factors and indices
    int nc;
    double scale;
    arma::uvec idx;

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

    // Instantiate the differential wrt eta, u and v
    dEta deta(n, m);
    dPar du(n, q+d);
    dPar dv(m, p+d);
    dCube dV(m, p+d, nepochs);

    // Save the optimization history
    arma::vec state(6);
    arma::mat trace(maxiter, 6);

    // Check if and where there are some NA values
    bool anyna = !Y.is_finite();
    arma::uvec isna = arma::find_nonfinite(Y);

    // Get the linear predictor, the mean and the variance matrices
    arma::mat eta(n, m), mu(n, m), var(n, m);
    eta = get_eta(u, v, etalo, etaup);
    mu = family->linkinv(eta);

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
        std::printf("--------------------------------------------\n");
        std::printf(" Iteration    Deviance    Change   Exe-Time \n");
        print_state(0, dev / nm, 1., time);
    }

    // Optimization loop
    int iter = 0; 
    double rate = this->rate0;
    for (iter = 1; iter < this->maxiter; iter++) {

        // Fill the missing values with the current predictions
        if (iter > 0 && anyna && iter % this->nafill == 0) {
            Y.elem(isna) = mu.elem(isna);
        }

        // Update the learning rate
        this->update_rate(rate, iter);

        // Cycle over the epochs
        for (int epoch = 0; epoch < nepochs; epoch++) {
            // Get the minibatch indices
            idx = chunks.get_chunk(epoch);
            nc = idx.n_elem;
            scale = n / nc;

            dPar dvt(m, p+d);
            dvt.dpar = dv.dpar;
            dvt.ddpar = dv.ddpar;

            // Update the linear predictor and the mean matrix
            arma::mat etat = get_eta(ut.rows(idx), vt, etalo, etaup);
            arma::mat mut = family->linkinv(etat);
            arma::mat Yt = Y.rows(idx);

            eta.rows(idx) = etat;
            mu.rows(idx) = mut;

            // Update the log-likelihood differentials of u
            this->update_deta(deta, idx, Yt, etat, mut, family);
            this->update_dpar(du, deta, idx, ut.cols(idu), vt.cols(idu), penu(idu), scale, false);

            // Update the log-likelihood differentials of v
            this->update_dpar(dvt, deta, idx, vt.cols(idv), ut.cols(idv), penv(idv), scale, true);
            dV.dpar.slice(epoch) = dvt.dpar;
            dV.ddpar.slice(epoch) = dvt.ddpar;
            
            // Update and smooth the estimate of u
            this->update_par(ut, du, rate, idx, idu);
            this->smooth_par(u, ut, iter, idx, idu);
        }

        // Accumulate al the differentials wrt v
        dv.dpar = arma::mean(dV.dpar, 2);
        dv.ddpar = arma::mean(dV.ddpar, 2);

        // Update and smooth the estimate of v
        this->update_par(vt, dv, rate, idv);
        this->smooth_par(v, vt, iter, idv);
        
        if (iter % frequency == 0) {
            // Update the deviance, penalty and objective functions
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

            if (this->verbose) {
                print_state(iter, dev / nm, change, time);
            }
        }
        
        // Check for convergence
        if (change < this->tol) {break;}
    }

    // Get the estimated predictions and variances
    eta = get_eta(u, v, etalo, etaup);
    mu = family->linkinv(eta);
    var = family->variance(mu);

    // The the deviance, penalty and objective function     
    dev = arma::accu(deviance(Y, mu, family));
    pen = penalty(u, penu) + penalty(v, penv);
    obj = dev + 0.5 * pen;
    
    // Get the final execution time
    end = clock();
    time = exetime(start, end);

    if (this->verbose) {
        print_state(iter, dev / nm, change, time);
        std::printf("--------------------------------------------\n");
    }
    
    // Get the final output
    Rcpp::List output;
    output["method"] = std::string("B-SGD");
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


Rcpp::List MSGD::fit2 (
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

    // Get the number of cores and threads
    unsigned int ncores = std::thread::hardware_concurrency();
    unsigned int nthreads = parallel ? ncores-1 : 1;

    // Set the number of threads for openMP
    omp_set_num_threads(nthreads);

    // Get the range of the data, and the lower and upper bounds
    double mulo, muup, etalo, etaup;
    set_data_bounds(mulo, muup, etalo, etaup, this->eps, Y.min(), Y.max(), family);

    // Get the row and column minibatch partition of the data
    this->size = std::min(n, this->size);
    Chunks chunks(n, this->size, true);
    int nepochs = chunks.nchunks;

    // Set the chunk dimensions, scale factors and indices
    int nc;
    double scale;
    arma::uvec idx;

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

    // Instantiate the differential wrt eta, u and v
    dEta deta(n, m);
    dPar du(n, q+d);
    dPar dv(m, p+d);
    dCube dV(m, p+d, nepochs);

    // Save the optimization history
    arma::vec state(6);
    arma::mat trace(maxiter, 6);

    // Check if and where there are some NA values
    bool anyna = !Y.is_finite();
    arma::uvec isna = arma::find_nonfinite(Y);

    // Get the linear predictor, the mean and the variance matrices
    arma::mat eta(n, m), mu(n, m), var(n, m);
    eta = get_eta(u, v, etalo, etaup);
    mu = family->linkinv(eta);

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
        std::printf("--------------------------------------------\n");
        std::printf(" Iteration    Deviance    Change   Exe-Time \n");
        print_state(0, dev / nm, 1., time);
    }

    // Optimization loop
    int iter = 0; 
    double rate = this->rate0;
    for (iter = 1; iter < this->maxiter; iter++) {

        // Fill the missing values with the current predictions
        if (iter > 0 && anyna && iter % this->nafill == 0) {
            Y.elem(isna) = mu.elem(isna);
        }

        // Update the learning rate
        this->update_rate(rate, iter);

        // Cycle over the epochs
        #pragma omp parallel for
        for (int epoch = 0; epoch < nepochs; epoch++) {
            // Get the minibatch indices
            idx = chunks.get_chunk(epoch);
            nc = idx.n_elem;
            scale = n / nc;

            dPar dvt(m, p+d);
            dvt.dpar = dv.dpar;
            dvt.ddpar = dv.ddpar;

            // Update the linear predictor and the mean matrix
            arma::mat etat = get_eta(ut.rows(idx), vt, etalo, etaup);
            arma::mat mut = family->linkinv(etat);
            arma::mat Yt = Y.rows(idx);

            eta.rows(idx) = etat;
            mu.rows(idx) = mut;

            // Update the log-likelihood differentials of u
            this->update_deta(deta, idx, Yt, etat, mut, family);
            this->update_dpar(du, deta, idx, ut.cols(idu), vt.cols(idu), penu(idu), scale, false);

            // Update the log-likelihood differentials of v
            this->update_dpar(dvt, deta, idx, vt.cols(idv), ut.cols(idv), penv(idv), scale, true);
            dV.dpar.slice(epoch) = dvt.dpar;
            dV.ddpar.slice(epoch) = dvt.ddpar;
            
            // Update and smooth the estimate of u
            this->update_par(ut, du, rate, idx, idu);
            this->smooth_par(u, ut, iter, idx, idu);
        }

        // Accumulate al the differentials wrt v
        dv.dpar = arma::mean(dV.dpar, 2);
        dv.ddpar = arma::mean(dV.ddpar, 2);

        // Update and smooth the estimate of v
        this->update_par(vt, dv, rate, idv);
        this->smooth_par(v, vt, iter, idv);
        
        if (iter % frequency == 0) {
            // Update the deviance, penalty and objective functions
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

            if (this->verbose) {
                print_state(iter, dev / nm, change, time);
            }
        }
        
        // Check for convergence
        if (change < this->tol) {break;}
    }

    // Get the estimated predictions and variances
    eta = get_eta(u, v, etalo, etaup);
    mu = family->linkinv(eta);
    var = family->variance(mu);

    // The the deviance, penalty and objective function 
    dev = arma::accu(deviance(Y, mu, family));
    pen = penalty(u, penu) + penalty(v, penv);
    obj = dev + 0.5 * pen;
    
    // Get the final execution time
    end = clock();
    time = exetime(start, end);

    if (this->verbose) {
        print_state(iter, dev / nm, change, time);
        std::printf("--------------------------------------------\n");
    }
    
    // Get the final output
    Rcpp::List output;
    output["method"] = std::string("B-SGD");
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
