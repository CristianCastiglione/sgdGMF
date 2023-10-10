// msgd.cpp
// author: Cristian Castiglione
// creation: 05/10/2023
// last change: 07/10/2023

#include "msgd.h"

void MSGD::summary () {
    std::printf("------------------\n");
    std::printf(" maxiter = %i \n", this->maxiter);
    std::printf(" epochs = %i \n", this->epochs);
    std::printf(" eps = %.5f \n", this->eps);
    std::printf(" nafill = %i \n", this->nafill);
    std::printf(" tol = %.5f \n", this->tol);
    std::printf(" size = %i \n", this->size);
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
    dEta & deta, const arma::uvec & idx, const arma::mat & Y, 
    const arma::mat & eta, const arma::mat & mu, 
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
        grad = - deta.deta.rows(idx).t() * v.rows(idx) + u * pen;
        hess = deta.ddeta.rows(idx).t() * arma::square(v.rows(idx)) + arma::ones(m,d) * pen + this->damping;
        dpar.dpar = (1 - this->rate1) * dpar.dpar + this->rate1 * grad;
        dpar.ddpar = (1 - this->rate2) * dpar.ddpar + this->rate2 * hess;
    } else {
        // This update is for the factor scores
        const int n = idx.n_elem;
        const int d = u.n_cols;
        arma::mat grad(n,d), hess(n,d), pen(d,d);
        pen = arma::diagmat(penalty);
        grad = - deta.deta.rows(idx) * v + u.rows(idx) * pen;
        hess = deta.ddeta.rows(idx) * (v % v) + arma::ones(n,d) * pen + this->damping;
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
    arma::mat & u, const arma::mat & ut, const int & iter
) {
    int thr = floor(double(this->maxiter) * this->burn);
    if (iter > thr) {
        double rate = 1 / (iter - thr);
        u = (1 - rate) * u + rate * ut;
    } else {
        u = ut;
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

    // Build the left and right decomposition matrices
    arma::mat u, v, ut, vt;
    u = arma::join_rows(X, A, U);
    v = arma::join_rows(B, Z, V);
    ut = u;
    vt = v;

    // Get the minibatch partition of the data
    Chunks chunks(n, this->size, true);
    const int k = chunks.nchunks;

    // Build the penalization vectors for the U and V columns
    arma::vec penu, penv;
    set_uv_penalty(penu, penv, lambda, p, q, d);

    // Get the column indices of (A,U) and (B,V)
    arma::uvec idu, idv;
    set_uv_indices(idu, idv, p, q, d);

    // Instantiate the linear predictor differentials
    dEta deta(n,m);

    // Instantiate the parameter differentials
    dPar du(n,d);
    dPar dv(m,d);
    dCube dvv(m,d,k);
    
    // Save the optimization history
    arma::vec state(6);
    arma::mat trace(maxiter, 6);

    // Check if and where there are some NA values
    bool anyna = !Y.is_finite();
    arma::uvec isna = arma::find_nonfinite(Y);

    // Get the linear predictor, the mean and the variance matrices
    arma::mat eta(n,m), mu(n,m), var(n,m);
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
    int iter; double rate = this->rate0;
    for (iter = 1; iter < this->maxiter; iter++) {

        // Fill the missing values with the current predictions
        if (iter > 0 && anyna && iter % this->nafill == 0) {
            Y.elem(isna) = mu.elem(isna);
        }

        // Update the learning rate
        this->update_rate(rate, iter);

        // Cycle over the minibatch chunks for updating u
        // Notice: this cycle can be parallelized via openMP
        for (int epoch = 0; epoch < chunks.nchunks; epoch++) {
            // Get the current minibatch chunk
            arma::uvec idx = chunks.get_chunk(epoch);

            // Get the minibatch normalizetion factor
            double scale = n / idx.n_elem;

            // Get the chunk-specific linear predictor and mean matrix 
            arma::mat etat = ut.rows(idx) * vt.t();
            arma::mat mut = family->linkinv(etat);
            arma::mat Yt = Y.rows(idx);

            utils::trim(mut, mulo, muup);
            utils::trim(etat, etalo, etaup);
            
            eta.rows(idx) = etat;
            mu.rows(idx) = mut;

            // Update the differentials wrt to eta
            this->update_deta(deta, idx, Yt, etat, mut, family);

            // Update the differentials wrt to u and v
            this->update_dpar(du, deta, idx, ut.cols(idu), vt.cols(idu), penu(idu), scale, false);
            this->update_dpar(dv, deta, idx, vt.cols(idv), ut.cols(idv), penv(idv), scale, true);

            // Accumulate the differentials wrt v
            dvv.dpar.slice(epoch) = dv.dpar;
            dvv.ddpar.slice(epoch) = dv.ddpar;

            // Update u via averaged stochastic gradient
            this->update_par(ut, du, rate, idx, idu);
            this->smooth_par(u, ut, iter);
        }

        // Update the averaged differentials wrt to v
        dv.dpar = arma::mean(dvv.dpar, 2);
        dv.ddpar = arma::mean(dvv.ddpar, 2);

        // Update v via averaged stochastic gradient
        this->update_par(vt, dv, rate, idv);
        this->smooth_par(v, vt, iter);

        if (iter % frequency == 0) {
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
    }

    // Get the final predictions, deviance and penalty
    eta = u * v.t();
    mu = family->linkinv(eta);
    var = family->variance(mu);
    dev = arma::accu(deviance(Y, mu, family));
    pen = penalty(u, penu) + penalty(v, penv);
    obj = dev + 0.5 * pen;

    // Get the final execution time
    end = clock();
    time = exetime(start, end);

    if (this->verbose) {
        print_state(iter, dev / nm, change, time);
        std::printf("-------------------------------------------\n");
    }
    
    // Get the final output
    Rcpp::List output;
    output["method"] = std::string("M-SGD");
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
    output["trace"] = trace;

    // Return the estimated model
    return output;
}




