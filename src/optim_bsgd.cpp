// optim_bsgd.cpp
// author: Cristian Castiglione
// creation: 07/10/2023
// last change: 13/11/2023

#include "optim.h"

using namespace glm;

void BSGD::summary () {
    std::printf("------------------\n");
    std::printf(" maxiter = %i \n", this->maxiter);
    std::printf(" eps = %.5f \n", this->eps);
    std::printf(" nafill = %i \n", this->nafill);
    std::printf(" tol = %.5f \n", this->tol);
    std::printf(" size1 = %i \n", this->size1);
    std::printf(" size2 = %i \n", this->size2);
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

void BSGD::update_rate (double & rate, const int & iter) {
    rate = this->rate0 / std::pow(1 + this->decay * this->rate0 * iter, .75);
}

void BSGD::update_deta (
    dEta & deta, const arma::uvec & idx, const arma::uvec & idy, 
    const arma::mat & Y,  const arma::mat & eta, const arma::mat & mu, 
    const std::unique_ptr<Family> & family
) {
    arma::mat var = family->variance(mu);
    arma::mat mueta = family->mueta(eta);
    deta.deta(idx, idy) = (Y - mu) % mueta / var;
    deta.ddeta(idx, idy) = (mueta % mueta) / var;
}

void BSGD::update_dpar (
    dPar & dpar, const dEta & deta, 
    const arma::uvec & idx, const arma::uvec & idy,
    const arma::mat & u, const arma::mat & v, 
    const arma::vec & penalty, const double & scale,
    const bool & transp
) {
    if (transp) {
        // This update is for the loading matrix
        const int m = idy.n_elem;
        const int d = u.n_cols;
        arma::mat grad(m,d), hess(m,d), pen(d,d);
        pen = arma::diagmat(penalty);
        grad = - scale * (deta.deta(idx, idy).t() * v.rows(idx)) + u.rows(idy) * pen;
        hess = scale * (deta.ddeta(idx, idy).t() * arma::square(v.rows(idx))) + arma::ones(m, d) * pen + this->damping;
        dpar.dpar.rows(idy) = (1 - this->rate1) * dpar.dpar.rows(idy) + this->rate1 * grad;
        dpar.ddpar.rows(idy) = (1 - this->rate2) * dpar.ddpar.rows(idy) + this->rate2 * hess;
    } else {
        // This update is for the factor scores
        const int n = idx.n_elem;
        const int d = u.n_cols;
        arma::mat grad(n,d), hess(n,d), pen(d,d);
        pen = arma::diagmat(penalty);
        grad = - scale * (deta.deta(idx, idy) * v.rows(idy)) + u.rows(idx) * pen;
        hess = scale * (deta.ddeta(idx, idy) * arma::square(v.rows(idy))) + arma::ones(n, d) * pen + this->damping;
        dpar.dpar.rows(idx) = (1 - this->rate1) * dpar.dpar.rows(idx) + this->rate1 * grad;
        dpar.ddpar.rows(idx) = (1 - this->rate2) * dpar.ddpar.rows(idx) + this->rate2 * hess;
    }
}

void BSGD::update_par (
    arma::mat & par, const dPar & dpar, const double & rate,
    const arma::uvec & idx, const arma::uvec & idy
) {
    par(idx, idy) = par(idx, idy) - rate * (dpar.dpar.rows(idx) / dpar.ddpar.rows(idx));
}

void BSGD::smooth_par (
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

// Initialize the dispersion parameter estimate
void BSGD::init_phi (
    double & phi, const int & df, 
    const arma::mat & Y, const arma::mat & mu, 
    const std::unique_ptr<Family> & family
) {
    double ssq;
    arma::mat var;
    if (family->estdisp()) {
        if (family->getfamily() == "NegativeBinomial") {
            ssq = arma::accu(arma::square(Y - mu) - arma::accu(mu));
            phi = std::max(1e-08, ssq / arma::accu(mu % mu));
            family->setdisp(1 / phi);
        } else {
            var = family->variance(mu);
            ssq = arma::accu(arma::square(Y - mu) / var);
            phi = std::max(1e-08, ssq / df);
            family->setdisp(phi);
        }
    }
}

// Update and smooth the dispersion parameter estimate
void BSGD::update_phi (
    double & phi, const double & rate, 
    const int & nm, const int & df, 
    const arma::mat & Y, const arma::mat & mu, 
    const arma::uvec & idx, const arma::uvec & idy, 
    const std::unique_ptr<Family> & family
) {
    const int ni = idx.n_elem;
    const int mi = idy.n_elem;
    const int nmi = ni * mi;
    double ssq, phit;
    arma::mat yi = Y(idx, idy);
    arma::mat mui = mu(idx, idy);
    arma::mat vari(ni, mi);
    if (family->estdisp()) {
        if (family->getfamily() == "NegativeBinomial") {
            ssq = arma::accu(arma::square(yi - mui) - mui);
            phit = ssq / arma::accu(mui % mui);
            phit = std::max(1e-08, phit);
            phi = (1 - rate) * phi + rate * phit;
            family->setdisp(1 / phi);
        } else {
            vari = family->variance(mui);
            ssq = arma::accu(arma::square(yi - mui) / vari) / nmi;
            phit = std::max(1e-08, ssq * (nm / df));
            phi = (1 - rate) * phi + rate * phit;
            family->setdisp(phi);
        }
    }
}

Rcpp::List BSGD::fit (
    arma::mat & Y, 
    const arma::mat & X, const arma::mat & B, 
    const arma::mat & A, const arma::mat & Z,
    const arma::mat & U, const arma::mat & V,
    const std::unique_ptr<Family> & family,
    const int & ncomp, const arma::vec & lambda
) {
    /*
    In this implementation we divide the dataset in chunks subsampling both the 
    rows and the columns, obtainig a rectangular tassellization of the original 
    data matrix. At each iteration of the algorthm, we pick a particular chunk
    and we update only the rows of u (score matrix) and the columns of v (loading matrix) 
    corresponding to that chunk
    This may require many iterations of the algorithm to complete one epoch, where
    one epoch correspond to a complete scan of the data matrix.
    */


    // Get the initial CPU time
    clock_t start = clock();

    // Get the data dimensions
    const int n = Y.n_rows; 
    const int m = Y.n_cols; 
    const int d = ncomp; 
    const int p = X.n_cols; 
    const int q = Z.n_cols;
    const double nm = n * m;

    // Get the range of the data, and the lower and upper bounds
    double mulo, muup, etalo, etaup;
    set_data_bounds(mulo, muup, etalo, etaup, this->eps, Y.min(), Y.max(), family);

    // Get the row and column minibatch partition of the data
    this->size1 = std::min(n, this->size1);
    this->size2 = std::min(m, this->size2);
    Chunks rowchunks(n, this->size1, true);
    Chunks colchunks(m, this->size2, true);

    // Get the row and column chunk piles, which permits us to
    // efficiently chose the new chunk to use at the next iteration 
    // of the algorithm. The new minibatch is sampled in such a way
    // that the same chunk is re-visited only when all the other 
    // chunks have already been used
    ChunkPile rowpile(rowchunks.nchunks, true);
    ChunkPile colpile(colchunks.nchunks, true);

    // Set the chunk dimensions, scale factors and indices
    int nc, mc;
    double scaler, scalec;
    arma::uvec idr, idc;

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

    // Save the optimization history
    arma::vec state(6);
    arma::mat trace(0, 6);

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
    trace = arma::join_cols(trace, state.t());

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

        // Sample the minibatch indices
        rowpile.update();
        colpile.update();
        idr = rowchunks.get_chunk(rowpile.idx);
        idc = colchunks.get_chunk(colpile.idx);
        nc = idr.n_elem;
        mc = idc.n_elem;

        // Get the minibatch normalization factors
        scaler = n / nc;
        scalec = m / mc;

        // Update the linear predictor and the mean matrix
        // arma::mat etat = ut.rows(idr) * vt.rows(idc).t();
        arma::mat etat = get_eta(ut.rows(idr), vt.rows(idc), etalo, etaup);
        arma::mat mut = family->linkinv(etat);
        arma::mat Yt = Y(idr, idc);

        eta(idr, idc) = etat;
        mu(idr, idc) = mut;

        // Update the log-likelihood differentials
        this->update_deta(deta, idr, idc, Yt, etat, mut, family);
        this->update_dpar(du, deta, idr, idc, ut.cols(idu), vt.cols(idu), penu(idu), scalec, false);
        this->update_dpar(dv, deta, idr, idc, vt.cols(idv), ut.cols(idv), penv(idv), scaler, true);
        
        // Update the parameter estimates
        this->update_par(ut, du, rate, idr, idu);
        this->update_par(vt, dv, rate, idc, idv);
        
        // Smooth the parameter estimates
        this->smooth_par(u, ut, iter, idr, idu);
        this->smooth_par(v, vt, iter, idc, idv);
        
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
            trace = arma::join_cols(trace, state.t());

            if (this->verbose) {
                print_state(iter, dev / nm, change, time);
            }
        }
        
        // Check for convergence
        if (change < this->tol) {break;}
    }

    // Get the estimated predictions
    eta = get_eta(u, v, etalo, etaup);
    mu = family->linkinv(eta);

    // The the estimated variances, deviance and penalty 
    var = family->variance(mu);
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
    output["family"] = family->getfamily();
    output["link"] = family->getlink();
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


Rcpp::List BSGD::fit2 (
    arma::mat & Y, 
    const arma::mat & X, const arma::mat & B, 
    const arma::mat & A, const arma::mat & Z,
    const arma::mat & U, const arma::mat & V,
    const std::unique_ptr<Family> & family,
    const int & ncomp, const arma::vec & lambda
) {
    /*

    In this implementation we divide the dataset in chunks subsampling both the 
    rows and the columns, obtainig a rectangular tassellization of the original 
    data matrix. At each iteration of the algorthm, we cycle over all the column
    chunks to update entirely v (loading matrix), while we update u (score matrix) 
    only for the rows corresponding to the selected row-chunk.
    Then, in this implementation, we need less iterations to reach the convergence.
    */

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

    // Get the row and column minibatch partition of the data
    this->size1 = std::min(n, this->size1);
    this->size2 = std::min(m, this->size2);
    Chunks rowchunks(n, this->size1, true);
    Chunks colchunks(m, this->size2, true);

    // Set the number of epochs needed to perform a comple column scan
    int nepochs = colchunks.nchunks;

    // Set the chunk dimensions, scale factors and indices
    int nc, mc;
    double scaler, scalec;
    arma::uvec idr, idc;

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

    // Save the optimization history
    arma::vec state(6);
    arma::mat trace(0, 6);

    // Check if and where there are some NA values
    bool anyna = !Y.is_finite();
    arma::uvec isna = arma::find_nonfinite(Y);

    // Get the linear predictor, the mean and the variance matrices
    arma::mat eta(n, m), mu(n, m), var(n, m);
    eta = get_eta(u, v, etalo, etaup);
    mu = family->linkinv(eta);

    // Fill the missing values with the initial predictions
    Y.elem(isna) = mu.elem(isna);

    // Get the initial dispersion parameter
    double phi = 1;
    this->init_phi(phi, df, Y, mu, family);

    // Get the initial deviance, penalty and objective function
    double dev, pen, obj;
    double devt, objt;
    double change, scanned;
    dev = arma::accu(deviance(Y, mu, family));
    pen = penalty(u, penu) + penalty(v, penv);
    obj = dev + 0.5 * pen;
    devt = dev;
    objt = obj;
    change = INFINITY;
    scanned = 0;

    // Get the current execution time
    clock_t end = clock();
    double time = exetime(start, end);

    // Save the objective status before starting the optimization loop
    state = arma::vec{.0, dev, pen, obj, change, time};
    trace = arma::join_cols(trace, state.t());

    // Print the optimization state
    if (verbose) {
        std::printf("------------------------------------------------------\n");
        std::printf(" Iteration    Deviance    Change   Scanned   Exe-Time \n");
        print_state(0, dev / nm, 1., time, scanned);
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

        // Set the row-chunk
        idr = rowchunks.get_chunk(iter);
        nc = idr.n_elem;
        scaler = n / nc;

        // Cycle over the column minibatch for a complete update of v 
        for (int epoch = 0; epoch < nepochs; epoch++) {

            // Set the column-chunk
            idc = colchunks.get_chunk(epoch);
            mc = idc.n_elem;
            scalec = m / mc;

            // Update the linear predictor and the mean matrix
            arma::mat etat = get_eta(ut.rows(idr), vt.rows(idc), etalo, etaup);
            arma::mat mut = family->linkinv(etat);
            arma::mat Yt = Y(idr, idc);

            eta(idr, idc) = etat;
            mu(idr, idc) = mut;

            // Update the log-likelihood differentials
            this->update_deta(deta, idr, idc, Yt, etat, mut, family);
            this->update_dpar(du, deta, idr, idc, ut.cols(idu), vt.cols(idu), penu(idu), scalec, false);
            this->update_dpar(dv, deta, idr, idc, vt.cols(idv), ut.cols(idv), penv(idv), scaler, true);
            
            // Update the parameter estimates
            this->update_par(ut, du, rate, idr, idu);
            this->update_par(vt, dv, rate, idc, idv);
            
            // Smooth the parameter estimates
            this->smooth_par(u, ut, iter, idr, idu);
            this->smooth_par(v, vt, iter, idc, idv);

            // Update the dispersion estimate
            this->update_phi(phi, rate, nm, df, Y, mu, idr, idc, family);
        }

        if (iter % frequency == 0) {
            // Update the deviance, penalty and objective function
            dev = arma::accu(deviance(Y, mu, family));
            pen = penalty(ut, penu) + penalty(vt, penv);
            // objt = obj; 
            obj = dev + 0.5 * pen;
            change = std::abs(obj - objt) / (std::abs(objt) + 1e-04);
            scanned = (iter * this->size1) / n;

            // Save the current values of the daviance, penalty and objective function
            devt = dev;
            objt = obj;

            // Get the current execution time
            end = clock();
            time = exetime(start, end);

            // Store the optimization state at the current iteration
            state = arma::vec{double(iter), dev, pen, obj, change, time};
            trace = arma::join_cols(trace, state.t());

            if (this->verbose) {
                print_state(iter, dev / nm, change, time, scanned);
            }
        }
        
        // Check for convergence
        if (change < this->tol) {break;}
    }

    // Save the final deviance, penalty and objective function
    devt = dev;
    objt = obj;

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
        scanned = (this->maxiter * this->size1) / n;
        print_state(iter, devt / nm, change, time, scanned);
        std::printf("------------------------------------------------------\n");
    }
    
    // Get the final output
    Rcpp::List output;
    output["method"] = std::string("B-SGD");
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
    output["trace"] = trace;

    // Return the estimated model
    return output;
}
