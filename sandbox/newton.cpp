// deviance.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 29/09/2023

#include "newton.h"

void Newton::update (
    arma::mat & U, const arma::mat & V, 
    const arma::vec & pen, const arma::uvec & idx,
    const arma::mat & deta, const arma::mat & ddeta, 
    const double & stepsize, const double & damping) {

    unsigned int n = U.n_rows;
    unsigned int m = idx.n_rows;
    arma::mat dU(n,m), ddU(n,m);
    dU = - deta * V.cols(idx) + U.cols(idx) * arma::diagmat(pen(idx));
    ddU = ddeta * (V.cols(idx) % V.cols(idx)) + arma::ones(n, m) * arma::diagmat(pen(idx)) + damping;
    U.cols(idx) = U.cols(idx) - stepsize * (dU / ddU);
}

template<class F>
Rcpp::List Newton::fit (
    arma::mat Y, 
    const arma::mat & X, const arma::mat & Z,
    const arma::mat & B, const arma::mat & A,
    const arma::mat & U, const arma::mat & V,
    const F & family, const int & ncomp, 
    const double & penU, const double & penV,
    const double & penB, const double & penA,
    const int & maxiter, const double & stepsize, const double & eps,
    const int & nafill, const double & tol, const double & damping, 
    const bool & verbose, const int & frequency
) {

    // Get the initial CPU time
    clock_t start = clock();

    // Get the data dimensions
    const unsigned int n, m, nm, d, p, q;
    n = Y.n_rows; m = Y.n_cols; nm = n * m;
    d = ncomp; p = X.n_rows; q = Z.n_rows;

    // Get the range of the data, and the lower and upper bounds
    double ymin, ymax, mulo, muup, etalo, etaup;
    set_data_bounds(mulo, muup, etalo, etaup, eps, ymin, ymax, family);

    // Build the left and right decomposition matrices
    arma::mat u, v, ut, vt;
    u = arma::join_rows(X, A, U);
    v = arma::join_rows(B, Z, V);
    ut = u;
    vt = v;

    // Build the penalization vectors for the U and V columns
    arma::vec penu, penv;
    penu = join_cols(arma::zeros(p), penA * arma::ones(q), penU * arma::ones(d));
    penu = join_cols(penB * arma::zeros(p), arma::ones(q), penV * arma::ones(d));

    // Get the column indices of (A,U) and (B,V)
    arma::uvec idu, idv;
    idu = join_cols(arma::linspace<arma::uvec>(p, p+q-1, q), arma::linspace<arma::uvec>(p+q, p+q+d, d));
    idv = join_cols(arma::linspace<arma::uvec>(0, p, p), arma::linspace<arma::uvec>(p+q, p+q+d, d));

    // Save the optimization history
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
    // arma::mat dratio(n,m), ddratio(n,m), ddiff(n,m);

    // Get the current execution time
    clock_t end = clock();
    double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;

    // Save the objective status before starting the optimization loop
    trace.row(0) = arma::vec{0, dev, pen, obj, change, time};

    // Optimization loop
    int iter;
    for (iter = 1; iter < maxiter; iter++) {

        // Fill the missing values with the current predictions
        if (iter > 0 && iter % nafill == 0) {
            Y.elem(isna) = mu.elem(isna);
        }

        // Set up helper matrices for computing differentials
        deta = (Y - mu) % mueta / var;
        ddeta = (mueta % mueta) / var;
        // dratio = mueta / var;
        // ddratio = mueta % dratio;
        // ddiff = (Y - mu) % dratio;

        // Update U and V elementwise via quasi-Newton
        Newton::update(ut, v, penu, idu, deta, ddeta, stepsize, damping);
        Newton::update(vt, u, penv, idv, deta.t(), ddeta.t(), stepsize, damping);
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
        time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        
        // Store the optimization state at the current iteration
        trace.row(iter) = arma::vec{iter, dev, pen, obj, change, time};

        // Check for convergence
        if (change < tol) {break;}
    }

    // Get the final execution time
    end = clock();
    time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    
    // Get the final output
    Rcpp::List output;
    output["U"] = u;
    output["V"] = v;
    output["eta"] = eta;
    output["mu"] = mu;
    output["var"] = var;
    output["penalty"] = pen;
    output["deviance"] = dev;
    output["objective"] = obj;
    output["exe.time"] = time;
    output["trace"] = trace.rows(0, iter);
    output["method"] = std::string("Newton");

    // Return the estimated model
    return output;
}

template<class F>
void set_data_bounds (
    double & mulo, double & muup, double & etalo, double & etaup, 
    const double & eps, const double & ymin, const double & ymax, F family
) {
    // We compute the lower and upper bounds on a matrices of dim 1x1,
    // since we need to back-transform them using the family->linkfun()
    // method, which is defined only for arma::mat objects
    arma::mat mulot(1,1), muupt(1,1);
    arma::mat etalot(1,1), etaupt(1,1);
    
    // Mean boundas
    muupt(0,0) = ymax - eps * (ymax - ymin);
    mulot(0,0) = ymin + eps * (ymax - ymin);
    
    // Linear predictor boundas
    etalot = family->linkfun(mulot);
    etaupt = family->linkfun(muupt);

    // Inplace sssignment
    mulo = mulot(0,0); muup = muupt(0,0);
    etalo = etalot(0,0); etaup = etaupt(0,0);
}