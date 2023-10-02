// airwls.cpp
// author: Cristian Castiglione
// creation: 02/10/2023
// last change: 02/10/2023

#include "airwls.h"

void AIRWLS::wlsfit (
    arma::vec & beta, const arma::vec & y, const arma::mat & X,
    const std::unique_ptr<Family::Family> & family, 
    const arma::vec & offset, const double & penalty
) {
    const int p = X.n_cols;
    // const double thr = 1e+10;

    arma::vec eta = offset + X * y;
    arma::vec mu = family->linkinv(eta);
    arma::vec mueta = family->mueta(eta);
    arma::vec varmu = family->variance(mu);
    arma::vec invw = varmu / (mueta % mueta);
    arma::vec z = (eta - offset) + invw % (y - mu);

    // We still have to implement a safety control for extreme gradient and weight values
    // ...

    arma::mat xtwx = X.t() * arma::diagmat(invw) * X + penalty * arma::eye(p, p);
    arma::vec xtwz = X.t() * arma::diagmat(invw) * z;
    arma::vec betat = arma::solve(xtwx, xtwz);
    beta = (1 - this->stepsize) * beta + this->stepsize * betat;
}

void AIRWLS::glmfit (
    arma::vec & beta, const arma::vec & y, const arma::mat & X,
    const std::unique_ptr<Family::Family> & family, 
    const arma::vec & offset, const double & penalty
) {
    double tol = 1e-05;
    const int p = beta.n_rows; 
    arma::vec betaold(p);
    for (int iter = 0; iter < this->nsteps; iter++) {
        betaold = beta;
        this->wlsfit(beta, y, X, family, offset, penalty);
        if (utils::absmax(beta, betaold) < tol) {break;}
    }
}

void AIRWLS::update (
    arma::mat & beta, const arma::mat & Y, const arma::mat & X,
    const std::unique_ptr<Family::Family> & family,
    const int & nslices, const arma::mat & offset, 
    const double & penalty
) {
    const int p = X.n_cols;
    arma::vec coef(p);
    for (int slice = 0; slice < nslices; slice++) {
        coef = beta.col(slice);
        this->glmfit(coef, Y.col(slice), X, family, offset.col(slice), penalty);
        beta.col(slice) = coef;
    }
}