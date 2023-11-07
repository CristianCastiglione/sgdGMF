// deviance.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 28/09/2023

#include "deviance.h"

using namespace glm;

// Pointwise deviance
void deviance (
    arma::mat & dev, const arma::mat & y, const arma::mat & mu, 
    const std::unique_ptr<Family> & family
) {
    bool anyna = !y.is_finite();
    if (anyna) {
        arma::uvec notna = arma::find_finite(y);
        dev.elem(notna) = family->devresid(y.elem(notna), mu.elem(notna));
    } else {
        dev = family->devresid(y, mu);
    }
};

arma::mat deviance (
    const arma::mat & y, const arma::mat & mu, 
    const std::unique_ptr<Family> & family
) {
    arma::mat dev(arma::size(y));
    deviance(dev, y, mu, family);
    return dev;
}

// Penalty matrix
void penalty (double & pen, const arma::mat & u, const arma::vec & p) {
    pen = arma::accu((u % u) * arma::diagmat(p));
};

double penalty (const arma::mat & u, const arma::vec & p) {
    double pen;
    penalty(pen, u, p);
    return pen;
};