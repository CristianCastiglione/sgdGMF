// deviance.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 28/09/2023

#include "deviance.h"

// Pointwise deviance
template<class F> void deviance (const arma::mat & mu, const amra::mat & y, F family);

template<class F> arma::mat deviance (const arma::mat & mu, const amra::mat & y, F family) {
    arma::mat dev(arma::size(y));
    bool anyna = !y.is_finite();
    if (anyna) {
        arma::uvec notna = arma::find_finite(y);
        dev.elem(notna) = family->devresid(y.elem(notna), mu.elem(notna));
    } else {
        dev = family->devresid(y, mu);
    }
    return dev;
}

