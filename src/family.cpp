// family.cpp
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 28/09/2023

#include "family.h"

// Gaussian family
arma::mat Family::Gaussian::variance (const arma::mat & mu) const {return arma::ones(arma::size(mu));}
arma::mat Family::Gaussian::initialize (const arma::mat & y) const {return y;}
bool Family::Gaussian::validmu (const arma::mat & mu) const {return true;}
bool Family::Gaussian::valideta (const arma::mat & eta) const {return true;}
arma::mat Family::Gaussian::devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const {
    return wt % arma::square(y - mu);
}

// Binomial family
arma::mat Family::Binomial::variance (const arma::mat & mu) const {return mu % (1 - mu);}
arma::mat Family::Binomial::initialize (const arma::mat & y) const {return (2 * y - 1) + 0.05 * arma::randu(size(y));}
bool Family::Binomial::validmu (const arma::mat & mu) const {return arma::all(arma::all(mu > 0)) && arma::all(arma::all(mu < 1));}
bool Family::Binomial::valideta (const arma::mat & eta) const {return true;}
arma::mat Family::Binomial::devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const {
    return - 2 * wt % (y % arma::log(mu) + (1 - y) % arma::log(1 - mu));
}

// Poisson family
arma::mat Family::Poisson::variance (const arma::mat & mu) const {return mu;}
arma::mat Family::Poisson::initialize (const arma::mat & y) const {return this->linkfun(y + 0.1);}
bool Family::Poisson::validmu (const arma::mat & mu) const {return arma::all(arma::all(mu > 0));}
bool Family::Poisson::valideta (const arma::mat & eta) const {return true;}
arma::mat Family::Poisson::devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const {
    return 2 * wt % (y % arma::log(y / mu) - (y - mu));
}

// Gamma family
arma::mat Family::Gamma::variance (const arma::mat & mu) const {return arma::square(mu);}
arma::mat Family::Gamma::initialize (const arma::mat & y) const {return this->linkfun(y);}
bool Family::Gamma::validmu (const arma::mat & mu) const {return arma::all(arma::all(mu > 0));}
bool Family::Gamma::valideta (const arma::mat & eta) const {return true;}
arma::mat Family::Gamma::devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const {
    return - 2 * wt % (arma::log(y / mu) - (y - mu) / mu);
}
