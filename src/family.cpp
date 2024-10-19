// family.cpp
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 30/09/2023

#include "family.h"

using namespace glm;

// Gaussian family
arma::mat Gaussian::variance (const arma::mat & mu) const {return arma::ones(size(mu));}
arma::mat Gaussian::initialize (const arma::mat & y) const {return y;}
arma::mat Gaussian::devresid (const arma::mat & y, const arma::mat & mu) const {return arma::square(y - mu);}
bool Gaussian::validmu (const arma::mat & mu) const {return true;}
bool Gaussian::valideta (const arma::mat & eta) const {return true;}

// Binomial family
arma::mat Binomial::variance (const arma::mat & mu) const {return mu % (1 - mu);}
arma::mat Binomial::initialize (const arma::mat & y) const {return 2 * y - 1;}
arma::mat Binomial::devresid (const arma::mat & y, const arma::mat & mu) const {return - 2 * (y % arma::log(mu) + (1 - y) % arma::log1p(-mu));}
bool Binomial::validmu (const arma::mat & mu) const {return arma::all(arma::all(mu > 0)) && arma::all(arma::all(mu < 1));}
bool Binomial::valideta (const arma::mat & eta) const {return true;}

// Poisson family
arma::mat Poisson::variance (const arma::mat & mu) const {return mu;}
arma::mat Poisson::initialize (const arma::mat & y) const {return this->linkfun(y + 0.1);}
arma::mat Poisson::devresid (const arma::mat & y, const arma::mat & mu) const {return 2 * (utils::xlogx(y) - y % arma::log(mu) - (y - mu));}
bool Poisson::validmu (const arma::mat & mu) const {return arma::all(arma::all(mu > 0));}
bool Poisson::valideta (const arma::mat & eta) const {return true;}

// Gamma family
arma::mat Gamma::variance (const arma::mat & mu) const {return arma::square(mu);}
arma::mat Gamma::initialize (const arma::mat & y) const {return this->linkfun(y);}
arma::mat Gamma::devresid (const arma::mat & y, const arma::mat & mu) const {return - 2 * (arma::log(y / mu) - (y - mu) / mu);}
bool Gamma::validmu (const arma::mat & mu) const {return arma::all(arma::all(mu > 0));}
bool Gamma::valideta (const arma::mat & eta) const {return true;}

// Negative-Binomial family
arma::mat NegativeBinomial::variance (const arma::mat & mu) const {return mu + (mu % mu) / this->dispersion;}
arma::mat NegativeBinomial::initialize (const arma::mat & y) const {return this->linkfun(y + 0.1);}
arma::mat NegativeBinomial::devresid (const arma::mat & y, const arma::mat & mu) const {
    double phi = this->dispersion;
    return 2 * (utils::xlogx(y) - y % arma::log(mu) - (y + phi) % (arma::log(y + phi) - arma::log(mu + phi)));
}
bool NegativeBinomial::validmu (const arma::mat & mu) const {return arma::all(arma::all(mu > 0));}
bool NegativeBinomial::valideta (const arma::mat & eta) const {return true;}

// Quasi-Binomial family
arma::mat QuasiBinomial::variance (const arma::mat & mu) const {return mu % (1 - mu);}
arma::mat QuasiBinomial::initialize (const arma::mat & y) const {return 2 * y - 1;}
arma::mat QuasiBinomial::devresid (const arma::mat & y, const arma::mat & mu) const {return - 2 * (y % arma::log(mu) + (1 - y) % arma::log1p(-mu));}
bool QuasiBinomial::validmu (const arma::mat & mu) const {return arma::all(arma::all(mu > 0)) && arma::all(arma::all(mu < 1));}
bool QuasiBinomial::valideta (const arma::mat & eta) const {return true;}

// Quasi-Poisson family
arma::mat QuasiPoisson::variance (const arma::mat & mu) const {return mu;}
arma::mat QuasiPoisson::initialize (const arma::mat & y) const {return this->linkfun(y + 0.1);}
arma::mat QuasiPoisson::devresid (const arma::mat & y, const arma::mat & mu) const {return 2 * (utils::xlogx(y) - y % arma::log(mu) - (y - mu));}
bool QuasiPoisson::validmu (const arma::mat & mu) const {return arma::all(arma::all(mu > 0));}
bool QuasiPoisson::valideta (const arma::mat & eta) const {return true;}
