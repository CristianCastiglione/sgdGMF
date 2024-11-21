// family.cpp
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 21/11/2024

#include "family.h"

using namespace glm;

// Gaussian family
arma::mat Gaussian::variance (const arma::mat & mu) const {return arma::ones(size(mu));}
arma::mat Gaussian::initialize (const arma::mat & y) const {return y;}
arma::mat Gaussian::devresid (const arma::mat & y, const arma::mat & mu) const {return arma::square(y - mu);}

// Binomial family
arma::mat Binomial::variance (const arma::mat & mu) const {return mu % (1 - mu);}
arma::mat Binomial::initialize (const arma::mat & y) const {return 2 * y - 1;}
arma::mat Binomial::devresid (const arma::mat & y, const arma::mat & mu) const {
    return - 2 * (y % arma::log(mu) + (1 - y) % arma::log1p(-mu));
}

// Poisson family
arma::mat Poisson::variance (const arma::mat & mu) const {return mu;}
arma::mat Poisson::initialize (const arma::mat & y) const {return this->linkfun(arma::clamp(y, 0.1, infty));}
arma::mat Poisson::devresid (const arma::mat & y, const arma::mat & mu) const {
    return 2 * (utils::xlogx(y) - y % arma::log(mu) - (y - mu));
}

// Gamma family
arma::mat Gamma::variance (const arma::mat & mu) const {return arma::square(mu);}
arma::mat Gamma::initialize (const arma::mat & y) const {return this->linkfun(y);}
arma::mat Gamma::devresid (const arma::mat & y, const arma::mat & mu) const {
    return - 2 * (arma::log(y / mu) - (y - mu) / mu);
}

// Inverse-Gaussian family
arma::mat InverseGaussian::variance (const arma::mat & mu) const {return mu % mu % mu;}
arma::mat InverseGaussian::initialize (const arma::mat & y) const {return this->linkfun(y);}
arma::mat InverseGaussian::devresid (const arma::mat & y, const arma::mat & mu) const {
    return arma::square(y - mu) / (y % mu % mu);
}

// Negative-Binomial family
arma::mat NegativeBinomial::variance (const arma::mat & mu) const {return mu + (mu % mu) / this->dispersion;}
arma::mat NegativeBinomial::initialize (const arma::mat & y) const {return this->linkfun(arma::clamp(y, 0.1, infty));}
arma::mat NegativeBinomial::devresid (const arma::mat & y, const arma::mat & mu) const {
    const double phi = this->dispersion;
    return 2 * (utils::xlogx(y) - y % arma::log(mu) - (y + phi) % (arma::log(y + phi) - arma::log(mu + phi)));
}

// Quasi-Binomial family
arma::mat QuasiBinomial::variance (const arma::mat & mu) const {return mu % (1 - mu);}
arma::mat QuasiBinomial::initialize (const arma::mat & y) const {return 2 * y - 1;}
arma::mat QuasiBinomial::devresid (const arma::mat & y, const arma::mat & mu) const {
    return - 2 * (y % arma::log(mu) + (1 - y) % arma::log1p(-mu));
}

// Quasi-Poisson family
arma::mat QuasiPoisson::variance (const arma::mat & mu) const {return mu;}
arma::mat QuasiPoisson::initialize (const arma::mat & y) const {return this->linkfun(arma::clamp(y, 0.1, infty));}
arma::mat QuasiPoisson::devresid (const arma::mat & y, const arma::mat & mu) const {
    return 2 * (utils::xlogx(y) - y % arma::log(mu) - (y - mu));
}

// Quasi family
arma::mat Quasi::variance (const arma::mat & mu) const {return this->varfun(mu);}
arma::mat Quasi::initialize (const arma::mat & y) const {return this->linkfun(this->initfun(y));}
arma::mat Quasi::devresid (const arma::mat & y, const arma::mat & mu) const {return this->devfun(y, mu);}
