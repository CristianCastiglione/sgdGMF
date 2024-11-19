// variance.cpp
// author: Cristian Castiglione
// creation: 08/11/2023
// last change: 19/11/2024

#include "variance.h"

using namespace glm;

// Constant variance
bool Constant::validmu (const arma::mat & mu) {return true;}
arma::mat Constant::initfun (const arma::mat & y) {return y;}
arma::mat Constant::varfun (const arma::mat & mu, const double & phi) {return arma::ones(arma::size(mu));}
arma::mat Constant::devfun (const arma::mat & y, const arma::mat & mu, const double & phi) {
    return arma::square(y - mu);
}

// Linear variance
bool Linear::validmu (const arma::mat & mu) {return utils::all(mu > 0);}
arma::mat Linear::initfun (const arma::mat & y) {return arma::clamp(y, 0.1, infty);}
arma::mat Linear::varfun (const arma::mat & mu, const double & phi) {return mu;}
arma::mat Linear::devfun (const arma::mat & y, const arma::mat & mu, const double & phi) {
    return 2 * (utils::xlogx(y) - y % arma::log(mu) - (y - mu));
}

// Squared variance
bool Squared::validmu (const arma::mat & mu) {return utils::all(mu > 0);}
arma::mat Squared::initfun (const arma::mat & y) {return arma::clamp(y, 0.1, infty);}
arma::mat Squared::varfun (const arma::mat & mu, const double & phi) {return mu % mu;}
arma::mat Squared::devfun (const arma::mat & y, const arma::mat & mu, const double & phi) {
    return - 2 * (arma::log(y / mu) - (y - mu) / mu);
}

// Cubic variance
bool Cubic::validmu (const arma::mat & mu) {return utils::all(mu > 0);}
arma::mat Cubic::initfun (const arma::mat & y) {return arma::clamp(y, 0.1, infty);}
arma::mat Cubic::varfun (const arma::mat & mu, const double & phi) {return mu % mu % mu;}
arma::mat Cubic::devfun (const arma::mat & y, const arma::mat & mu, const double & phi) {
    return arma::square(y - mu) / (y % mu % mu);
}

// cSquared variance
bool cSquared::validmu (const arma::mat & mu) {return utils::all(mu > 0) && utils::all(mu < 1);}
arma::mat cSquared::initfun (const arma::mat & y) {return 0.90 * (y - 0.5) + 0.5;}
arma::mat cSquared::varfun (const arma::mat & mu, const double & phi) {return mu % (1 - mu);}
arma::mat cSquared::devfun (const arma::mat & y, const arma::mat & mu, const double & phi) {
    return - 2 * (y % arma::log(mu) + (1 - y) % arma::log1p(-mu));
}

// Negative-Binomial variance
bool NBVariance::validmu (const arma::mat & mu) {return utils::all(mu > 0);}
arma::mat NBVariance::initfun (const arma::mat & y) {return arma::clamp(y, 0.1, infty);}
arma::mat NBVariance::varfun (const arma::mat & mu, const double & phi) {return mu % (1 + mu / phi);}
arma::mat NBVariance::devfun (const arma::mat & y, const arma::mat & mu, const double & phi) {
    return 2 * (utils::xlogx(y) - y % arma::log(mu) - (y + phi) % (arma::log(y + phi) - arma::log(mu + phi)));
}
