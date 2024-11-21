// link.cpp
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 21/11/2024

#include "link.h"

using namespace glm;

// Identity link
bool Identity::valideta (const arma::mat & eta){return true;}
arma::mat Identity::linkfun (const arma::mat & mu) {return mu;}
arma::mat Identity::linkinv (const arma::mat & eta) {return eta;}
arma::mat Identity::mueta (const arma::mat & eta) {return arma::ones(arma::size(eta));}

// Logit link
bool Logit::valideta (const arma::mat & eta){return true;}
arma::mat Logit::linkfun (const arma::mat & mu) {return arma::log(mu) - arma::log1p(-mu);}
arma::mat Logit::linkinv (const arma::mat & eta) {return arma::exp(eta - arma::log1p(arma::exp(eta)));}
arma::mat Logit::mueta (const arma::mat & eta) {return arma::exp(eta - 2 * arma::log1p(arma::exp(eta)));}

// Probit link
bool Probit::valideta (const arma::mat & eta){return true;}
arma::mat Probit::linkfun (const arma::mat & mu) {
    // This code should be replaced with a vectorized implementation
    arma::mat eta = mu;
    eta.transform([](double & x) {return R::qnorm(x, 0, 1, true, false);});
    return eta;
}
arma::mat Probit::linkinv (const arma::mat & eta) {return arma::normcdf(eta);}
arma::mat Probit::mueta (const arma::mat & eta) {return arma::normpdf(eta);}

// Cauchit link
bool Cauchit::valideta (const arma::mat & eta){return true;}
arma::mat Cauchit::linkfun (const arma::mat & mu) {return arma::tan(pi * (mu - 0.5));}
arma::mat Cauchit::linkinv (const arma::mat & eta) {return 0.5 + arma::atan(eta) / pi;}
arma::mat Cauchit::mueta (const arma::mat & eta) {return invpi / (eta % eta + 1);}

// cLogLog link
bool cLogLog::valideta (const arma::mat & eta){return true;}
arma::mat cLogLog::linkfun (const arma::mat & mu) {return arma::log(- arma::log1p(-mu));}
arma::mat cLogLog::linkinv (const arma::mat & eta) {return 1 - arma::exp(- arma::exp(eta));}
arma::mat cLogLog::mueta (const arma::mat & eta) {return arma::exp(- eta - arma::exp(-eta));}

// Log link
bool Log::valideta (const arma::mat & eta){return true;}
arma::mat Log::linkfun (const arma::mat & mu) {return arma::log(mu);}
arma::mat Log::linkinv (const arma::mat & eta) {return arma::exp(eta);}
arma::mat Log::mueta (const arma::mat & eta) {return arma::exp(eta);}

// Inverse link
bool Inverse::valideta (const arma::mat & eta){return utils::all(eta > 0);}
arma::mat Inverse::linkfun (const arma::mat & mu) {return 1 / mu;}
arma::mat Inverse::linkinv (const arma::mat & eta) {return 1 / eta;}
arma::mat Inverse::mueta (const arma::mat & eta) {return - 1 / (eta % eta);}

// Squared inverse link
bool SquaredInverse::valideta (const arma::mat & eta){return utils::all(eta > 0);}
arma::mat SquaredInverse::linkfun (const arma::mat & mu) {return 1 / arma::square(mu);}
arma::mat SquaredInverse::linkinv (const arma::mat & eta) {return 1 / arma::sqrt(eta);}
arma::mat SquaredInverse::mueta (const arma::mat & eta) {return - 1 / (2 * arma::pow(eta, 1.5));}

// Sqrt link
bool Sqrt::valideta (const arma::mat & eta){return utils::all(eta > 0);}
arma::mat Sqrt::linkfun (const arma::mat & mu) {return arma::sqrt(mu);}
arma::mat Sqrt::linkinv (const arma::mat & eta) {return arma::square(eta);}
arma::mat Sqrt::mueta (const arma::mat & eta) {return 2 * eta;}

