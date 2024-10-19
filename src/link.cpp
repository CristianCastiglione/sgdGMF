// link.cpp
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 30/09/2023

#include "link.h"

using namespace glm;

// Identity link
arma::mat Identity::linkfun (const arma::mat & mu) {return mu;}
arma::mat Identity::linkinv (const arma::mat & eta) {return eta;}
arma::mat Identity::mueta (const arma::mat & eta) {return arma::ones(arma::size(eta));}

// Logit link
arma::mat Logit::linkfun (const arma::mat & mu) {return arma::log(mu) - arma::log1p(-mu);}
arma::mat Logit::linkinv (const arma::mat & eta) {return arma::exp(eta - arma::log1p(arma::exp(eta)));}
arma::mat Logit::mueta (const arma::mat & eta) {return arma::exp(eta - 2 * arma::log1p(arma::exp(eta)));}

// Probit link
arma::mat Probit::linkfun (const arma::mat & mu) {
    // This code should be replaced with a vectorized implementation
    arma::mat eta = mu;
    eta.transform([](double & x) {return R::qnorm(x, 0, 1, true, false);});
    return eta;
}
arma::mat Probit::linkinv (const arma::mat & eta) {return arma::normcdf(eta);}
arma::mat Probit::mueta (const arma::mat & eta) {return arma::normpdf(eta);}

// Cauchy link
arma::mat Cauchy::linkfun (const arma::mat & mu) {return arma::tan(pi * (mu - 0.5));}
arma::mat Cauchy::linkinv (const arma::mat & eta) {return 0.5 + arma::atan(eta) / pi;}
arma::mat Cauchy::mueta (const arma::mat & eta) {return invpi / (eta % eta + 1);}

// cLogLog link
arma::mat cLogLog::linkfun (const arma::mat & mu) {return arma::log(- arma::log1p(-mu));}
arma::mat cLogLog::linkinv (const arma::mat & eta) {return 1 - arma::exp(- arma::exp(eta));}
arma::mat cLogLog::mueta (const arma::mat & eta) {return arma::exp(- eta - arma::exp(-eta));}

// Log link
arma::mat Log::linkfun (const arma::mat & mu) {return arma::log(mu);}
arma::mat Log::linkinv (const arma::mat & eta) {return arma::exp(eta);}
arma::mat Log::mueta (const arma::mat & eta) {return arma::exp(eta);}

// Inverse link
arma::mat Inverse::linkfun (const arma::mat & mu) {return 1 / mu;}
arma::mat Inverse::linkinv (const arma::mat & eta) {return 1 / eta;}
arma::mat Inverse::mueta (const arma::mat & eta) {return - 1 / (eta % eta);}

// Sqrt link
arma::mat Sqrt::linkfun (const arma::mat & mu) {return arma::sqrt(mu);}
arma::mat Sqrt::linkinv (const arma::mat & eta) {return arma::square(eta);}
arma::mat Sqrt::mueta (const arma::mat & eta) {return 2 * eta;}

