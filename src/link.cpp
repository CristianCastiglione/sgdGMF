// link.cpp
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 28/09/2023

#include "link.h"

// Identity link
arma::mat Link::Identity::linkfun (const arma::mat & mu) const {return mu;}
arma::mat Link::Identity::linkinv (const arma::mat & eta) const {return eta;}
arma::mat Link::Identity::mueta (const arma::mat & eta) const {return arma::ones(arma::size(eta));}

// Logit link
arma::mat Link::Logit::linkfun (const arma::mat & mu) const {return arma::log(mu) - arma::log(1 - mu);}
arma::mat Link::Logit::linkinv (const arma::mat & eta) const {return arma::exp(eta - arma::log1p(arma::exp(eta)));}
arma::mat Link::Logit::mueta (const arma::mat & eta) const {return arma::exp(eta - 2 * arma::log1p(arma::exp(eta)));}

// Probit link
arma::mat Link::Probit::linkfun (const arma::mat & mu) const {
    // This code should be replaced with a vectorized implementation
    arma::mat eta = mu;
    eta.transform([](double & x) {return R::qnorm(x, 0, 1, true, false);});
    return eta;
}
arma::mat Link::Probit::linkinv (const arma::mat & eta) const {return arma::normcdf(eta);}
arma::mat Link::Probit::mueta (const arma::mat & eta) const {return arma::normpdf(eta);}

// Cauchy link
arma::mat Link::Cauchy::linkfun (const arma::mat & mu) const {return arma::tan(pi * (mu - 0.5));}
arma::mat Link::Cauchy::linkinv (const arma::mat & eta) const {return 0.5 + arma::atan(eta) / pi;}
arma::mat Link::Cauchy::mueta (const arma::mat & eta) const {return invpi / (eta % eta + 1);}

// cLogLog link
arma::mat Link::cLogLog::linkfun (const arma::mat & mu) const {return arma::log(- arma::log(1 - mu));}
arma::mat Link::cLogLog::linkinv (const arma::mat & eta) const {return 1 - arma::exp(- arma::exp(eta));}
arma::mat Link::cLogLog::mueta (const arma::mat & eta) const {return arma::exp(- eta - arma::exp(-eta));}

// Log link
arma::mat Link::Log::linkfun (const arma::mat & mu) const {return arma::log(mu);}
arma::mat Link::Log::linkinv (const arma::mat & eta) const {return arma::exp(eta);}
arma::mat Link::Log::mueta (const arma::mat & eta) const {return arma::exp(eta);}

// Inverse link
arma::mat Link::Inverse::linkfun (const arma::mat & mu) const {return 1 / mu;}
arma::mat Link::Inverse::linkinv (const arma::mat & eta) const {return 1 / eta;}
arma::mat Link::Inverse::mueta (const arma::mat & eta) const {return - 1 / (eta % eta);}

// Sqrt link
arma::mat Link::Sqrt::linkfun (const arma::mat & mu) const {return arma::sqrt(mu);}
arma::mat Link::Sqrt::linkinv (const arma::mat & eta) const {return arma::square(eta);}
arma::mat Link::Sqrt::mueta (const arma::mat & eta) const {return 2 * eta;}

