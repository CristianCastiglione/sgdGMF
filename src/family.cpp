// family.cpp
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 28/09/2023

#include "family.h"

// Gaussian family
template<class L>
arma::mat Family::Gaussian<L>::variance (const arma::mat & mu) const {return arma::ones(arma::size(mu));}

template<class L>
arma::mat Family::Gaussian<L>::initialize (const arma::mat & y) const {return y;}

template<class L>
bool Family::Gaussian<L>::validmu (const arma::mat & mu) const {return true;}

template<class L>
bool Family::Gaussian<L>::valideta (const arma::mat & eta) const {return true;}

template<class L>
arma::mat Family::Gaussian<L>::devresid (const arma::mat & y, const arma::mat & mu) const {
    return arma::square(y - mu);
}

template<class L>
arma::mat Family::Gaussian<L>::devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const {
    return wt % arma::square(y - mu);
}

// Binomial family
template<class L>
arma::mat Family::Binomial<L>::variance (const arma::mat & mu) const {return mu % (1 - mu);}

template<class L>
arma::mat Family::Binomial<L>::initialize (const arma::mat & y) const {return (2 * y - 1) + 0.25 * arma::randu(size(y));}

template<class L>
bool Family::Binomial<L>::validmu (const arma::mat & mu) const {return arma::all(arma::all(mu > 0)) && arma::all(arma::all(mu < 1));}

template<class L>
bool Family::Binomial<L>::valideta (const arma::mat & eta) const {return true;}

template<class L>
arma::mat Family::Binomial<L>::devresid (const arma::mat & y, const arma::mat & mu) const {
    return - 2 * (y % arma::log(mu) + (1 - y) % arma::log(1 - mu));
}

template<class L>
arma::mat Family::Binomial<L>::devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const {
    return - 2 * wt % (y % arma::log(mu) + (1 - y) % arma::log(1 - mu));
}

// Poisson family
template<class L>
arma::mat Family::Poisson<L>::variance (const arma::mat & mu) const {return mu;}

template<class L>
arma::mat Family::Poisson<L>::initialize (const arma::mat & y) const {return this->linkfun(y + 0.1);}

template<class L>
bool Family::Poisson<L>::validmu (const arma::mat & mu) const {return arma::all(arma::all(mu > 0));}

template<class L>
bool Family::Poisson<L>::valideta (const arma::mat & eta) const {return true;}

template<class L>
arma::mat Family::Poisson<L>::devresid (const arma::mat & y, const arma::mat & mu) const {
    return 2 * (y % arma::log(y / mu) - (y - mu));
}

template<class L>
arma::mat Family::Poisson<L>::devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const {
    return 2 * wt % (y % arma::log(y / mu) - (y - mu));
}

// Gamma family
template<class L>
arma::mat Family::Gamma<L>::variance (const arma::mat & mu) const {return arma::square(mu);}

template<class L>
arma::mat Family::Gamma<L>::initialize (const arma::mat & y) const {return this->linkfun(y);}

template<class L>
bool Family::Gamma<L>::validmu (const arma::mat & mu) const {return arma::all(arma::all(mu > 0));}

template<class L>
bool Family::Gamma<L>::valideta (const arma::mat & eta) const {return true;}

template<class L>
arma::mat Family::Gamma<L>::devresid (const arma::mat & y, const arma::mat & mu) const {
    return - 2 * (arma::log(y / mu) - (y - mu) / mu);
}

template<class L>
arma::mat Family::Gamma<L>::devresid (const arma::mat & y, const arma::mat & mu, const arma::mat & wt) const {
    return - 2 * wt % (arma::log(y / mu) - (y - mu) / mu);
}






// :::: TEST FUNCTIONS ::::


// [[Rcpp::export]]
arma::vec c_gaussian_variance (const arma::vec & mu) {
    Link::Identity l; Family::Gaussian f(l); return f.variance(mu);}
// [[Rcpp::export]]
arma::vec c_gaussian_initialize (const arma::vec & y) {
    Link::Identity l; Family::Gaussian f(l); return f.initialize(y);}
// [[Rcpp::export]]
arma::vec c_gaussian_devresid (const arma::vec & y, const arma::vec & mu) {
    Link::Identity l; Family::Gaussian f(l); return f.devresid(y, mu);
}

// [[Rcpp::export]]
arma::vec c_binomial_variance (const arma::vec & mu) {
    Link::Logit l; Family::Binomial f(l); return f.variance(mu);}
// [[Rcpp::export]]
arma::vec c_binomial_initialize (const arma::vec & y) {
    Link::Logit l; Family::Binomial f(l); return f.initialize(y);}
// [[Rcpp::export]]
arma::vec c_binomial_devresid (const arma::vec & y, const arma::vec & mu) {
    Link::Logit l; Family::Binomial f(l); return f.devresid(y, mu);
}

// [[Rcpp::export]]
arma::vec c_poisson_variance (const arma::vec & mu) {
    Link::Log l; Family::Poisson f(l); return f.variance(mu);}
// [[Rcpp::export]]
arma::vec c_poisson_initialize (const arma::vec & y) {
    Link::Log l; Family::Poisson f(l); return f.initialize(y);}
// [[Rcpp::export]]
arma::vec c_poisson_devresid (const arma::vec & y, const arma::vec & mu) {
    Link::Log l; Family::Poisson f(l); return f.devresid(y, mu);
}

// [[Rcpp::export]]
arma::vec c_gamma_variance (const arma::vec & mu) {
    Link::Log l; Family::Gamma f(l); return f.variance(mu);
}
// [[Rcpp::export]]
arma::vec c_gamma_initialize (const arma::vec & y) {
    Link::Log l; Family::Gamma f(l); return f.initialize(y);
}
// [[Rcpp::export]]
arma::vec c_gamma_devresid (const arma::vec & y, const arma::vec & mu) {
    Link::Log l; Family::Gamma f(l); return f.devresid(y, mu);
}
