// test-family.cpp
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 30/09/2023

#include "family.h"

// [[Rcpp::export]]
arma::vec c_gaussian_variance (const arma::vec & mu) {
    std::unique_ptr<Link::Link> ptr = std::make_unique<Link::Identity>();
    Family::Gaussian f(ptr);
    return f.variance(mu);
}

// [[Rcpp::export]]
arma::vec c_gaussian_initialize (const arma::vec & y) {
    std::unique_ptr<Link::Link> ptr = std::make_unique<Link::Identity>();
    Family::Gaussian f(ptr);
    return f.initialize(y);
}

// [[Rcpp::export]]
arma::vec c_gaussian_devresid (const arma::vec & y, const arma::vec & mu) {
    std::unique_ptr<Link::Link> ptr = std::make_unique<Link::Identity>();
    Family::Gaussian f(ptr);
    return f.devresid(y, mu);
}

// [[Rcpp::export]]
arma::vec c_binomial_variance (const arma::vec & mu) {
    std::unique_ptr<Link::Link> ptr = std::make_unique<Link::Logit>();
    Family::Binomial f(ptr);
    return f.variance(mu);
}

// [[Rcpp::export]]
arma::vec c_binomial_initialize (const arma::vec & y) {
    std::unique_ptr<Link::Link> ptr = std::make_unique<Link::Logit>();
    Family::Binomial f(ptr);
    return f.initialize(y);
}

// [[Rcpp::export]]
arma::vec c_binomial_devresid (const arma::vec & y, const arma::vec & mu) {
    std::unique_ptr<Link::Link> ptr = std::make_unique<Link::Logit>();
    Family::Binomial f(ptr);
    return f.devresid(y, mu);
}

// [[Rcpp::export]]
arma::vec c_poisson_variance (const arma::vec & mu) {
    std::unique_ptr<Link::Link> ptr = std::make_unique<Link::Log>();
    Family::Poisson f(ptr);
    return f.variance(mu);
}

// [[Rcpp::export]]
arma::vec c_poisson_initialize (const arma::vec & y) {
    std::unique_ptr<Link::Link> ptr = std::make_unique<Link::Log>();
    Family::Poisson f(ptr);
    return f.initialize(y);
}

// [[Rcpp::export]]
arma::vec c_poisson_devresid (const arma::vec & y, const arma::vec & mu) {
    std::unique_ptr<Link::Link> ptr = std::make_unique<Link::Log>();
    Family::Poisson f(ptr);
    return f.devresid(y, mu);
}

// [[Rcpp::export]]
arma::vec c_gamma_variance (const arma::vec & mu) {
    std::unique_ptr<Link::Link> ptr = std::make_unique<Link::Log>();
    Family::Gamma f(ptr);
    return f.variance(mu);
}

// [[Rcpp::export]]
arma::vec c_gamma_initialize (const arma::vec & y) {
    std::unique_ptr<Link::Link> ptr = std::make_unique<Link::Log>();
    Family::Gamma f(ptr);
    return f.initialize(y);
}

// [[Rcpp::export]]
arma::vec c_gamma_devresid (const arma::vec & y, const arma::vec & mu) {
    std::unique_ptr<Link::Link> ptr = std::make_unique<Link::Log>();
    Family::Gamma f(ptr);
    return f.devresid(y, mu);
}