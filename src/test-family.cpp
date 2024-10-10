// test-family.cpp
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 30/09/2023

#include "family.h"

using namespace glm;

//' @keywords internal
// [[Rcpp::export("cpp.family.gaussian.variance")]]
arma::vec cpp_gaussian_variance (const arma::vec & mu) {
    std::unique_ptr<Link> ptr = std::make_unique<Identity>();
    Gaussian f(ptr);
    return f.variance(mu);
}

//' @keywords internal
// [[Rcpp::export("cpp.family.gaussian.initialize")]]
arma::vec cpp_gaussian_initialize (const arma::vec & y) {
    std::unique_ptr<Link> ptr = std::make_unique<Identity>();
    Gaussian f(ptr);
    return f.initialize(y);
}

//' @keywords internal
// [[Rcpp::export("cpp.family.gaussian.devresid")]]
arma::vec cpp_gaussian_devresid (const arma::vec & y, const arma::vec & mu) {
    std::unique_ptr<Link> ptr = std::make_unique<Identity>();
    Gaussian f(ptr);
    return f.devresid(y, mu);
}

//' @keywords internal
// [[Rcpp::export("cpp.family.binomial.variance")]]
arma::vec cpp_binomial_variance (const arma::vec & mu) {
    std::unique_ptr<Link> ptr = std::make_unique<Logit>();
    Binomial f(ptr);
    return f.variance(mu);
}

//' @keywords internal
// [[Rcpp::export("cpp.family.binomial.initialize")]]
arma::vec cpp_binomial_initialize (const arma::vec & y) {
    std::unique_ptr<Link> ptr = std::make_unique<Logit>();
    Binomial f(ptr);
    return f.initialize(y);
}

//' @keywords internal
// [[Rcpp::export("cpp.family.binomial.devresid")]]
arma::vec cpp_binomial_devresid (const arma::vec & y, const arma::vec & mu) {
    std::unique_ptr<Link> ptr = std::make_unique<Logit>();
    Binomial f(ptr);
    return f.devresid(y, mu);
}

//' @keywords internal
// [[Rcpp::export("cpp.family.poisson.variance")]]
arma::vec cpp_poisson_variance (const arma::vec & mu) {
    std::unique_ptr<Link> ptr = std::make_unique<Log>();
    Poisson f(ptr);
    return f.variance(mu);
}

//' @keywords internal
// [[Rcpp::export("cpp.family.poisson.initialize")]]
arma::vec cpp_poisson_initialize (const arma::vec & y) {
    std::unique_ptr<Link> ptr = std::make_unique<Log>();
    Poisson f(ptr);
    return f.initialize(y);
}

//' @keywords internal
// [[Rcpp::export("cpp.family.poisson.devresid")]]
arma::vec cpp_poisson_devresid (const arma::vec & y, const arma::vec & mu) {
    std::unique_ptr<Link> ptr = std::make_unique<Log>();
    Poisson f(ptr);
    return f.devresid(y, mu);
}

//' @keywords internal
// [[Rcpp::export("cpp.family.gamma.variance")]]
arma::vec cpp_gamma_variance (const arma::vec & mu) {
    std::unique_ptr<Link> ptr = std::make_unique<Log>();
    Gamma f(ptr);
    return f.variance(mu);
}

//' @keywords internal
// [[Rcpp::export("cpp.family.gamma.initialize")]]
arma::vec cpp_gamma_initialize (const arma::vec & y) {
    std::unique_ptr<Link> ptr = std::make_unique<Log>();
    Gamma f(ptr);
    return f.initialize(y);
}

//' @keywords internal
// [[Rcpp::export("cpp.family.gamma.devresid")]]
arma::vec cpp_gamma_devresid (const arma::vec & y, const arma::vec & mu) {
    std::unique_ptr<Link> ptr = std::make_unique<Log>();
    Gamma f(ptr);
    return f.devresid(y, mu);
}

//' @keywords internal
// [[Rcpp::export("cpp.family.negbinom.variance")]]
arma::vec cpp_negbinom_variance (const arma::vec & mu) {
    std::unique_ptr<Link> ptr = std::make_unique<Log>();
    NegativeBinomial f(ptr);
    return f.variance(mu);
}

//' @keywords internal
// [[Rcpp::export("cpp.family.negbinom.initialize")]]
arma::vec cpp_negbinom_initialize (const arma::vec & y) {
    std::unique_ptr<Link> ptr = std::make_unique<Log>();
    NegativeBinomial f(ptr);
    return f.initialize(y);
}

//' @keywords internal
// [[Rcpp::export("cpp.family.negbinom.devresid")]]
arma::vec cpp_negbinom_devresid (const arma::vec & y, const arma::vec & mu) {
    std::unique_ptr<Link> ptr = std::make_unique<Log>();
    NegativeBinomial f(ptr);
    return f.devresid(y, mu);
}