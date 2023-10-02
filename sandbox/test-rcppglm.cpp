// rcppglm.cpp
// author: Cristian Castiglione
// creation: 01/10/2023
// last change: 01/10/2023

#include "rcppglm.h"

// [[Rcpp::export]]
void test_make_gaussian (std::string linkname) {
    Rcpp::XPtr<Family::Gaussian> family = make_gaussian(linkname);
    Rcpp::Rcout << "family: " << family->family << "\n";
    Rcpp::Rcout << "link: " << family->link << "\n";
    Rcpp::Rcout << "dispersion: " << family->dispersion << "\n";
}

// [[Rcpp::export]]
void test_make_binomial (std::string linkname) {
    Rcpp::XPtr<Family::Binomial> family = make_binomial(linkname);
    Rcpp::Rcout << "family: " << family->family << "\n";
    Rcpp::Rcout << "link: " << family->link << "\n";
    Rcpp::Rcout << "dispersion: " << family->dispersion << "\n";
}

// [[Rcpp::export]]
void test_make_poisson (std::string linkname) {
    Rcpp::XPtr<Family::Poisson> family = make_poisson(linkname);
    Rcpp::Rcout << "family: " << family->family << "\n";
    Rcpp::Rcout << "link: " << family->link << "\n";
    Rcpp::Rcout << "dispersion: " << family->dispersion << "\n";
}

// [[Rcpp::export]]
void test_make_gamma (std::string linkname) {
    Rcpp::XPtr<Family::Gamma> family = make_gamma(linkname);
    Rcpp::Rcout << "family: " << family->family << "\n";
    Rcpp::Rcout << "link: " << family->link << "\n";
    Rcpp::Rcout << "dispersion: " << family->dispersion << "\n";
}