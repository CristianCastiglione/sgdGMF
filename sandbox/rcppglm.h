// rcppglm.cpp
// author: Cristian Castiglione
// creation: 30/09/2023
// last change: 30/09/2023

#ifndef RCPPGLM_H
#define RCPPGLM_H

#include <RcppArmadillo.h>
#include "link.h"
#include "family.h"

// Create a Link::Link object from a string identifying the link function
// std::unique_ptr<Link::Link> make_link(const std::string & link);
std::unique_ptr<Link::Link> make_link(const std::string & linkname);

// Create a parametrized Family::Family object with link function identified by a string
template<class F> inline Rcpp::XPtr<F> make_family(std::string linkname);

// Create a Family::Gaussian/Binomial/Poisson/Gamma object from a link function string
Rcpp::XPtr<Family::Gaussian> make_gaussian(std::string linkname);
Rcpp::XPtr<Family::Binomial> make_binomial(std::string linkname);
Rcpp::XPtr<Family::Poisson> make_poisson(std::string linkname);
Rcpp::XPtr<Family::Gamma> make_gamma(std::string linkname);

#endif