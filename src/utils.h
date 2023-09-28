// utils.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 28/09/2023

#ifndef UTILS_H
#define UTILS_H

// #include <math.h>
#include <RcppArmadillo.h>

const double pi = M_PI;
const double invpi = 1.0 / M_PI;
const double log2pi = std::log(2.0 * M_PI);
const double sqrt2 = std::sqrt(2.0);
const double sqrtpi = std::sqrt(M_PI);
const double sqrt2pi = std::sqrt(2.0 * M_PI);

namespace utils {

// Maximum relative difference between two scalars/vectors
double absmax (const double & u, const double & v);
double absmax (const arma::vec & u, const arma::vec & v);

// Truncated representation of a vertor x, such that a <= x[i] <= b
void trim (arma::mat & x, double a, double b);
// arma::mat trim (const arma::mat & x, double a, double b);

// Stable calculation of log(1 + exp(x))
arma::mat log1pexp (const arma::mat & x);

// Stable calculation of log(1 - exp(-x))
arma::mat log1mexp (const arma::mat & x);

// Logistic transformation
arma::mat logit (const arma::mat & x);

// Inverse of the logistic transformation
arma::mat expit (const arma::mat & x);
arma::mat expit2 (const arma::mat & x);
arma::mat expitn (const arma::mat & x, double n);

// Complementary log-log and exp-exp transformation
arma::mat cloglog (const arma::mat & x);
arma::mat cexpexp (const arma::mat & x);

// Log-log and exp-exp transformations
arma::mat loglog (const arma::mat & x);
arma::mat expexp (const arma::mat & x);

// Standard Gaussian probability and cumulative density function
arma::mat pdfn (const arma::mat & x);
arma::mat cdfn (const arma::mat & x);

// Standard Gaussian log-probability and cumulative density function
arma::mat logpdfn (const arma::mat & x);
arma::mat logcdfn (const arma::mat & x);

// Gamma function
arma::mat gamma (const arma::mat & x);
arma::mat loggamma (const arma::mat & x);
arma::mat digamma (const arma::mat & x);
arma::mat trigamma (const arma::mat & x);

// Beta function
arma::mat beta (const arma::mat & x, const arma::mat & y);
arma::mat logbeta (const arma::mat & x, const arma::mat & y);
arma::mat dibeta (const arma::mat & x, const arma::mat & y);
arma::mat tribeta (const arma::mat & x, const arma::mat & y);

// Hinge loss function
arma::mat hinge (const arma::mat & x);

// Delta function
arma::mat dirac (const arma::mat & x, double a);

// Step function
arma::mat step (const arma::mat & x, double a, bool lower);

// Extract the half-vectorization of the square matrix M
arma::vec vech (const arma::mat & A);

}

#endif
