// utils.h
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 19/11/2024

#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>

// Standard constants used in numerical computations
const double pi = M_PI;
const double invpi = 1.0 / M_PI;
const double log2pi = std::log(2.0 * M_PI);
const double sqrt2 = std::sqrt(2.0);
const double sqrtpi = std::sqrt(M_PI);
const double sqrt2pi = std::sqrt(2.0 * M_PI);
const double infty = arma::datum::inf;

// The following vectors of coefficients serves to approximate the 
// quantile function of a standard normal distribution.
// For more details see: https://ar5iv.labs.arxiv.org/html/1002.0567
const arma::vec qn_inner_coef = {
    + 0.195740115269792, - 0.652871358365296, + 1.246899760652504,
    + 0.155331081623168, - 0.839293158122257};
const arma::vec qn_tails_coef = {
    +16.682320830719986527, + 4.120411523939115059, + 0.029814187308200211, 
    - 1.000182518730158122, + 7.173787663925508066, + 8.759693508958633869};

namespace utils {

// Maximum relative difference between two scalars/vectors
double absmax (const double & u, const double & v);
double absmax (const arma::vec & u, const arma::vec & v);

// Truncated representation of a vector/matrix x, such that a <= x[i,j] <= b
void trim (arma::mat & x, const double & a, const double & b);
void trim (arma::mat & x, const double & a, const double & b, const arma::uvec & idx);
void trim (arma::mat & x, const double & a, const double & b, const arma::uvec & idx, const arma::uvec & idy);

// All and any operator for boolean matrices
bool all(const arma::umat & x);
bool any(const arma::umat & x);

// Lp norm of a vector/matrix
double norm (const arma::mat & x);
double norm (const arma::mat & x, const double & p);

// Pointwise maximum between 0 and x (to the power of p)
arma::mat max0 (const arma::mat & x);
arma::mat max0 (const arma::mat & x, const double & p);

// Stable calculation of x*log(x), with 0*log(0) = 0;
arma::mat xlogx (const arma::mat & x);

// Stable calculation of log(1 + exp(x))
arma::mat log1pexp (const arma::mat & x);

// Stable calculation of log(1 - exp(-x))
arma::mat log1mexp (const arma::mat & x);

// Logistic transformation
arma::mat logit (const arma::mat & x);

// Inverse of the logistic transformation
arma::mat expit (const arma::mat & x);
arma::mat expit2 (const arma::mat & x);
arma::mat expitn (const arma::mat & x, const double & n);

// Complementary log-log and exp-exp transformation
arma::mat cloglog (const arma::mat & x);
arma::mat cexpexp (const arma::mat & x);

// Log-log and exp-exp transformations
arma::mat loglog (const arma::mat & x);
arma::mat expexp (const arma::mat & x);

// Standard Gaussian probability and cumulative density function
arma::mat pdfn (const arma::mat & x);
arma::mat cdfn (const arma::mat & x);
arma::mat qdfn (const arma::mat & p);

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
arma::mat dirac (const arma::mat & x, const double & a);

// Step function
arma::mat step (const arma::mat & x, const double & a, const bool & lower);

// Extract the half-vectorization of the square matrix M
arma::vec vech (const arma::mat & A);

}

#endif
