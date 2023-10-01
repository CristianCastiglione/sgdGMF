// misc.h
// author: Cristian Castiglione
// creation: 30/09/2023
// last change: 30/09/2023

#include <RcppArmadillo.h>
#include <time.h>
#include "link.h"
#include "family.h"

// Set the lower and upper bounds for mu and eta based on the observed data range
// so as to avoid to produce prediction with too extreme values 
template<class F>
void set_data_bounds (
    double & mulo, double & muup, double & etalo, double & etaup, 
    const double & eps, const double & ymin, const double & ymax, F family);

// Set the augmented u and v matrices merging by column the fixed and latent effect matrices
void set_uv_matrices (
    arma::mat & u, arma::mat & v,
    const arma::mat & A, const arma::mat & Z,
    const arma::mat & X, const arma::mat & B,
    const arma::mat & U, const arma::mat & V);

// Set the indices of the parameters to be update along the optimization in u and v
void set_uv_indices (
    arma::uvec & idu, arma::uvec & idv, 
    const int & p, const int & q, const int & d);

// Set the vectors of penalty parameters to be multiplied columnwise to u and v
void set_uv_penalty (
    arma::vec & penu, arma::vec & penv, const arma::vec & pen,
    const int & p, const int & q, const int & d);

// Convert the cpu clock-time in the elapsed execution time (seconds)
double exetime (const clock_t & start, const clock_t & end);
