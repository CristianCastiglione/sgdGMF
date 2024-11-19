// misc.h
// author: Cristian Castiglione
// creation: 30/09/2023
// last change: 16/11/2024

#include <RcppArmadillo.h>
#include <time.h>
#include "utils.h"
#include "link.h"
#include "variance.h"
#include "family.h"
#include <memory>

using namespace glm;

// Create a dynamic pointer to an appropriate link/family class starting  
// from a string identifying the correct link/family to chose
std::unique_ptr<Link> make_link (const std::string & linkname);
std::unique_ptr<Variance> make_varf (const std::string & varname);
std::unique_ptr<Family> make_family (
    const std::string & familyname, 
    const std::string & linkname, 
    const std::string & varfname);

// Set the lower and upper bounds for mu and eta based on the observed data range
// so as to avoid to produce prediction with too extreme values 
// template<class F, class L>
void set_data_bounds (
    double & mulo, double & muup, double & etalo, double & etaup, 
    const double & eps, const double & ymin, const double & ymax, 
    const std::unique_ptr<Family> & family);

// Set the linear predictor trimming the extreme values
void set_eta (
    arma::mat & eta, const arma::mat & offset,
    const arma::mat & u, const arma::mat & v, 
    const double & etamin, const double & etamax);

// Get the linear predictor trimming the extreme values
arma::mat get_eta (
    const arma::mat & offset,
    const arma::mat & u, const arma::mat & v, 
    const double & etamin, const double & etamax);

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

// Print the optimization state
void print_state (
    const int & iter, const double & div, 
    const double & change, const double & time);

// Print the optimization state
void print_state (
    const int & iter, const double & div, 
    const double & change, const double & time,
    const double & scanned);

// Divide the data indices in random chunks
std::list<arma::uvec> sample_chunks (
    const int & n, const int & size, const bool & randomize);

// Select the appropriate chunk for the current iteration
int select_chunk (const int & iter, const int & nchunks);