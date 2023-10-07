// minibatch.cpp
// author: Cristian Castiglione
// creation: 06/10/2023
// last change: 06/10/2023

#include "minibatch.h"

// [[Rcpp::export]]
arma::uvec c_get_chunk (
    const int & iter, const int & n, 
    const int & size, const bool & randomize
) {
    Chunks chunks;
    chunks.set_chunks(n, size, randomize);
    return chunks.get_chunk(iter);
}

// [[Rcpp::export]]
std::list<arma::uvec> c_get_chunks (
    const arma::uvec & iters, const int & n, 
    const int & size, const bool & randomize
) {
    Chunks chunks;
    chunks.set_chunks(n, size, randomize);
    return chunks.get_chunks(iters);
}
