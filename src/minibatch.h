// minibatch.h
// author: Cristian Castiglione
// creation: 06/10/2023
// last change: 06/10/2023

#ifndef MINIBATCH_H
#define MINIBATCH_H

#include <RcppArmadillo.h>

class Chunks {
    public:
        int nidx;   // number of observations
        int nchunks;    // number of chunks
        bool randomize; // are the indices randomized?
        arma::uvec idx;     // data index vector
        arma::uvec start;   // vector of starting indices (of idx) for each chunk
        arma::uvec end;     // vector of ending indices (of idx) for each chunk
        arma::uvec range;   // vector of lengths for each chunk

        // Get the data indices corresponding the chunk at iteration 'iter'
        arma::uvec get_chunk (const int & iter);

        // Get the list of data indices corresponding each chunk in the partition
        std::list<arma::uvec> get_chunks (const arma::uvec & iters);

        // Set all the chunks via index partition
        void set_chunks (const int & n, const int & size, const bool & randomize);

        // Class constructor
        Chunks () {}
        Chunks (const int & n, const int & size, const bool & randomize) {
            this->set_chunks(n, size, randomize);
        }
};

#endif
