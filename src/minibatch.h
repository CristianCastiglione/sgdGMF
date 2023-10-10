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
        bool randomize; // should we reshuffle the indices?
        arma::uvec idx;     // data index vector
        arma::uvec start;   // vector of starting indices (of idx) for each chunk
        arma::uvec end;     // vector of ending indices (of idx) for each chunk
        arma::uvec range;   // vector of lengths for each chunk

        // Get the data indices corresponding to the chunk at iteration 'iter'
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

class ChunkPile {
    public:
        int idx;
        bool random;
        arma::uvec tovisit;
        arma::uvec visited;

        void fill_tovisit ();
        void empty_visited ();
        void pop_tovisit (const int & id);
        void push_visited (const int & id);
        void sample_idx ();
        void update ();

        // Class constructor
        ChunkPile () {}
        ChunkPile (const int & n, const bool & rnd)  {
            this->idx = -1;
            this->random = rnd;
            this->tovisit = arma::linspace<arma::uvec>(0, n-1, n);
            this->visited = {};
        }
};


#endif
