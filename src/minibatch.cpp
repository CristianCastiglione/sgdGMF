// minibatch.cpp
// author: Cristian Castiglione
// creation: 06/10/2023
// last change: 06/10/2023

#include "minibatch.h"

void Chunks::set_chunks (const int & n, const int & size, const bool & randomize) {
    this->nidx = n;
    this->nchunks = ceil(double(n) / size);
    this->randomize = randomize;
    this->idx = arma::linspace<arma::uvec>(0, n-1, n);
    this->start = arma::zeros<arma::uvec>(this->nchunks);
    this->end = arma::zeros<arma::uvec>(this->nchunks);
    this->range = arma::zeros<arma::uvec>(this->nchunks);
    if (this->randomize) {
        this->idx = arma::shuffle(this->idx);
    }
    for (int i = 0; i < this->nchunks; i++) {
        this->start(i) = i * size;
        this->end(i) = std::min((i + 1) * size, n);
        this->range(i) = this->end(i) - this->start(i);
    }
}

arma::uvec Chunks::get_chunk (const int & iter) {
    int mod = iter % this->nchunks;
    // int i = (mod == 0) ? this->nchunks : mod;
    int i;
    if (iter == 0 && mod == 0) {i = iter;}
    if (iter != 0 && mod != 0) {i = mod;}
    if (iter != 0 && mod == 0) {i = 0;}
    int a = this->start(i);
    int b = this->end(i)-1;
    int c = this->range(i);
    arma::uvec which = arma::linspace<arma::uvec>(a, b, c);
    arma::uvec chunk = this->idx(which);
    return chunk;
}

std::list<arma::uvec> Chunks::get_chunks (const arma::uvec & iters) {
    arma::uvec chunk;
    std::list<arma::uvec> chunks;
    for (int iter : iters) {
        chunk = this->get_chunk(iter);
        chunks.push_back(chunk);
    }
    return chunks;
}
