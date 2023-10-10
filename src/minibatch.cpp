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


void ChunkPile::fill_tovisit () {
    this->tovisit = this->visited;
}

void ChunkPile::empty_visited () {
    this->visited = {};
}

void ChunkPile::pop_tovisit (const int & id) {
    int n = this->tovisit.n_elem;
    arma::uvec h = arma::find(this->tovisit == id);
    int i = h(0);
    if (i ==   0) {this->tovisit = this->tovisit.tail(n-1);}
    if (i == n-1) {this->tovisit = this->tovisit.head(n-1);}
    if (i > 0 && i < n-1) {
        arma::uvec head = this->tovisit.head(i);
        arma::uvec tail = this->tovisit.tail(n-i-1);
        this->tovisit = arma::join_cols(head, tail);
    }
}

void ChunkPile::push_visited (const int & id) {
    arma::uword i = id;
    this->visited = arma::join_cols(this->visited, arma::uvec{i});
}

void ChunkPile::sample_idx () {
    int n = this->tovisit.n_elem;
    int which;
    if (this->random) {
        which = arma::randi<int>(arma::distr_param(0, n-1));
    } else {
        which = 0;
    }
    this->idx = this->tovisit(which);
}

void ChunkPile::update () {
    // If tovisit is empty, fill it using visited and empty the later
    int n = this->tovisit.n_elem;
    if (n == 0) {
        this->fill_tovisit();
        this->empty_visited();
    }

    // Sample a random index, pop it from to visit and push it to visited
    this->sample_idx();
    this->pop_tovisit(idx);
    this->push_visited(idx);
}