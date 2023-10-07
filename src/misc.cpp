// misc.cpp
// author: Cristian Castiglione
// creation: 30/09/2023
// last change: 01/10/2023

#include "misc.h"


std::unique_ptr<Link::Link> make_link (const std::string & linkname) {
    std::unique_ptr<Link::Link> ptr;
    if (linkname == "identity") { ptr = std::make_unique<Link::Identity>();
    } else if (linkname == "logit") { ptr = std::make_unique<Link::Logit>();
    } else if (linkname == "probit") { ptr = std::make_unique<Link::Probit>();
    } else if (linkname == "cauchit") { ptr = std::make_unique<Link::Cauchy>();
    } else if (linkname == "cloglog") { ptr = std::make_unique<Link::cLogLog>();
    } else if (linkname == "log") { ptr = std::make_unique<Link::Log>();
    } else if (linkname == "inverse") { ptr = std::make_unique<Link::Inverse>();
    } else if (linkname == "sqrt") { ptr = std::make_unique<Link::Sqrt>();
    } else { Rcpp::stop("Link function not available."); }
    return ptr;
}

std::unique_ptr<Family::Family> make_family (const std::string & familyname, const std::string & linkname) {
    std::unique_ptr<Link::Link> link = make_link(linkname);
    std::unique_ptr<Family::Family> family;
    if (familyname == "gaussian") { family = std::make_unique<Family::Gaussian>(link);
    } else if (familyname == "binomial") { family = std::make_unique<Family::Binomial>(link);
    } else if (familyname == "poisson") { family = std::make_unique<Family::Poisson>(link);
    } else if (familyname == "gamma") { family = std::make_unique<Family::Gamma>(link);
    } else { Rcpp::stop("Family not available."); }
    return family;
}

void set_data_bounds (
    double & mulo, double & muup, double & etalo, double & etaup, 
    const double & eps, const double & ymin, const double & ymax, 
    const std::unique_ptr<Family::Family> & family
) {
    // We compute the lower and upper bounds on matrices of dim 1x1,
    // since we need to back-transform them using the family->linkfun()
    // method, which is defined only for arma::mat objects
    arma::mat mulot(1,1), muupt(1,1);
    arma::mat etalot(1,1), etaupt(1,1);
    
    // Mean bounds
    muupt(0,0) = ymax - eps * (ymax - ymin);
    mulot(0,0) = ymin + eps * (ymax - ymin);
    
    // Linear predictor bounds
    etalot = family->linkfun(mulot);
    etaupt = family->linkfun(muupt);

    // Inplace assignment
    mulo = mulot(0,0); muup = muupt(0,0);
    etalo = etalot(0,0); etaup = etaupt(0,0);
}

void set_uv_matrices (
    arma::mat & u, arma::mat & v,
    const arma::mat & A, const arma::mat & Z,
    const arma::mat & X, const arma::mat & B,
    const arma::mat & U, const arma::mat & V
) {
    u = arma::join_rows(X, A, U);
    v = arma::join_rows(B, Z, V);
}

void set_uv_indices (
    arma::uvec & idu, arma::uvec & idv, 
    const int & p, const int & q, const int & d
) {
    idu = join_cols(arma::linspace<arma::uvec>(p, p+q-1, q), arma::linspace<arma::uvec>(p+q, p+q+d-1, d));
    idv = join_cols(arma::linspace<arma::uvec>(0, p-1, p), arma::linspace<arma::uvec>(p+q, p+q+d-1, d));
}

void set_uv_penalty (
    arma::vec & penu, arma::vec & penv, const arma::vec & pen,
    const int & p, const int & q, const int & d
) {
    double penA, penB, penU, penV;
    penA = pen(0); penB = pen(1); penU = pen(2); penV = pen(3);
    penu = join_cols(arma::zeros(p), penA * arma::ones(q), penU * arma::ones(d));
    penv = join_cols(penB * arma::ones(p), arma::zeros(q), penV * arma::ones(d));
}

double exetime (const clock_t & start, const clock_t & end) {
    return static_cast<double>(end - start) / CLOCKS_PER_SEC;
}

void print_state (
    const int & iter, const double & dev, 
    const double & change, const double & time
) {
    if (time < 60) {
        std::printf(" %9i %11.2f %8.4f %8.2f s \n", iter, dev, change, time);
    } else {
        std::printf(" %9i %11.2f %8.4f %8.2f m \n", iter, dev, change, time/60);
    }
}



std::list<arma::uvec> sample_chunks (
    const int & n, const int & size, const bool & randomize
) {
    // If the data are ordered, it might be useful to shuffle the observations
    arma::uvec idx(n);
    idx = arma::linspace<arma::uvec>(0, n-1, n);
    if (randomize) {idx = arma::shuffle(idx);}
    
    // Instantiate the chunk-specific starting and ending indices
    int nchunks, start, end, range;
    nchunks = ceil(double(n) / size);
    
    // Partition the data into chunks of (almost) the same size
    // The last chunk could have a different size depending on 'n' and 'size'
    arma::uvec chunk;
    std::list<arma::uvec> chunks;
    for (int i = 0; i < nchunks; i++) {
        start = i * size;
        end = std::min((i + 1) * size, n);
        range = end - start;
        chunk = arma::linspace<arma::uvec>(start, end-1, range);
        chunks.push_back(idx(chunk));
    }

    // Return the chunk partition
    return chunks;
}

int select_chunk (const int & iter, const int & nchunks) {
    int mod = iter % nchunks;
    int idx = (mod == 0) ? nchunks : mod;
    return idx;
}