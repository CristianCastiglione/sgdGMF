// misc.cpp
// author: Cristian Castiglione
// creation: 30/09/2023
// last change: 21/11/2024

#include "misc.h"

using namespace glm;

std::unique_ptr<Link> make_link (
    const std::string & linkname
) {
    bool flag = true;
    std::unique_ptr<Link> link;
    if (linkname == "identity") { flag = false; link = std::make_unique<Identity>(); }
    if (linkname == "logit") { flag = false; link = std::make_unique<Logit>(); }
    if (linkname == "probit") { flag = false; link = std::make_unique<Probit>(); }
    if (linkname == "cauchit") { flag = false; link = std::make_unique<Cauchit>(); }
    if (linkname == "cloglog") { flag = false; link = std::make_unique<cLogLog>(); }
    if (linkname == "log") { flag = false; link = std::make_unique<Log>(); }
    if (linkname == "inverse") { flag = false; link = std::make_unique<Inverse>(); }
    if (linkname == "1/mu^2") { flag = false; link = std::make_unique<SquaredInverse>(); }
    if (linkname == "sqrt") { flag = false; link = std::make_unique<Sqrt>(); }
    if (flag) { throw std::string("Link function not available"); }
    return link;
}

std::unique_ptr<Variance> make_varf (
    const std::string & varfname
) {
    bool flag = true;
    std::unique_ptr<Variance> varf;
    if (varfname == "const") { flag = false; varf = std::make_unique<Constant>(); }
    if (varfname == "mu") { flag = false; varf = std::make_unique<Linear>(); }
    if (varfname == "mu^2") { flag = false; varf = std::make_unique<Squared>(); }
    if (varfname == "mu^3") { flag = false; varf = std::make_unique<Cubic>(); }
    if (varfname == "mu(1-mu)") { flag = false; varf = std::make_unique<cSquared>(); }
    if (varfname == "mu(1+t*mu)") { flag = false; varf = std::make_unique<NBVariance>(); }
    if (flag) { throw std::string("Variance function not available"); }
    return varf;
}

std::unique_ptr<Family> make_family (
    const std::string & familyname, 
    const std::string & linkname,
    const std::string & varfname
) {
    bool flag = true;
    std::unique_ptr<Link> link = make_link(linkname);
    std::unique_ptr<Variance> varf;
    std::unique_ptr<Family> family;
    if (familyname == "gaussian") { 
        flag = false; 
        varf = make_varf("const");
        family = std::make_unique<Gaussian>(link, varf);
    }
    if (familyname == "binomial") { 
        flag = false; 
        varf = make_varf("mu(1-mu)");
        family = std::make_unique<Binomial>(link, varf);
    }
    if (familyname == "poisson") { 
        flag = false; 
        varf = make_varf("mu");
        family = std::make_unique<Poisson>(link, varf);
    }
    if (familyname == "gamma") { 
        flag = false; 
        varf = make_varf("mu^2");
        family = std::make_unique<Gamma>(link, varf);
    }
    if (familyname == "invgaussian") { 
        flag = false; 
        varf = make_varf("mu^3");
        family = std::make_unique<Gamma>(link, varf);
    }
    if (familyname == "negbinom") { 
        flag = false; 
        varf = make_varf("mu(1+t*mu)");
        family = std::make_unique<NegativeBinomial>(link, varf);
    }
    if (familyname == "quasibinomial") { 
        flag = false; 
        varf = make_varf("mu(1-mu)");
        family = std::make_unique<QuasiBinomial>(link, varf);
    }
    if (familyname == "quasipoisson") { 
        flag = false; 
        varf = make_varf("mu");
        family = std::make_unique<QuasiPoisson>(link, varf);
    }
    if (familyname == "quasi") { 
        flag = false; 
        varf = make_varf(varfname);
        family = std::make_unique<Quasi>(link, varf);
    }
    if (flag) { 
        throw std::string("Family not available"); 
    }
    return family;
}

void set_data_bounds (
    double & mulo, double & muup, double & etalo, double & etaup, 
    const double & eps, const double & ymin, const double & ymax, 
    const std::unique_ptr<Family> & family
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
    mulo = mulot(0,0); 
    muup = muupt(0,0);
    etalo = etalot(0,0); 
    etaup = etaupt(0,0);
}

void set_eta (
    arma::mat & eta, const arma::mat & offset, 
    const arma::mat & u, const arma::mat & v, 
    const double & etamin, const double & etamax
) {
    eta = offset + u * v.t();
    eta.clamp(etamin, etamax);
    // utils::trim(eta, etamin, etamax);
}

arma::mat get_eta (
    const arma::mat & offset,
    const arma::mat & u, const arma::mat & v, 
    const double & etamin, const double & etamax
) {
    arma::mat eta(u.n_rows, v.n_rows);
    set_eta(eta, offset, u, v, etamin, etamax);
    return eta;
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
    idu = join_cols(arma::linspace<arma::uvec>(p, p+q-1, q), 
                    arma::linspace<arma::uvec>(p+q, p+q+d-1, d));
    idv = join_cols(arma::linspace<arma::uvec>(0, p-1, p), 
                    arma::linspace<arma::uvec>(p+q, p+q+d-1, d));
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
        Rprintf(" %9i %11.2f %9.5f %8.2f s \n", iter, 100*dev, change, time);
    } else {
        Rprintf(" %9i %11.2f %9.5f %8.2f m \n", iter, 100*dev, change, time/60);
    }
}

void print_state (
    const int & iter, const double & dev, 
    const double & change, const double & time,
    const double & scanned
) {
    if (time < 60) {
        Rprintf(
            " %9i %11.2f %9.5f %7.0f %% %8.2f s \n", 
            iter, 100*dev, change, 100*scanned, time);
    } else {
        Rprintf(
            " %9i %11.2f %9.5f %7.0f %% %8.2f m \n", 
            iter, 100*dev, change, 100*scanned, time/60);
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