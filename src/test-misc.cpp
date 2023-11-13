// test-misc.cpp
// author: Cristian Castiglione
// creation: 01/10/2023
// last change: 13/10/2023

#include "misc.h"

using namespace glm;

void c_print_link_family (const std::unique_ptr<Family> & family) {
    Rcpp::Rcout << "Family: " << family->getfamily() << "\n";
    Rcpp::Rcout << "Link: " << family->getlink() << "\n";
    Rcpp::Rcout << "Mu: " << arma::vec{0.25, 0.5, 0.75} << "\n";
    Rcpp::Rcout << "Eta: " << family->linkfun(arma::vec{0.25, 0.5, 0.75}) << "\n";
}

// [[Rcpp::export]]
void c_make_link_family (const std::string & familyname, const std::string & linkname) {
    std::unique_ptr<Family> family = make_family(familyname, linkname);
    c_print_link_family(family);
}

// [[Rcpp::export]]
Rcpp::List c_get_data_bounds (
    const double & eps, const double & ymin, const double & ymax, 
    const std::string & familyname, const std::string & linkname
) {
    std::unique_ptr<Family> family = make_family(familyname, linkname);

    double mulo, muup, etalo, etaup;
    set_data_bounds(mulo, muup, etalo, etaup, eps, ymin, ymax, family);

    Rcpp::List out;
    out["family"] = family->getfamily();
    out["link"] = family->getlink();
    out["ylim"] = arma::vec{ymin, ymax};
    out["mulim"] = arma::vec{mulo, muup};
    out["etalim"] = arma::vec{etalo, etaup};

    return out;
}

// [[Rcpp::export]]
Rcpp::List c_get_uv_penalty (
    const arma::vec & pen, 
    const int & p, const int & q, const int & d
) {
    arma::vec penu(p+q+d), penv(p+q+d);
    set_uv_penalty(penu, penv, pen, p, q, d);
    
    Rcpp::List out;
    out["penu"] = penu;
    out["penv"] = penv;

    return out;
}

// [[Rcpp::export]]
Rcpp::List c_get_uv_indices (
    const int & p, const int & q, const int & d
) {
    arma::uvec idu, idv;
    set_uv_indices(idu, idv, p, q, d);

    Rcpp::List out;
    out["idu"] = idu;
    out["idv"] = idv;

    return out;
}

// [[Rcpp::export]]
std::list<arma::uvec> c_sample_minibatch (
    const int & n, const int & size, const bool & randomize
) {
    return sample_chunks(n, size, randomize);
}

// [[Rcpp::export]]
int c_select_minibatch (const int & iter, const int & nchunks) {
    return select_chunk(iter, nchunks);
}