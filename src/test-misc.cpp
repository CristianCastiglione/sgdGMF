// test-misc.cpp
// author: Cristian Castiglione
// creation: 01/10/2023
// last change: 01/10/2023

#include "misc.h"


// template<class F, class L>
// void c_print_link_family (const L & link, const F & family) {
//     Rcpp::Rcout << "Family: " << family->family << "\n";
//     Rcpp::Rcout << "Link: " << link->link << "\n";
// }

void c_print_link_family (
    const std::unique_ptr<Link::Link> & link, 
    const std::unique_ptr<Family::Family> & family
) {
    Rcpp::Rcout << "Family: " << family->family << "\n";
    Rcpp::Rcout << "Link: " << link->link << "\n";
    Rcpp::Rcout << "Mu: " << arma::vec{0.25, 0.5, 0.75} << "\n";
    Rcpp::Rcout << "Eta: " << link->linkfun(arma::vec{0.25, 0.5, 0.75}) << "\n";
}

// [[Rcpp::export]]
void c_make_link_family (
    const std::string & linkname, const std::string & familyname
) {
    std::unique_ptr<Link::Link> link = make_link(linkname);
    std::unique_ptr<Family::Family> family = make_family(familyname);
    c_print_link_family(link, family);
}

// [[Rcpp::export]]
Rcpp::List c_get_data_bounds (
    const double & eps, const double & ymin, const double & ymax, 
    const std::string & familyname, const std::string & linkname
) {
    std::unique_ptr<Family::Family> family = make_family(familyname);
    std::unique_ptr<Link::Link> link = make_link(linkname);

    double mulo, muup, etalo, etaup;
    set_data_bounds(mulo, muup, etalo, etaup, eps, ymin, ymax, family, link);

    // Rcpp::Rcout << "Family: " << family->family << "\n";
    // Rcpp::Rcout << "Link: " << link->link << "\n";
    // Rcpp::Rcout << "Y: [" << ymin << ", " << ymax << "] \n";
    // Rcpp::Rcout << "Mu: [" << mulo << ", " << muup << "] \n";
    // Rcpp::Rcout << "Eta: [" << etalo << ", " << etaup << "] \n";

    Rcpp::List out;
    out["family"] = family->family;
    out["link"] = link->link;
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
    // arma::vec pen = {1,2,3,4};
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