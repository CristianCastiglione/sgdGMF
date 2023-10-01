// rcppglm.cpp
// author: Cristian Castiglione
// creation: 30/09/2023
// last change: 30/09/2023

#include "rcppglm.h"

// std::unique_ptr<> is a smart pointer managing exclusice ownership and dynamic object allocation
// std::make_unique<> is a helper function returning manaaging object allocation and object construction
// they both require the inclusion of the <memory> library 
std::unique_ptr<Link::Link> make_link (const std::string & linkname) {
    std::unique_ptr<Link::Link> ptr;
    if (linkname == "identity") { ptr = std::make_unique<Link::Identity>();
    } else if (linkname == "logit") { ptr = std::make_unique<Link::Logit>();
    } else if (linkname == "probit") { ptr = std::make_unique<Link::Probit>();
    } else if (linkname == "cauchy") { ptr = std::make_unique<Link::Cauchy>();
    } else if (linkname == "cloglog") { ptr = std::make_unique<Link::cLogLog>();
    } else if (linkname == "log") { ptr = std::make_unique<Link::Log>();
    } else if (linkname == "inverse") { ptr = std::make_unique<Link::Inverse>();
    } else if (linkname == "sqrt") { ptr = std::make_unique<Link::Sqrt>();
    } else { Rcpp::stop("Link function not available."); }
    return ptr;
}

template<class F>
inline Rcpp::XPtr<F> make_family (std::string linkname) {
    auto link = make_link(linkname);
    auto* family = new F(link);
    Rcpp::XPtr<F> pointer(family, true);
    return pointer;
}

Rcpp::XPtr<Family::Gaussian> make_gaussian (std::string linkname) {
    return make_family<Family::Gaussian>(linkname);
}

Rcpp::XPtr<Family::Binomial> make_binomial (std::string linkname) {
    return make_family<Family::Binomial>(linkname);
}

Rcpp::XPtr<Family::Poisson> make_poisson (std::string linkname) {
    return make_family<Family::Poisson>(linkname);
}

Rcpp::XPtr<Family::Gamma> make_gamma (std::string linkname) {
    return make_family<Family::Gamma>(linkname);
}

