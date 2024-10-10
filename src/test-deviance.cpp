// test-deviance.h
// author: Cristian Castiglione
// creation: 02/10/2023
// last change: 02/10/2023

#include "deviance.h"
#include "misc.h"
#include <memory>

using namespace glm;

//' @keywords internal
// [[Rcpp::export("cpp.deviance")]]
arma::mat cpp_deviance (const arma::mat & y, const arma::mat & mu, const std::string & familyname) {
    std::unique_ptr<Family> family = make_family(familyname, std::string("identity"));
    return deviance(y, mu, family);
}

//' @keywords internal
// [[Rcpp::export("cpp.penalty")]]
double cpp_penalty (const arma::mat & u, const arma::vec & p) {
    return penalty(u, p);
}
