// test-link.cpp
// author: Cristian Castiglione
// creation: 29/09/2023
// last change: 29/09/2023

#include "link.h"

using namespace glm;

//' @keywords internal
// [[Rcpp::export("cpp.link.identity.linkfun")]]
arma::vec cpp_link_identity_linkfun (const arma::vec & mu) {Identity link; return link.linkfun(mu);}
//' @keywords internal
// [[Rcpp::export("cpp.link.identity.linkinv")]]
arma::vec cpp_link_identity_linkinv (const arma::vec & eta) {Identity link; return link.linkinv(eta);}
//' @keywords internal
// [[Rcpp::export("cpp.link.identity.mueta")]]
arma::vec cpp_link_identity_mueta (const arma::vec & eta) {Identity link; return link.mueta(eta);}

//' @keywords internal
// [[Rcpp::export("cpp.link.logit.linkfun")]]
arma::vec cpp_link_logit_linkfun (const arma::vec & mu) {Logit link; return link.linkfun(mu);}
//' @keywords internal
// [[Rcpp::export("cpp.link.logit.linkinv")]]
arma::vec cpp_link_logit_linkinv (const arma::vec & eta) {Logit link; return link.linkinv(eta);}
//' @keywords internal
// [[Rcpp::export("cpp.link.logit.mueta")]]
arma::vec cpp_link_logit_mueta (const arma::vec & eta) {Logit link; return link.mueta(eta);}

//' @keywords internal
// [[Rcpp::export("cpp.link.probit.linkfun")]]
arma::vec cpp_link_probit_linkfun (const arma::vec & mu) {Probit link; return link.linkfun(mu);}
//' @keywords internal
// [[Rcpp::export("cpp.link.probit.linkinv")]]
arma::vec cpp_link_probit_linkinv (const arma::vec & eta) {Probit link; return link.linkinv(eta);}
//' @keywords internal
// [[Rcpp::export("cpp.link.probit.mueta")]]
arma::vec cpp_link_probit_mueta (const arma::vec & eta) {Probit link; return link.mueta(eta);}

//' @keywords internal
// [[Rcpp::export("cpp.link.cauchy.linkfun")]]
arma::vec cpp_link_cauchy_linkfun (const arma::vec & mu) {Cauchy link; return link.linkfun(mu);}
//' @keywords internal
// [[Rcpp::export("cpp.link.cauchy.linkinv")]]
arma::vec cpp_link_cauchy_linkinv (const arma::vec & eta) {Cauchy link; return link.linkinv(eta);}
//' @keywords internal
// [[Rcpp::export("cpp.link.cauchy.mueta")]]
arma::vec cpp_link_cauchy_mueta (const arma::vec & eta) {Cauchy link; return link.mueta(eta);}

//' @keywords internal
// [[Rcpp::export("cpp.link.cloglog.linkfun")]]
arma::vec cpp_link_cloglog_linkfun (const arma::vec & mu) {cLogLog link; return link.linkfun(mu);}
//' @keywords internal
// [[Rcpp::export("cpp.link.cloglog.linkinv")]]
arma::vec cpp_link_cloglog_linkinv (const arma::vec & eta) {cLogLog link; return link.linkinv(eta);}
//' @keywords internal
// [[Rcpp::export("cpp.link.cloglog.mueta")]]
arma::vec cpp_link_cloglog_mueta (const arma::vec & eta) {cLogLog link; return link.mueta(eta);}

//' @keywords internal
// [[Rcpp::export("cpp.link.log.linkfun")]]
arma::vec cpp_link_log_linkfun (const arma::vec & mu) {Log link; return link.linkfun(mu);}
//' @keywords internal
// [[Rcpp::export("cpp.link.log.linkinv")]]
arma::vec cpp_link_log_linkinv (const arma::vec & eta) {Log link; return link.linkinv(eta);}
//' @keywords internal
// [[Rcpp::export("cpp.link.log.mueta")]]
arma::vec cpp_link_log_mueta (const arma::vec & eta) {Log link; return link.mueta(eta);}

//' @keywords internal
// [[Rcpp::export("cpp.link.inverse.linkfun")]]
arma::vec cpp_link_inverse_linkfun (const arma::vec & mu) {Inverse link; return link.linkfun(mu);}
//' @keywords internal
// [[Rcpp::export("cpp.link.inverse.linkinv")]]
arma::vec cpp_link_inverse_linkinv (const arma::vec & eta) {Inverse link; return link.linkinv(eta);}
//' @keywords internal
// [[Rcpp::export("cpp.link.inverse.mueta")]]
arma::vec cpp_link_inverse_mueta (const arma::vec & eta) {Inverse link; return link.mueta(eta);}

//' @keywords internal
// [[Rcpp::export("cpp.link.sqrt.linkfun")]]
arma::vec cpp_link_sqrt_linkfun (const arma::vec & mu) {Sqrt link; return link.linkfun(mu);}
//' @keywords internal
// [[Rcpp::export("cpp.link.sqrt.linkinv")]]
arma::vec cpp_link_sqrt_linkinv (const arma::vec & eta) {Sqrt link; return link.linkinv(eta);}
//' @keywords internal
// [[Rcpp::export("cpp.link.sqrt.mueta")]]
arma::vec cpp_link_sqrt_mueta (const arma::vec & eta) {Sqrt link; return link.mueta(eta);}
