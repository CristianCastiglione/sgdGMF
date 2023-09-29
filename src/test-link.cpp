// test-link.cpp
// author: Cristian Castiglione
// creation: 29/09/2023
// last change: 29/09/2023

#include "link.h"


// [[Rcpp::export]]
arma::vec c_link_identity_linkfun (const arma::vec & mu) {Link::Identity link; return link.linkfun(mu);}
// [[Rcpp::export]]
arma::vec c_link_identity_linkinv (const arma::vec & eta) {Link::Identity link; return link.linkinv(eta);}
// [[Rcpp::export]]
arma::vec c_link_identity_mueta (const arma::vec & eta) {Link::Identity link; return link.mueta(eta);}

// [[Rcpp::export]]
arma::vec c_link_logit_linkfun (const arma::vec & mu) {Link::Logit link; return link.linkfun(mu);}
// [[Rcpp::export]]
arma::vec c_link_logit_linkinv (const arma::vec & eta) {Link::Logit link; return link.linkinv(eta);}
// [[Rcpp::export]]
arma::vec c_link_logit_mueta (const arma::vec & eta) {Link::Logit link; return link.mueta(eta);}

// [[Rcpp::export]]
arma::vec c_link_probit_linkfun (const arma::vec & mu) {Link::Probit link; return link.linkfun(mu);}
// [[Rcpp::export]]
arma::vec c_link_probit_linkinv (const arma::vec & eta) {Link::Probit link; return link.linkinv(eta);}
// [[Rcpp::export]]
arma::vec c_link_probit_mueta (const arma::vec & eta) {Link::Probit link; return link.mueta(eta);}

// [[Rcpp::export]]
arma::vec c_link_cauchy_linkfun (const arma::vec & mu) {Link::Cauchy link; return link.linkfun(mu);}
// [[Rcpp::export]]
arma::vec c_link_cauchy_linkinv (const arma::vec & eta) {Link::Cauchy link; return link.linkinv(eta);}
// [[Rcpp::export]]
arma::vec c_link_cauchy_mueta (const arma::vec & eta) {Link::Cauchy link; return link.mueta(eta);}

// [[Rcpp::export]]
arma::vec c_link_cloglog_linkfun (const arma::vec & mu) {Link::cLogLog link; return link.linkfun(mu);}
// [[Rcpp::export]]
arma::vec c_link_cloglog_linkinv (const arma::vec & eta) {Link::cLogLog link; return link.linkinv(eta);}
// [[Rcpp::export]]
arma::vec c_link_cloglog_mueta (const arma::vec & eta) {Link::cLogLog link; return link.mueta(eta);}

// [[Rcpp::export]]
arma::vec c_link_log_linkfun (const arma::vec & mu) {Link::Log link; return link.linkfun(mu);}
// [[Rcpp::export]]
arma::vec c_link_log_linkinv (const arma::vec & eta) {Link::Log link; return link.linkinv(eta);}
// [[Rcpp::export]]
arma::vec c_link_log_mueta (const arma::vec & eta) {Link::Log link; return link.mueta(eta);}

// [[Rcpp::export]]
arma::vec c_link_inverse_linkfun (const arma::vec & mu) {Link::Inverse link; return link.linkfun(mu);}
// [[Rcpp::export]]
arma::vec c_link_inverse_linkinv (const arma::vec & eta) {Link::Inverse link; return link.linkinv(eta);}
// [[Rcpp::export]]
arma::vec c_link_inverse_mueta (const arma::vec & eta) {Link::Inverse link; return link.mueta(eta);}

// [[Rcpp::export]]
arma::vec c_link_sqrt_linkfun (const arma::vec & mu) {Link::Sqrt link; return link.linkfun(mu);}
// [[Rcpp::export]]
arma::vec c_link_sqrt_linkinv (const arma::vec & eta) {Link::Sqrt link; return link.linkinv(eta);}
// [[Rcpp::export]]
arma::vec c_link_sqrt_mueta (const arma::vec & eta) {Link::Sqrt link; return link.mueta(eta);}
