// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// c_gaussian_variance
arma::vec c_gaussian_variance(const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_gaussian_variance(SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_gaussian_variance(mu));
    return rcpp_result_gen;
END_RCPP
}
// c_gaussian_initialize
arma::vec c_gaussian_initialize(const arma::vec& y);
RcppExport SEXP _sgdGMF_c_gaussian_initialize(SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(c_gaussian_initialize(y));
    return rcpp_result_gen;
END_RCPP
}
// c_gaussian_devresid
arma::vec c_gaussian_devresid(const arma::vec& y, const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_gaussian_devresid(SEXP ySEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_gaussian_devresid(y, mu));
    return rcpp_result_gen;
END_RCPP
}
// c_binomial_variance
arma::vec c_binomial_variance(const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_binomial_variance(SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_binomial_variance(mu));
    return rcpp_result_gen;
END_RCPP
}
// c_binomial_initialize
arma::vec c_binomial_initialize(const arma::vec& y);
RcppExport SEXP _sgdGMF_c_binomial_initialize(SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(c_binomial_initialize(y));
    return rcpp_result_gen;
END_RCPP
}
// c_binomial_devresid
arma::vec c_binomial_devresid(const arma::vec& y, const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_binomial_devresid(SEXP ySEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_binomial_devresid(y, mu));
    return rcpp_result_gen;
END_RCPP
}
// c_poisson_variance
arma::vec c_poisson_variance(const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_poisson_variance(SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_poisson_variance(mu));
    return rcpp_result_gen;
END_RCPP
}
// c_poisson_initialize
arma::vec c_poisson_initialize(const arma::vec& y);
RcppExport SEXP _sgdGMF_c_poisson_initialize(SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(c_poisson_initialize(y));
    return rcpp_result_gen;
END_RCPP
}
// c_poisson_devresid
arma::vec c_poisson_devresid(const arma::vec& y, const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_poisson_devresid(SEXP ySEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_poisson_devresid(y, mu));
    return rcpp_result_gen;
END_RCPP
}
// c_gamma_variance
arma::vec c_gamma_variance(const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_gamma_variance(SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_gamma_variance(mu));
    return rcpp_result_gen;
END_RCPP
}
// c_gamma_initialize
arma::vec c_gamma_initialize(const arma::vec& y);
RcppExport SEXP _sgdGMF_c_gamma_initialize(SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(c_gamma_initialize(y));
    return rcpp_result_gen;
END_RCPP
}
// c_gamma_devresid
arma::vec c_gamma_devresid(const arma::vec& y, const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_gamma_devresid(SEXP ySEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_gamma_devresid(y, mu));
    return rcpp_result_gen;
END_RCPP
}
// c_binomial_logit_loglik
arma::vec c_binomial_logit_loglik(const arma::vec& y, const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_binomial_logit_loglik(SEXP ySEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_binomial_logit_loglik(y, eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_identity_linkfun
arma::vec c_link_identity_linkfun(const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_link_identity_linkfun(SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_identity_linkfun(mu));
    return rcpp_result_gen;
END_RCPP
}
// c_link_identity_linkinv
arma::vec c_link_identity_linkinv(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_identity_linkinv(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_identity_linkinv(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_identity_mueta
arma::vec c_link_identity_mueta(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_identity_mueta(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_identity_mueta(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_logit_linkfun
arma::vec c_link_logit_linkfun(const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_link_logit_linkfun(SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_logit_linkfun(mu));
    return rcpp_result_gen;
END_RCPP
}
// c_link_logit_linkinv
arma::vec c_link_logit_linkinv(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_logit_linkinv(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_logit_linkinv(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_logit_mueta
arma::vec c_link_logit_mueta(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_logit_mueta(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_logit_mueta(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_probit_linkfun
arma::vec c_link_probit_linkfun(const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_link_probit_linkfun(SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_probit_linkfun(mu));
    return rcpp_result_gen;
END_RCPP
}
// c_link_probit_linkinv
arma::vec c_link_probit_linkinv(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_probit_linkinv(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_probit_linkinv(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_probit_mueta
arma::vec c_link_probit_mueta(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_probit_mueta(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_probit_mueta(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_cauchy_linkfun
arma::vec c_link_cauchy_linkfun(const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_link_cauchy_linkfun(SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_cauchy_linkfun(mu));
    return rcpp_result_gen;
END_RCPP
}
// c_link_cauchy_linkinv
arma::vec c_link_cauchy_linkinv(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_cauchy_linkinv(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_cauchy_linkinv(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_cauchy_mueta
arma::vec c_link_cauchy_mueta(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_cauchy_mueta(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_cauchy_mueta(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_cloglog_linkfun
arma::vec c_link_cloglog_linkfun(const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_link_cloglog_linkfun(SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_cloglog_linkfun(mu));
    return rcpp_result_gen;
END_RCPP
}
// c_link_cloglog_linkinv
arma::vec c_link_cloglog_linkinv(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_cloglog_linkinv(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_cloglog_linkinv(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_cloglog_mueta
arma::vec c_link_cloglog_mueta(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_cloglog_mueta(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_cloglog_mueta(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_log_linkfun
arma::vec c_link_log_linkfun(const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_link_log_linkfun(SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_log_linkfun(mu));
    return rcpp_result_gen;
END_RCPP
}
// c_link_log_linkinv
arma::vec c_link_log_linkinv(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_log_linkinv(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_log_linkinv(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_log_mueta
arma::vec c_link_log_mueta(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_log_mueta(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_log_mueta(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_inverse_linkfun
arma::vec c_link_inverse_linkfun(const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_link_inverse_linkfun(SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_inverse_linkfun(mu));
    return rcpp_result_gen;
END_RCPP
}
// c_link_inverse_linkinv
arma::vec c_link_inverse_linkinv(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_inverse_linkinv(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_inverse_linkinv(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_inverse_mueta
arma::vec c_link_inverse_mueta(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_inverse_mueta(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_inverse_mueta(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_sqrt_linkfun
arma::vec c_link_sqrt_linkfun(const arma::vec& mu);
RcppExport SEXP _sgdGMF_c_link_sqrt_linkfun(SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_sqrt_linkfun(mu));
    return rcpp_result_gen;
END_RCPP
}
// c_link_sqrt_linkinv
arma::vec c_link_sqrt_linkinv(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_sqrt_linkinv(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_sqrt_linkinv(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_link_sqrt_mueta
arma::vec c_link_sqrt_mueta(const arma::vec& eta);
RcppExport SEXP _sgdGMF_c_link_sqrt_mueta(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_link_sqrt_mueta(eta));
    return rcpp_result_gen;
END_RCPP
}
// c_dabsmax
double c_dabsmax(const double& u, const double& v);
RcppExport SEXP _sgdGMF_c_dabsmax(SEXP uSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dabsmax(u, v));
    return rcpp_result_gen;
END_RCPP
}
// c_vabsmax
double c_vabsmax(const arma::vec& u, const arma::vec& v);
RcppExport SEXP _sgdGMF_c_vabsmax(SEXP uSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(c_vabsmax(u, v));
    return rcpp_result_gen;
END_RCPP
}
// c_trim
arma::vec c_trim(const arma::vec& x, double a, double b);
RcppExport SEXP _sgdGMF_c_trim(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(c_trim(x, a, b));
    return rcpp_result_gen;
END_RCPP
}
// c_log1pexp
arma::vec c_log1pexp(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_log1pexp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_log1pexp(x));
    return rcpp_result_gen;
END_RCPP
}
// c_log1mexp
arma::vec c_log1mexp(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_log1mexp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_log1mexp(x));
    return rcpp_result_gen;
END_RCPP
}
// c_logit
arma::vec c_logit(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_logit(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_logit(x));
    return rcpp_result_gen;
END_RCPP
}
// c_expit
arma::vec c_expit(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_expit(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_expit(x));
    return rcpp_result_gen;
END_RCPP
}
// c_expit2
arma::vec c_expit2(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_expit2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_expit2(x));
    return rcpp_result_gen;
END_RCPP
}
// c_expitn
arma::vec c_expitn(const arma::vec& x, double n);
RcppExport SEXP _sgdGMF_c_expitn(SEXP xSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(c_expitn(x, n));
    return rcpp_result_gen;
END_RCPP
}
// c_cloglog
arma::vec c_cloglog(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_cloglog(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_cloglog(x));
    return rcpp_result_gen;
END_RCPP
}
// c_cexpexp
arma::vec c_cexpexp(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_cexpexp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_cexpexp(x));
    return rcpp_result_gen;
END_RCPP
}
// c_loglog
arma::vec c_loglog(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_loglog(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_loglog(x));
    return rcpp_result_gen;
END_RCPP
}
// c_expexp
arma::vec c_expexp(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_expexp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_expexp(x));
    return rcpp_result_gen;
END_RCPP
}
// c_pdfn
arma::vec c_pdfn(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_pdfn(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_pdfn(x));
    return rcpp_result_gen;
END_RCPP
}
// c_cdfn
arma::vec c_cdfn(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_cdfn(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_cdfn(x));
    return rcpp_result_gen;
END_RCPP
}
// c_logpdfn
arma::vec c_logpdfn(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_logpdfn(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_logpdfn(x));
    return rcpp_result_gen;
END_RCPP
}
// c_logcdfn
arma::vec c_logcdfn(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_logcdfn(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_logcdfn(x));
    return rcpp_result_gen;
END_RCPP
}
// c_gamma
arma::vec c_gamma(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_gamma(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_gamma(x));
    return rcpp_result_gen;
END_RCPP
}
// c_loggamma
arma::vec c_loggamma(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_loggamma(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_loggamma(x));
    return rcpp_result_gen;
END_RCPP
}
// c_digamma
arma::vec c_digamma(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_digamma(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_digamma(x));
    return rcpp_result_gen;
END_RCPP
}
// c_trigamma
arma::vec c_trigamma(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_trigamma(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_trigamma(x));
    return rcpp_result_gen;
END_RCPP
}
// c_beta
arma::vec c_beta(const arma::vec& x, const arma::vec& y);
RcppExport SEXP _sgdGMF_c_beta(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(c_beta(x, y));
    return rcpp_result_gen;
END_RCPP
}
// c_logbeta
arma::vec c_logbeta(const arma::vec& x, const arma::vec& y);
RcppExport SEXP _sgdGMF_c_logbeta(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(c_logbeta(x, y));
    return rcpp_result_gen;
END_RCPP
}
// c_dibeta
arma::vec c_dibeta(const arma::vec& x, const arma::vec& y);
RcppExport SEXP _sgdGMF_c_dibeta(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(c_dibeta(x, y));
    return rcpp_result_gen;
END_RCPP
}
// c_tribeta
arma::vec c_tribeta(const arma::vec& x, const arma::vec& y);
RcppExport SEXP _sgdGMF_c_tribeta(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(c_tribeta(x, y));
    return rcpp_result_gen;
END_RCPP
}
// c_hinge
arma::vec c_hinge(const arma::vec& x);
RcppExport SEXP _sgdGMF_c_hinge(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_hinge(x));
    return rcpp_result_gen;
END_RCPP
}
// c_dirac
arma::vec c_dirac(const arma::vec& x, double a);
RcppExport SEXP _sgdGMF_c_dirac(SEXP xSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dirac(x, a));
    return rcpp_result_gen;
END_RCPP
}
// c_step
arma::vec c_step(const arma::vec& x, double a, bool lower);
RcppExport SEXP _sgdGMF_c_step(SEXP xSEXP, SEXP aSEXP, SEXP lowerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< bool >::type lower(lowerSEXP);
    rcpp_result_gen = Rcpp::wrap(c_step(x, a, lower));
    return rcpp_result_gen;
END_RCPP
}
// c_vech
arma::vec c_vech(const arma::mat& A);
RcppExport SEXP _sgdGMF_c_vech(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(c_vech(A));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sgdGMF_c_gaussian_variance", (DL_FUNC) &_sgdGMF_c_gaussian_variance, 1},
    {"_sgdGMF_c_gaussian_initialize", (DL_FUNC) &_sgdGMF_c_gaussian_initialize, 1},
    {"_sgdGMF_c_gaussian_devresid", (DL_FUNC) &_sgdGMF_c_gaussian_devresid, 2},
    {"_sgdGMF_c_binomial_variance", (DL_FUNC) &_sgdGMF_c_binomial_variance, 1},
    {"_sgdGMF_c_binomial_initialize", (DL_FUNC) &_sgdGMF_c_binomial_initialize, 1},
    {"_sgdGMF_c_binomial_devresid", (DL_FUNC) &_sgdGMF_c_binomial_devresid, 2},
    {"_sgdGMF_c_poisson_variance", (DL_FUNC) &_sgdGMF_c_poisson_variance, 1},
    {"_sgdGMF_c_poisson_initialize", (DL_FUNC) &_sgdGMF_c_poisson_initialize, 1},
    {"_sgdGMF_c_poisson_devresid", (DL_FUNC) &_sgdGMF_c_poisson_devresid, 2},
    {"_sgdGMF_c_gamma_variance", (DL_FUNC) &_sgdGMF_c_gamma_variance, 1},
    {"_sgdGMF_c_gamma_initialize", (DL_FUNC) &_sgdGMF_c_gamma_initialize, 1},
    {"_sgdGMF_c_gamma_devresid", (DL_FUNC) &_sgdGMF_c_gamma_devresid, 2},
    {"_sgdGMF_c_binomial_logit_loglik", (DL_FUNC) &_sgdGMF_c_binomial_logit_loglik, 2},
    {"_sgdGMF_c_link_identity_linkfun", (DL_FUNC) &_sgdGMF_c_link_identity_linkfun, 1},
    {"_sgdGMF_c_link_identity_linkinv", (DL_FUNC) &_sgdGMF_c_link_identity_linkinv, 1},
    {"_sgdGMF_c_link_identity_mueta", (DL_FUNC) &_sgdGMF_c_link_identity_mueta, 1},
    {"_sgdGMF_c_link_logit_linkfun", (DL_FUNC) &_sgdGMF_c_link_logit_linkfun, 1},
    {"_sgdGMF_c_link_logit_linkinv", (DL_FUNC) &_sgdGMF_c_link_logit_linkinv, 1},
    {"_sgdGMF_c_link_logit_mueta", (DL_FUNC) &_sgdGMF_c_link_logit_mueta, 1},
    {"_sgdGMF_c_link_probit_linkfun", (DL_FUNC) &_sgdGMF_c_link_probit_linkfun, 1},
    {"_sgdGMF_c_link_probit_linkinv", (DL_FUNC) &_sgdGMF_c_link_probit_linkinv, 1},
    {"_sgdGMF_c_link_probit_mueta", (DL_FUNC) &_sgdGMF_c_link_probit_mueta, 1},
    {"_sgdGMF_c_link_cauchy_linkfun", (DL_FUNC) &_sgdGMF_c_link_cauchy_linkfun, 1},
    {"_sgdGMF_c_link_cauchy_linkinv", (DL_FUNC) &_sgdGMF_c_link_cauchy_linkinv, 1},
    {"_sgdGMF_c_link_cauchy_mueta", (DL_FUNC) &_sgdGMF_c_link_cauchy_mueta, 1},
    {"_sgdGMF_c_link_cloglog_linkfun", (DL_FUNC) &_sgdGMF_c_link_cloglog_linkfun, 1},
    {"_sgdGMF_c_link_cloglog_linkinv", (DL_FUNC) &_sgdGMF_c_link_cloglog_linkinv, 1},
    {"_sgdGMF_c_link_cloglog_mueta", (DL_FUNC) &_sgdGMF_c_link_cloglog_mueta, 1},
    {"_sgdGMF_c_link_log_linkfun", (DL_FUNC) &_sgdGMF_c_link_log_linkfun, 1},
    {"_sgdGMF_c_link_log_linkinv", (DL_FUNC) &_sgdGMF_c_link_log_linkinv, 1},
    {"_sgdGMF_c_link_log_mueta", (DL_FUNC) &_sgdGMF_c_link_log_mueta, 1},
    {"_sgdGMF_c_link_inverse_linkfun", (DL_FUNC) &_sgdGMF_c_link_inverse_linkfun, 1},
    {"_sgdGMF_c_link_inverse_linkinv", (DL_FUNC) &_sgdGMF_c_link_inverse_linkinv, 1},
    {"_sgdGMF_c_link_inverse_mueta", (DL_FUNC) &_sgdGMF_c_link_inverse_mueta, 1},
    {"_sgdGMF_c_link_sqrt_linkfun", (DL_FUNC) &_sgdGMF_c_link_sqrt_linkfun, 1},
    {"_sgdGMF_c_link_sqrt_linkinv", (DL_FUNC) &_sgdGMF_c_link_sqrt_linkinv, 1},
    {"_sgdGMF_c_link_sqrt_mueta", (DL_FUNC) &_sgdGMF_c_link_sqrt_mueta, 1},
    {"_sgdGMF_c_dabsmax", (DL_FUNC) &_sgdGMF_c_dabsmax, 2},
    {"_sgdGMF_c_vabsmax", (DL_FUNC) &_sgdGMF_c_vabsmax, 2},
    {"_sgdGMF_c_trim", (DL_FUNC) &_sgdGMF_c_trim, 3},
    {"_sgdGMF_c_log1pexp", (DL_FUNC) &_sgdGMF_c_log1pexp, 1},
    {"_sgdGMF_c_log1mexp", (DL_FUNC) &_sgdGMF_c_log1mexp, 1},
    {"_sgdGMF_c_logit", (DL_FUNC) &_sgdGMF_c_logit, 1},
    {"_sgdGMF_c_expit", (DL_FUNC) &_sgdGMF_c_expit, 1},
    {"_sgdGMF_c_expit2", (DL_FUNC) &_sgdGMF_c_expit2, 1},
    {"_sgdGMF_c_expitn", (DL_FUNC) &_sgdGMF_c_expitn, 2},
    {"_sgdGMF_c_cloglog", (DL_FUNC) &_sgdGMF_c_cloglog, 1},
    {"_sgdGMF_c_cexpexp", (DL_FUNC) &_sgdGMF_c_cexpexp, 1},
    {"_sgdGMF_c_loglog", (DL_FUNC) &_sgdGMF_c_loglog, 1},
    {"_sgdGMF_c_expexp", (DL_FUNC) &_sgdGMF_c_expexp, 1},
    {"_sgdGMF_c_pdfn", (DL_FUNC) &_sgdGMF_c_pdfn, 1},
    {"_sgdGMF_c_cdfn", (DL_FUNC) &_sgdGMF_c_cdfn, 1},
    {"_sgdGMF_c_logpdfn", (DL_FUNC) &_sgdGMF_c_logpdfn, 1},
    {"_sgdGMF_c_logcdfn", (DL_FUNC) &_sgdGMF_c_logcdfn, 1},
    {"_sgdGMF_c_gamma", (DL_FUNC) &_sgdGMF_c_gamma, 1},
    {"_sgdGMF_c_loggamma", (DL_FUNC) &_sgdGMF_c_loggamma, 1},
    {"_sgdGMF_c_digamma", (DL_FUNC) &_sgdGMF_c_digamma, 1},
    {"_sgdGMF_c_trigamma", (DL_FUNC) &_sgdGMF_c_trigamma, 1},
    {"_sgdGMF_c_beta", (DL_FUNC) &_sgdGMF_c_beta, 2},
    {"_sgdGMF_c_logbeta", (DL_FUNC) &_sgdGMF_c_logbeta, 2},
    {"_sgdGMF_c_dibeta", (DL_FUNC) &_sgdGMF_c_dibeta, 2},
    {"_sgdGMF_c_tribeta", (DL_FUNC) &_sgdGMF_c_tribeta, 2},
    {"_sgdGMF_c_hinge", (DL_FUNC) &_sgdGMF_c_hinge, 1},
    {"_sgdGMF_c_dirac", (DL_FUNC) &_sgdGMF_c_dirac, 2},
    {"_sgdGMF_c_step", (DL_FUNC) &_sgdGMF_c_step, 3},
    {"_sgdGMF_c_vech", (DL_FUNC) &_sgdGMF_c_vech, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_sgdGMF(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}