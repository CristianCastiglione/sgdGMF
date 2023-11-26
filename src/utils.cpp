// utils.cpp
// author: Cristian Castiglione
// creation: 28/09/2023
// last change: 07/10/2023

#include "utils.h"

namespace utils {

double absmax (const double & u, const double & v) {
  double a = std::abs(u - v);
  double b = std::max(std::abs(u), std::abs(v)) + 1e-04;
  return a / b;
}

double absmax (const arma::vec & u, const arma::vec & v) {
  double a = arma::max(arma::abs(u - v));
  double b = arma::max(arma::abs(arma::join_cols(u, v))) + 1e-04;
  return a / b;
}

double norm (const arma::mat & x) {
  return std::sqrt(arma::accu(arma::square(x)));
}

double norm (const arma::mat & x, const double & p) {
  return std::pow(arma::accu(arma::pow(x, p)), 1/p);
}

void trim (arma::mat & x, const double & a, const double & b) {
  arma::uvec above = arma::find(x > b);
  arma::uvec below = arma::find(x < a);
  x.elem(above).fill(b);
  x.elem(below).fill(a);
}

void trim (arma::mat & x, const double & a, const double & b, const arma::uvec & idx) {
  arma::mat xt = x.rows(idx);
  trim(xt, a, b);
  x.rows(idx) = xt;
}

void trim (arma::mat & x, const double & a, const double & b, const arma::uvec & idx, const arma::uvec & idy) {
  arma::mat xt = x(idx, idy);
  trim(xt, a, b);
  x(idx, idy) = xt;
}

arma::mat max0 (const arma::mat & x) {
  return .5 * (arma::abs(x) + x);
}

arma::mat max0 (const arma::mat & x, const double & p) {
  return arma::pow(.5 * (arma::abs(x) + x), p);
}

arma::mat xlogx (const arma::mat & x){
  arma::mat y = x;
  arma::uvec above = arma::find(x >  0);
  arma::uvec below = arma::find(x <= 0);
  y.elem(above) = x.elem(above) % arma::log(x.elem(above));
  y.elem(below).fill(0);
  return y;
}

arma::mat log1pexp (const arma::mat & x) {
  bool stable = false;
  if (stable) {
    // Thresholds
    double c0 = -37.0, c1 = 18.0, c2 = 33.3;
    
    // Test conditions
    arma::umat tst;
    arma::umat tst0 = (x >  c0);
    arma::umat tst1 = (x <= c1);
    arma::umat tst2 = (x <= c2);

    // Index matrix
    arma::umat idx;

    // Output for x <= c0
    arma::mat out = arma::exp(x);
    
    // Output for c0 < x <= c1
    tst = tst0 % tst1;
    if (arma::any(arma::any(tst))) {
      idx = arma::find(tst);
      out.elem(idx) = arma::log1p(out.elem(idx));
    }

    // Output for c1 < x <= c2
    tst = (1 - tst1) % tst2;
    if (arma::any(arma::any(tst))) {
      idx = arma::find(tst);
      out.elem(idx) = x.elem(idx) + 1 / out.elem(idx);
    }

    // Output for x > c2
    tst = (1 - tst2);
    if (arma::any(arma::any(tst))) {
      idx = arma::find(tst);
      out.elem(idx) = x.elem(idx);
    }

    return (out);
  } else {
    return arma::log1p(arma::exp(x));
  }
}

arma::mat log1mexp (const arma::mat & x) {
  arma::umat tst = (x > std::log(2.0));
  arma::vec out = x;
  arma::umat idx;
  idx = arma::find(tst);
  out.elem(idx) = arma::log(-arma::expm1(-x.elem(idx)));
  idx = arma::find(1-tst);
  out.elem(idx) = arma::log1p(-arma::exp(-x.elem(idx)));
  return out;
}

arma::mat logit (const arma::mat & x) {
  return arma::log(x) - arma::log(1 - x);
}

arma::mat expit (const arma::mat & x) {
  return arma::exp(x - arma::log1p(arma::exp(x)));
}

arma::mat expit2 (const arma::mat & x) {
  return arma::exp(x - 2.0 * arma::log1p(arma::exp(x)));
}

arma::mat expitn (const arma::mat & x, const double & n) {
  return arma::exp(x - n * arma::log1p(arma::exp(x)));
}

arma::mat cloglog (const arma::mat & x) {
  return arma::log(- arma::log(1 - x));
}

arma::mat cexpexp (const arma::mat & x) {
  return 1 - arma::exp(- arma::exp(x));
}

arma::mat loglog (const arma::mat & x) {
  return - arma::log(- arma::log(x));
}

arma::mat expexp (const arma::mat & x) {
  return arma::exp(- arma::exp(- x));
}

arma::mat pdfn (const arma::mat & x) {
  return arma::normpdf(x);
}

arma::mat cdfn (const arma::mat & x) {
  return arma::normcdf(x);
}

arma::mat logpdfn (const arma::mat & x) {
  return arma::log_normpdf(x);
}

arma::mat logcdfn (const arma::mat & x) {
  return arma::log(arma::normcdf(x));
}

arma::mat gamma (const arma::mat & x) {
  return arma::tgamma(x);
}

arma::mat loggamma (const arma::mat & x) {
  return arma::lgamma(x);
}

arma::mat digamma (const arma::mat & x) {
  int n = x.n_rows;
  int m = x.n_cols;
  arma::vec coef = {-2, -12, +120, -252, +240, -132, +32760/691, -12};
  arma::mat psi(n, m, arma::fill::zeros);
  if (x.min() > 6) {
    psi  = arma::log(x);
    psi += 1 / (coef(0) * arma::pow(x, 1)) + 1 / (coef(1) * arma::pow(x, 2)) + 1 / (coef(2) * arma::pow(x, 4));
    psi += 1 / (coef(3) * arma::pow(x, 6)) + 1 / (coef(4) * arma::pow(x, 8)) + 1 / (coef(5) * arma::pow(x, 10));
    // psi += 1 / (coef(6) * arma::pow(x, 12)) + 1 / (coef(7) * arma::pow(x, 14));
  } else {
    arma::mat z = x + 6;
    psi  = arma::log(z);
    psi += 1 / (coef(0) * arma::pow(z, 1)) + 1 / (coef(1) * arma::pow(z, 2)) + 1 / (coef(2) * arma::pow(z, 4));
    psi += 1 / (coef(3) * arma::pow(z, 6)) + 1 / (coef(4) * arma::pow(z, 8)) + 1 / (coef(5) * arma::pow(z, 10));
    // psi += 1 / (coef(6) * arma::pow(z, 12)) + 1 / (coef(7) * arma::pow(z, 14));
    psi -= 1 / x + 1 / (x + 1) + 1 / (x + 2) + 1 / (x + 3) + 1 / (x + 4) + 1 / (x + 5);
  }
  return psi;
}

arma::mat trigamma (const arma::mat & x) {
  int n = x.n_rows;
  int m = x.n_cols;
  arma::vec coef = {+1, +2, +6, -30, +42, -30, +66/5, -2730/691, +6/7};
  arma::mat psi(n, m, arma::fill::zeros);
  if (x.min() > 6) {
    psi += 1 / (coef(0) * arma::pow(x, 1)) + 1 / (coef(1) * arma::pow(x, 2)) + 1 / (coef(2) * arma::pow(x, 3));
    psi += 1 / (coef(3) * arma::pow(x, 5)) + 1 / (coef(4) * arma::pow(x, 7)) + 1 / (coef(5) * arma::pow(x, 9));
    // psi += 1 / (coef(6) * arma::pow(x, 11)) + 1 / (coef(7) * arma::pow(x, 13));
  } else {
    arma::mat z = x + 6;
    psi += 1 / (coef(0) * arma::pow(z, 1)) + 1 / (coef(1) * arma::pow(z, 2)) + 1 / (coef(2) * arma::pow(z, 3));
    psi += 1 / (coef(3) * arma::pow(z, 5)) + 1 / (coef(4) * arma::pow(z, 7)) + 1 / (coef(5) * arma::pow(z, 9));
    // psi += 1 / (coef(6) * arma::pow(z, 11)) + 1 / (coef(7) * arma::pow(z, 13));
    psi += 1 / arma::square(x) + 1 / arma::square(x + 1) + 1 / arma::square(x + 2);
    psi += 1 / arma::square(x + 3) + 1 / arma::square(x + 4) + 1 / arma::square(x + 5);
  }
  return psi;
}

arma::mat beta (const arma::mat & x, const arma::mat & y) {
  return arma::tgamma(x) % arma::tgamma(y) / arma::tgamma(x + y);
}

arma::mat logbeta (const arma::mat & x, const arma::mat & y) {
  return arma::lgamma(x) + arma::lgamma(y) - arma::lgamma(x + y);
}

arma::mat dibeta (const arma::mat & x, const arma::mat & y) {
  return digamma(x) + digamma(y) - digamma(x + y);
}

arma::mat tribeta (const arma::mat & x, const arma::mat & y) {
  return trigamma(x) + trigamma(y) - trigamma(x + y);
}

arma::mat hinge (const arma::mat & x) {
  return 0.5 * (arma::abs(x) + x);
}

arma::mat dirac (const arma::mat & x, const double & a = 0) {
  return arma::conv_to<arma::mat>::from(x == a);
}

arma::mat step (const arma::mat & x, const double & a = 0, const bool & lower = true) {
  arma::mat s(arma::size(x));
  if (lower) {
    s = arma::conv_to<arma::mat>::from(x < a);
  } else {
    s = arma::conv_to<arma::mat>::from(x > a);
  }
  return s;
}

arma::vec vech(const arma::mat & A) {
  int n = A.n_rows;
  int m = static_cast<int>(0.5 * n * (n + 1));
  int k = 0;
  
  arma::vec a(m);
  
  for (int i = 0; i < n; i++) {
    k += i;
    for (int j = 0; j <= i; j++) {
      a(k+j) = A(i,j);
    }
  }
  
  return a;
}

}
