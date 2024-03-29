# test-utils.R
# author: Cristian Castiglione
# creation: 29/09/2023
# last change: 29/09/2023

## Workspace setup ----
rm(list = ls())
graphics.off()

# Package compilation and import
devtools::load_all()

## Test: absmax() ----
{
  du = runif(1); dv = runif(1)
  vu = runif(5); vv = runif(5)
  r_absmax = function (u, v) max(abs(u - v)) / (max(abs(c(u, v))) + 1e-04)
  print(all.equal(r_absmax(du, dv), sgdGMF::cpp.dabsmax(du, dv)))
  print(all.equal(r_absmax(vu, vv), sgdGMF::cpp.vabsmax(vu, vv)))
}


## Test: trim() ----
{
  x = rnorm(10); a = 0; b = 0.5
  r_trim = function (x, a, b) a * (x <= a) + b * (x >= b) + x * (x > a & x < b)
  print(all.equal(r_trim(x, a, b), drop(sgdGMF::cpp.trim(x, a, b))))
}

## Test: log1pexp() ----
{
  x = seq(from = -3, to = +3, length = 201)
  r_log1pexp = function (x) Rmpfr::log1pexp(x)
  print(all.equal(r_log1pexp(x), drop(sgdGMF::cpp.log1pexp(x))))
}

## Test: log1mexp() ----
{
  x = seq(from = 1e-08, to = 6, length = 201)
  r_log1mexp = function (x) Rmpfr::log1mexp(x)
  print(all.equal(r_log1mexp(x), drop(sgdGMF::cpp.log1mexp(x))))
}

## Test: logit() ----
{
  x = seq(from = 0.001, to = 0.999, length = 201)
  r_logit = function (x) log(x) - log(1 - x)
  print(all.equal(r_logit(x), drop(sgdGMF::cpp.logit(x))))
}

## Test: expit() ----
{
  x = seq(from = -3, to = +3, length = 201)
  print(all.equal(plogis(x), drop(sgdGMF::cpp.expit(x))))
}

## Test: expit2() ----
{
  x = seq(from = -3, to = +3, length = 201)
  print(all.equal(dlogis(x), drop(sgdGMF::cpp.expit2(x))))
}

## Test: expitn() ----
{
  x = seq(from = -3, to = +3, length = 201)
  print(all.equal(plogis(x), drop(sgdGMF::cpp.expitn(x, 1))))
  print(all.equal(dlogis(x), drop(sgdGMF::cpp.expitn(x, 2))))
  print(all.equal(exp(x) / (1 + exp(x))^3, drop(sgdGMF::cpp.expitn(x, 3))))
}

## Test: cloglog() ----
{
  x = seq(from = 0.001, to = 0.999, length = 201)
  r_cloglog = function (x) log(-log(1 - x))
  print(all.equal(r_cloglog(x), drop(sgdGMF::cpp.cloglog(x))))
}

## Test: cexpexp() ----
{
  x = seq(from = -3, to = +3, length = 201)
  r_cexpexp = function (x) 1 - exp(-exp(x))
  print(all.equal(r_cexpexp(x), drop(sgdGMF::cpp.cexpexp(x))))
}

## Test: loglog() ----
{
  x = seq(from = 0.001, to = 0.999, length = 201)
  r_loglog = function (x) - log(-log(x))
  print(all.equal(r_loglog(x), drop(sgdGMF::cpp.loglog(x))))
}

## Test: cexpexp() ----
{
  x = seq(from = -3, to = +3, length = 201)
  r_expexp = function (x) exp(-exp(-x))
  print(all.equal(r_expexp(x), drop(sgdGMF::cpp.expexp(x))))
}

## Test: pdfn(), cdfn() ----
{
  x = seq(from = -3, to = +3, length = 201)
  print(all.equal(dnorm(x), drop(sgdGMF::cpp.pdfn(x))))
  print(all.equal(pnorm(x), drop(sgdGMF::cpp.cdfn(x))))
}

## Test: logpdfn(), logcdf() ----
{
  x = seq(from = -3, to = +3, length = 201)
  print(all.equal(dnorm(x, log = TRUE), drop(sgdGMF::cpp.logpdfn(x))))
  print(all.equal(pnorm(x, log = TRUE), drop(sgdGMF::cpp.logcdfn(x))))
}

## Test: gamma(), loggamma(), digamma(), trigamma() ----
{
  x = seq(from = 0.001, to = +10, length = 201)
  print(all.equal(gamma(x), drop(sgdGMF::cpp.gamma(x))))
  print(all.equal(lgamma(x), drop(sgdGMF::cpp.loggamma(x))))
  print(all.equal(digamma(x), drop(sgdGMF::cpp.digamma(x))))
  print(all.equal(trigamma(x), drop(sgdGMF::cpp.trigamma(x))))
}

## Test: beta(), logbeta() ----
{
  x = rexp(200)
  y = rexp(200)
  print(all.equal(beta(x, y), drop(sgdGMF::cpp.beta(x, y))))
  print(all.equal(lbeta(x, y), drop(sgdGMF::cpp.logbeta(x, y))))
}

## Test: hinge() ----
{
  x = seq(from = -3, to = +3, length = 201)
  r_hinge = function (x) 0.5 * (abs(x) + x)
  print(all.equal(r_hinge(x), drop(sgdGMF::cpp.hinge(x))))
}

## Test: dirac() ----
{
  x = seq(from = -3, to = +3, length = 201)
  print(all.equal(1 * (x == 0), drop(sgdGMF::cpp.dirac(x))))
  print(all.equal(1 * (x == x[10]), drop(sgdGMF::cpp.dirac(x, x[10]))))
  print(all.equal(1 * (x == x[77]), drop(sgdGMF::cpp.dirac(x, x[77]))))
}


## Test: step() ----
{
  x = seq(from = -3, to = +3, length = 201)
  print(all.equal(1 * (x < +1), drop(sgdGMF::cpp.step(x, +1))))
  print(all.equal(1 * (x < -1), drop(sgdGMF::cpp.step(x, -1, TRUE))))
  print(all.equal(1 * (x > -1), drop(sgdGMF::cpp.step(x, -1, FALSE))))
}

## Test: vech() ----
{
  a = matrix(1:9, nrow = 3, ncol = 3)
  v = c(1, 2, 5, 3, 6, 9)
  print(all.equal(v, drop(sgdGMF::cpp.vech(a))))
}

