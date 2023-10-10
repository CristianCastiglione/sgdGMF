# test-sgd.R
# author: Cristian Castiglione
# creation: 08/10/2023
# last change: 08/10/2023

## Workspace setup ----
rm(list = ls())
graphics.off()

# Package compilation and import
devtools::load_all()

# Family constructor
get.glm.family = function (familyname = "gaussian", linkname = "identity") {
  f = NULL
  if (familyname == "gaussian") f = gaussian(link = linkname)
  if (familyname == "binomial") f = binomial(link = linkname)
  if (familyname == "poisson") f = poisson(link = linkname)
  if (familyname == "gamma") f = Gamma(link = linkname)
  return (f)
}

# Initialization
get.glm.init = function (familyname = "gaussian", linkname = "identity") {
  f = NULL
  if (familyname == "gaussian") f = function (x) x
  if (familyname == "binomial") f = function (x) 2 * x - 1
  if (familyname == "poisson") f = function (x) log(x + 0.1)
  if (familyname == "gamma") f = function (x) log(x)
  return (f)
}

## Test: sgd fit ----

{
  n = 1000; m = 100; p = 3; q = 1; d = 5
  familyname = "poisson"
  linkname = "log"
  penalty = rep(0, length = d)
  offset = matrix(0, nrow = n, ncol = m)
  family = get.glm.family(familyname, linkname)
  init = get.glm.init(familyname, linkname)
  X = matrix(rnorm(n*p, sd = 0.9), nrow = n, ncol = p) / sqrt(p)
  B = matrix(rnorm(m*p, sd = 0.9), nrow = m, ncol = p) / sqrt(p)
  A = matrix(rnorm(n*q, sd = 0.9), nrow = n, ncol = q) / sqrt(q)
  Z = matrix(rnorm(m*q, sd = 0.9), nrow = m, ncol = q) / sqrt(q)
  U = matrix(rnorm(n*d, sd = 0.9), nrow = n, ncol = d) / sqrt(d)
  V = matrix(rnorm(m*d, sd = 0.9), nrow = m, ncol = d) / sqrt(d)
  eta = tcrossprod(cbind(X, A, U), cbind(B, Z, V)) / 4
  mu = family$linkinv(eta)
  if (familyname == "gaussian") Y = matrix(sapply(mu, FUN = function (x) rnorm(1, mean = x, sd = 0.25)), nrow = n, ncol = m)
  if (familyname == "binomial") Y = matrix(sapply(mu, FUN = function (x) rbinom(1, size = 1, prob = x)), nrow = n, ncol = m)
  if (familyname == "poisson") Y = matrix(sapply(mu, FUN = function (x) rpois(1, lambda = x)), nrow = n, ncol = m)
  if (familyname == "gamma") Y = matrix(sapply(mu, FUN = function (x) rgamma(1, shape = 1, rate = x)), nrow = n, ncol = m)

  R = init(Y)
  B0 = t(solve(crossprod(X), crossprod(X, R)))
  A0 = t(solve(crossprod(Z), crossprod(Z, t(R - tcrossprod(X, B0)))))
  UV = svd::propack.svd(R - tcrossprod(cbind(X, A0), cbind(B0, Z)), neig = d)
  U0 = UV$u %*% diag(sqrt(UV$d))
  V0 = UV$v %*% diag(sqrt(UV$d))

  r.sgdfit = sgdGMF::sgdgmf(
    Y, X, Z, family = family, ncomp = d, method = "airwls", init = list(niter = 0),
    control = list(maxiter = 200, stepsize = 0.9, eps = 1e-08, tol = 1e-05,
                   damping = 1e-03, frequency = 10))

  c.gmffit = sgdGMF::c_fit_airwls(
    Y, X, B0, A0, Z, U0, V0, familyname = familyname, linkname = linkname,
    ncomp = d, lambda = c(0,0,1,0), maxiter = 200, nsteps = 1, stepsize = 0.9,
    eps = 1e-08, nafill = 1, tol = 1e-05, damping = 1e-03, verbose = TRUE,
    frequency = 10, parallel = FALSE)

  c.wlsfit = sgdGMF::c_fit_airwls(
    Y, X, B0, A0, Z, U0, V0, familyname = familyname, linkname = linkname,
    ncomp = d, lambda = c(0,0,1,0), maxiter = 500, nsteps = 1, stepsize = 0.9,
    eps = 1e-08, nafill = 1, tol = 1e-05, damping = 1e-03, verbose = TRUE,
    frequency = 10, parallel = TRUE)

  c.ntnfit = sgdGMF::c_fit_newton(
    Y, X, B0, A0, Z, U0, V0, familyname = familyname, linkname = linkname,
    ncomp = d, lambda = c(0,0,1,0), maxiter = 100, stepsize = 0.5,
    eps = 1e-08, nafill = 1, tol = 1e-05, damping = 1e-03, verbose = TRUE,
    frequency = 10, parallel = TRUE)

  r.sgdfit = sgdGMF::sgdgmf(
    Y, X, Z, family = family, ncomp = d, method = "b-sgd", init = list(niter = 0),
    control = list(maxiter = 5000, rate0 = 0.01, size = c(100, 10), frequency = 500))

  c.sgdfit = sgdGMF::c_fit_bsgd(
    Y, X, B0, A0, Z, U0, V0, familyname = familyname, linkname = linkname,
    ncomp = d, lambda = c(0,0,1,0), maxiter = 5000, burn = 1, rate0 = 0.01,
    size1 = 100, size2 = 10, frequency = 500)

  c.sgdfit = sgdGMF::c_fit2_bsgd(
    Y, X, B0, A0, Z, U0, V0, familyname = familyname, linkname = linkname,
    ncomp = d, lambda = c(0,0,1,0), maxiter = 500, burn = 1, rate0 = 0.5, decay = 1.0,
    size1 = 100, size2 = 10, frequency = 50)

  c.sgdfit = sgdGMF::c_fit_csgd(
    Y, X, B0, A0, Z, U0, V0, familyname = familyname, linkname = linkname,
    ncomp = d, lambda = c(0,0,1,0), maxiter = 1000, burn = 1, rate0 = 0.01,
    size1 = 10, size2 = 10, frequency = 100)

  c.sgdfit = sgdGMF::c_fit_msgd(
    Y, X, B0, A0, Z, U0, V0, familyname = familyname, linkname = linkname,
    ncomp = d, lambda = c(0,0,1,0), maxiter = 100, burn = 1, rate0 = 0.1,
    size = 100, frequency = 10)

  print(all.equal(c.wlsfit$mu, c.sgdfit$mu))
  print(all.equal(c.wlsfit$eta, c.sgdfit$eta))

  print(cor(c(c.wlsfit$mu), c(c.sgdfit$mu)))
  print(cor(c(c.wlsfit$eta), c(c.sgdfit$eta)))

  plot(c.wlsfit$mu, c.sgdfit$mu)
  plot(c.wlsfit$eta, c.sgdfit$eta)

  plot(mu, c.wlsfit$mu)
  plot(eta, c.wlsfit$eta)

  plot(Y, c.sgdfit$mu)
  plot(mu, c.sgdfit$mu)
  plot(eta, c.sgdfit$eta)

}


