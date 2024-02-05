# test-airwls.R
# author: Cristian Castiglione
# creation: 03/10/2023
# last change: 05/10/2023

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

get.glm.init = function (familyname = "gaussian", linkname = "identity") {
  f = NULL
  if (familyname == "gaussian") f = function (x) x
  if (familyname == "binomial") f = function (x) 2 * x - 1
  if (familyname == "poisson") f = function (x) log(x + 0.1)
  if (familyname == "gamma") f = function (x) log(x)
  return (f)
}

# Fit a vector GLM via sequential IRLS
vglm.fit = function (Y, X, familyname, linkname, penalty, transp) {
  family = get.glm.family(familyname, linkname)
  n = nrow(Y); m = ncol(Y); p = ncol(X)
  coef = NULL
  if (transp) {
    coef = matrix(NA, nrow = n, ncol = p)
    for (slice in 1:n) {
      y = Y[slice,]
      coef[slice,] = glm.fit(X, y, family = family)$coef
    }
  } else {
    coef = matrix(NA, nrow = m, ncol = p)
    for (slice in 1:m) {
      y = Y[,slice]
      coef[slice,] = glm.fit(X, y, family = family)$coef
    }
  }
  return (coef)
}


# Plotting function
plot.coef.path = function (
    y, X, coef, familyname = "gaussian", linkname = "identity",
    maxsteps = 100, offset = rep(0, nrow(X)), penalty = rep(0, ncol(X)),
    stepsize = 0.95, print = FALSE
) {
  result = matrix(NA, nrow = maxsteps+1, ncol = p)
  for (nsteps in 0:maxsteps) {
    result[nsteps+1,] = drop(sgdGMF::cpp.airwls.glmfit(
      beta = beta0, y = y, X = X,
      familyname = familyname, linkname = linkname,
      offset = offset, penalty = penalty, nsteps = nsteps,
      stepsize = stepsize, print = print))
  }

  ylim = range(c(c(result), coef))
  matplot(result, type = "o", pch = 20, ylim = ylim)
  abline(h = coef, lty = 2)
}

## Test: syntheric data ----
n = 200; p = 5
X = matrix(rnorm(n*p), nrow = n, ncol = p)
beta = rnorm(p, mean = 1)
eta = X %*% beta

## Test: glmfit gaussian-identity ----
{
  y = apply(eta / 3, MARGIN = 1, FUN = function(x) rnorm(1, mean = x, sd = 0.25))
  beta0 = rep(0, p)
  offset = rep(0, n)
  penalty = rep(0, p)
  familyname = "gaussian"
  linkname = "identity"

  glm.coef = glm.fit(X, y, family = get.glm.family(familyname, linkname))$coef
  gmf.coef = drop(sgdGMF::cpp.airwls.glmfit(
    beta = beta0, y = y, X = X, familyname = familyname, linkname = linkname,
    offset = offset, penalty = penalty, nsteps = 50, stepsize = 0.99, print = FALSE))

  print(all.equal(glm.coef, gmf.coef))
  plot.coef.path(y, X, glm.coef, familyname = familyname, linkname = linkname)
}

## Test: glmfit binomial-logit  ----
{
  y = apply(eta / 3, MARGIN = 1, FUN = function(x) rbinom(1, size = 1, prob = plogis(x)))
  beta0 = rep(0, p)
  offset = rep(0, n)
  penalty = rep(0, p)
  familyname = "binomial"
  linkname = "logit"

  glm.coef = glm.fit(X, y, family = get.glm.family(familyname, linkname))$coef
  gmf.coef = drop(sgdGMF::cpp.airwls.glmfit(
    beta = beta0, y = y, X = X, familyname = familyname, linkname = linkname,
    offset = offset, penalty = penalty, nsteps = 50, stepsize = 0.99, print = FALSE))

  print(all.equal(glm.coef, gmf.coef))
  plot.coef.path(y, X, glm.coef, familyname = familyname, linkname = linkname)
}


## Test: glmfit binomial-probit  ----
{
  y = apply(eta / 3, MARGIN = 1, FUN = function(x) rbinom(1, size = 1, prob = pnorm(x)))
  beta0 = rep(0, p)
  offset = rep(0, n)
  penalty = rep(0, p)
  familyname = "binomial"
  linkname = "probit"

  glm.coef = glm.fit(X, y, family = get.glm.family(familyname, linkname))$coef
  gmf.coef = drop(sgdGMF::cpp.airwls.glmfit(
    beta = beta0, y = y, X = X, familyname = familyname, linkname = linkname,
    offset = offset, penalty = penalty, nsteps = 100, stepsize = 0.99, print = FALSE))

  print(all.equal(glm.coef, gmf.coef))
  plot.coef.path(y, X, glm.coef, familyname = familyname, linkname = linkname)
}


## Test: glmfit binomial-cauchit  ----
{
  y = apply(eta / 3, MARGIN = 1, FUN = function(x) rbinom(1, size = 1, prob = pnorm(x)))
  beta0 = rep(0, p)
  offset = rep(0, n)
  penalty = rep(0, p)
  familyname = "binomial"
  linkname = "cauchit"

  glm.coef = glm.fit(X, y, family = get.glm.family(familyname, linkname))$coef
  gmf.coef = drop(sgdGMF::cpp.airwls.glmfit(
    beta = beta0, y = y, X = X, familyname = familyname, linkname = linkname,
    offset = offset, penalty = penalty, nsteps = 50, stepsize = 0.99, print = FALSE))

  print(all.equal(glm.coef, gmf.coef))
  plot.coef.path(y, X, glm.coef, familyname = familyname, linkname = linkname)
}


## Test: glmfit poisson-log ----
{
  y = apply(eta / 3, MARGIN = 1, FUN = function(x) rpois(1, lambda = exp(x)))
  beta0 = rep(0, p)
  offset = rep(0, n)
  penalty = rep(0, p)
  familyname = "poisson"
  linkname = "log"

  glm.coef = glm.fit(X, y, family = get.glm.family(familyname, linkname))$coef
  gmf.coef = drop(sgdGMF::cpp.airwls.glmfit(
    beta = beta0, y = y, X = X, familyname = familyname, linkname = linkname,
    offset = offset, penalty = penalty, nsteps = 50, stepsize = 0.99, print = FALSE))

  print(all.equal(glm.coef, gmf.coef))
  plot.coef.path(y, X, glm.coef, familyname = familyname, linkname = linkname)
}


## Test: glmfit gamma-log ----
{
  y = apply(eta / 3, MARGIN = 1, FUN = function(x) rgamma(1, shape = 1, rate = exp(x)))
  beta0 = rep(0, p)
  offset = rep(0, n)
  penalty = rep(0, p)
  familyname = "gamma"
  linkname = "log"

  glm.coef = glm.fit(X, y, family = get.glm.family(familyname, linkname))$coef
  gmf.coef = drop(sgdGMF::cpp.airwls.glmfit(
    beta = beta0, y = y, X = X, familyname = familyname, linkname = linkname,
    offset = offset, penalty = penalty, nsteps = 50, stepsize = 0.99, print = FALSE))

  print(all.equal(glm.coef, gmf.coef))
  plot.coef.path(y, X, glm.coef, familyname = familyname, linkname = linkname)
}

## Test: multivariate update ----
{
  n = 1000; m = 100; d = 5
  familyname = "poisson"
  linkname = "log"
  penalty = rep(0, length = d)
  offset = matrix(0, nrow = n, ncol = m)
  family = get.glm.family(familyname, linkname)
  U = matrix(rnorm(n*d), nrow = n, ncol = d)
  V = matrix(rnorm(m*d), nrow = m, ncol = d)
  eta = tcrossprod(U, V) / 4
  mu = family$linkinv(eta)
  if (familyname == "gaussian") Y = matrix(sapply(mu, FUN = function (x) rnorm(1, mean = x, sd = 0.25)), nrow = n, ncol = m)
  if (familyname == "binomial") Y = matrix(sapply(mu, FUN = function (x) rbinom(1, size = 1, prob = x)), nrow = n, ncol = m)
  if (familyname == "poisson") Y = matrix(sapply(mu, FUN = function (x) rpois(1, lambda = x)), nrow = n, ncol = m)
  if (familyname == "gamma") Y = matrix(sapply(mu, FUN = function (x) rgamma(1, shape = 1, rate = x)), nrow = n, ncol = m)

  t0 = proc.time()
  r.coef.v = vglm.fit(Y, U, familyname, linkname, penalty, FALSE)
  print((proc.time() - t0)[3])

  t0 = proc.time()
  r.coef.u = vglm.fit(Y, V, familyname, linkname, penalty, TRUE)
  print((proc.time() - t0)[3])

  t0 = proc.time()
  c.coef.v = sgdGMF::cpp.airwls.update(
    beta = V, Y = Y, X = U, familyname = familyname, linkname = linkname,
    idx = 1:d-1, offset = offset, penalty = penalty, transp = FALSE,
    nsteps = 100, stepsize = 0.99, print = FALSE, parallel = TRUE, nthreads = 4)
  print((proc.time() - t0)[3])

  t0 = proc.time()
  c.coef.u = sgdGMF::cpp.airwls.update(
    beta = U, Y = Y, X = V, familyname = familyname, linkname = linkname,
    idx = 1:d-1, offset = offset, penalty = penalty, transp = TRUE,
    nsteps = 100, stepsize = 0.99, print = FALSE, parallel = TRUE, nthreads = 4)
  print((proc.time() - t0)[3])


  print(all.equal(r.coef.v, c.coef.v))
  print(all.equal(r.coef.u, c.coef.u))

  par(mfrow = c(1,2))
  plot(r.coef.v, c.coef.v)
  plot(r.coef.u, c.coef.u)
  par(mfrow = c(1,1))
}

## Test: airwls fit ----
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

  r.gmffit = sgdGMF::sgdgmf(
    Y, X, Z, family = family, ncomp = d, method = "airwls", init = list(niter = 0),
    control = list(maxiter = 500, stepsize = 0.9, eps = 1e-08, tol = 1e-05,
                   damping = 1e-03, verbose = TRUE, frequency = 10))

  c.gmffit = sgdGMF::cpp.fit.airwls(
    Y, X, B0, A0, Z, U0, V0, familyname = familyname, linkname = linkname,
    ncomp = d, lambda = c(0,0,1,0), maxiter = 500, nsteps = 1, stepsize = 0.9,
    eps = 1e-08, nafill = 1, tol = 1e-05, damping = 1e-03, verbose = TRUE,
    frequency = 10, parallel = FALSE)

  c.gmffit = sgdGMF::cpp.fit.airwls(
    Y, X, B0, A0, Z, U0, V0, familyname = familyname, linkname = linkname,
    ncomp = d, lambda = c(0,0,1,0), maxiter = 500, nsteps = 1, stepsize = 0.9,
    eps = 1e-08, nafill = 1, tol = 1e-05, damping = 1e-03, verbose = TRUE,
    frequency = 10, parallel = TRUE)

  print(all.equal(r.gmffit$pred$mu, c.gmffit$mu))
  print(all.equal(r.gmffit$pred$eta, c.gmffit$eta))

  print(cor(c(r.gmffit$pred$mu), c(c.gmffit$mu)))
  print(cor(c(r.gmffit$pred$eta), c(c.gmffit$eta)))

  plot(r.gmffit$pred$mu, c.gmffit$mu)
  plot(r.gmffit$pred$eta, c.gmffit$eta)

  plot(mu, c.gmffit$mu)
  plot(eta, c.gmffit$eta)

}

## Test:: splatter data ----
{

  suppressPackageStartupMessages({
    library(splatter)
    library(scater)
  })

  n = 2000
  m = 250
  params = splatter::newSplatParams()
  params = splatter::setParam(params, "batchCells", n)
  params = splatter::setParam(params, "nGenes", m)
  params = splatter::setParam(params, "group.prob", c(0.1, 0.2, 0.2, 0.2, 0.3))
  params = splatter::setParam(params, "de.prob", c(0.3, 0.1, 0.2, 0.01, 0.1))
  params = splatter::setParam(params, "de.downProb", c(0.1, 0.4, 0.9, 0.6, 0.5))
  params = splatter::setParam(params, "de.facLoc", c(0.6, 0.1, 0.1, 0.01, 0.2))
  params = splatter::setParam(params, "de.facScale", c(0.1, 0.4, 0.2, 0.5, 0.4))
  # params = splatter::setParam(params, "seed", 140275)

  sim = splatter::splatSimulateGroups(params, verbose = FALSE)

  ## DIMENSION REDUCTION
  sim = scater::logNormCounts(sim)
  sim = scater::runPCA(sim)
  scater::plotPCA(sim, colour_by = "Group") +
    labs(title = "Different DE factors")

  ## DATA EXTRACTION
  counts = as.data.frame(counts(sim))
  cells = as.data.frame(colData(sim))
  genes = as.data.frame(rowData(sim))
  meta = metadata(sim)

  ## LATENT FACTOR ESTIMATION
  Y = matrix(NA, nrow = n, ncol = m)
  Y[] = t(counts)
  X = matrix(1, nrow = nrow(Y), ncol = 1)
  Z = matrix(1, nrow = ncol(Y), ncol = 1)

  ncomp = 5
  family = poisson()

  familyname = "poisson"
  linkname = "log"

  init = get.glm.init(familyname, linkname)
  R = init(Y)
  B0 = t(solve(crossprod(X), crossprod(X, R)))
  A0 = t(solve(crossprod(Z), crossprod(Z, t(R - tcrossprod(X, B0)))))
  UV = svd::propack.svd(R - tcrossprod(cbind(X, A0), cbind(B0, Z)), neig = ncomp)
  U0 = UV$u %*% diag(sqrt(UV$d))
  V0 = UV$v %*% diag(sqrt(UV$d))



  r.gmffit = sgdGMF::sgdgmf(
    Y, X, Z, family = family, ncomp = ncomp, method = "airwls", init = list(niter = 0),
    control = list(maxiter = 500, stepsize = 0.9, eps = 1e-08, tol = 1e-05,
                   damping = 1e-03, verbose = TRUE, frequency = 10))

  c.gmffit = sgdGMF::cpp.fit.airwls(
    Y, X, B0, A0, Z, U0, V0, familyname = familyname, linkname = linkname,
    ncomp = ncomp, lambda = c(0,0,1,0), maxiter = 500, nsteps = 1, stepsize = 0.9,
    eps = 1e-08, nafill = 1, tol = 1e-05, damping = 1e-03, verbose = TRUE,
    frequency = 10, parallel = FALSE)

  c.gmffit = sgdGMF::cpp.fit.airwls(
    Y, X, B0, A0, Z, U0, V0, familyname = familyname, linkname = linkname,
    ncomp = ncomp, lambda = c(0,0,1,0), maxiter = 500, nsteps = 1, stepsize = 0.9,
    eps = 1e-08, nafill = 1, tol = 1e-05, damping = 1e-03, verbose = TRUE,
    frequency = 10, parallel = TRUE)
}



## End of file ----
