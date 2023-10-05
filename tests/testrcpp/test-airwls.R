# test-airwls.R
# author: Cristian Castiglione
# creation: 03/10/2023
# last change: 04/10/2023

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
    result[nsteps+1,] = drop(sgdGMF::c_airwls_glmfit(
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
  gmf.coef = drop(sgdGMF::c_airwls_glmfit(
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
  gmf.coef = drop(sgdGMF::c_airwls_glmfit(
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
  gmf.coef = drop(sgdGMF::c_airwls_glmfit(
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
  gmf.coef = drop(sgdGMF::c_airwls_glmfit(
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
  gmf.coef = drop(sgdGMF::c_airwls_glmfit(
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
  gmf.coef = drop(sgdGMF::c_airwls_glmfit(
    beta = beta0, y = y, X = X, familyname = familyname, linkname = linkname,
    offset = offset, penalty = penalty, nsteps = 50, stepsize = 0.99, print = FALSE))

  print(all.equal(glm.coef, gmf.coef))
  plot.coef.path(y, X, glm.coef, familyname = familyname, linkname = linkname)
}

## Test: multivariate gaussian-identity ----
{
  familyname = "poisson"
  linkname = "log"
  penalty = rep(0, length = d)
  offset = matrix(0, nrow = n, ncol = m)
  family = get.glm.family(familyname, linkname)
  n = 100; m = 10; d = 3
  U = matrix(rnorm(n*d), nrow = n, ncol = d)
  V = matrix(rnorm(m*d), nrow = m, ncol = d)
  eta = tcrossprod(U, V) / 4
  mu = family$linkinv(eta)
  if (familyname == "gaussian") Y = matrix(sapply(mu, FUN = function (x) rnorm(1, mean = x, sd = 0.25)), nrow = n, ncol = m)
  if (familyname == "binomial") Y = matrix(sapply(mu, FUN = function (x) rbinom(1, size = 1, prob = x)), nrow = n, ncol = m)
  if (familyname == "poisson") Y = matrix(sapply(mu, FUN = function (x) rpois(1, lambda = x)), nrow = n, ncol = m)
  if (familyname == "gamma") Y = matrix(sapply(mu, FUN = function (x) rgamma(1, shape = 1, rate = x)), nrow = n, ncol = m)

  r.coef.v = vglm.fit(Y, U, familyname, linkname, FALSE)
  r.coef.u = vglm.fit(Y, V, familyname, linkname, TRUE)

  c.coef.v = sgdGMF::c_airwls_update(
    beta = V, Y = Y, X = U, familyname = familyname, linkname = linkname,
    idx = 1:d-1, offset = offset, penalty = penalty, transp = FALSE,
    nsteps = 100, stepsize = 0.99, print = FALSE, parallel = TRUE)

  c.coef.u = sgdGMF::c_airwls_update(
    beta = U, Y = Y, X = V, familyname = familyname, linkname = linkname,
    idx = 1:d-1, offset = offset, penalty = penalty, transp = TRUE,
    nsteps = 100, stepsize = 0.99, print = FALSE, parallel = TRUE)

  print(all.equal(r.coef.v, c.coef.v))
  print(all.equal(r.coef.u, c.coef.u))

  par(mfrow = c(1,2))
  plot(r.coef.v, c.coef.v)
  plot(r.coef.u, c.coef.u)
  par(mfrow = c(1,1))
}







## End of file ----
