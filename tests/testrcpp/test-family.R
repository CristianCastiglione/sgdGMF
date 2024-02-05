# test-family.R
# author: Cristian Castiglione
# creation: 29/09/2023
# last change: 29/09/2023

## Workspace setup ----
rm(list = ls())
graphics.off()

# Package compilation and import
devtools::load_all()

plot.link <- function (x, y, main = "") {
  plot(x, y, type = "l", xlab = "x", ylab = "link", main = main)
}

par(mfrow = c(1, 3))

## Test: gaussian ----
{
  n = 100
  x = seq(from = -3, to = +3, length = n)
  y = seq(from = -3, to = +3, length = n)
  z = rep(0, length = n)
  plot.link(x, sgdGMF::cpp.gaussian.variance(x), main = "Gaussian \n variance")
  plot.link(x, sgdGMF::cpp.gaussian.initialize(y), main = "Gaussian \n initialize")
  plot.link(x, sgdGMF::cpp.gaussian.devresid(z, x), main = "Gaussian \n devresid")

  r.variance = gaussian()$variance(x)
  c.variance = drop(sgdGMF::cpp.gaussian.variance(x))
  print(all.equal(r.variance, c.variance))

  r.devresid = gaussian()$dev.resid(z, x, 1)
  c.devresid = drop(sgdGMF::cpp.gaussian.devresid(z, x))
  print(all.equal(r.devresid, c.devresid))
}

## Test: binomial ----
{
  n = 100
  x = seq(from = +0.001, to = +0.999, length = n)
  y = c(rep(0, length = n/2), rep(1, length = n/2))
  z = rep(0, length = n)
  plot.link(x, sgdGMF::cpp.binomial.variance(x), main = "Binomial \n variance")
  plot.link(x, sgdGMF::cpp.binomial.initialize(y), main = "Binomial \n initialize")
  plot.link(x, sgdGMF::cpp.binomial.devresid(z, x), main = "Binomial \n devresid")

  r.variance = binomial()$variance(x)
  c.variance = drop(sgdGMF::cpp.binomial.variance(x))
  print(all.equal(r.variance, c.variance))

  r.devresid = binomial()$dev.resid(z, x, 1)
  c.devresid = drop(sgdGMF::cpp.binomial.devresid(z, x))
  print(all.equal(r.devresid, c.devresid))
}

## Test: poisson ----
{
  n = 100
  x = seq(from = 1, to = 10, length = n)
  y = seq(from = 1, to = 20, by = 1)
  z = rep(3, length = n)
  plot.link(x, sgdGMF::cpp.poisson.variance(x), main = "Poisson \n variance")
  plot.link(y, sgdGMF::cpp.poisson.initialize(y), main = "Poisson \n initialize")
  plot.link(x, sgdGMF::cpp.poisson.devresid(z, x), main = "Poisson \n devresid")

  r.variance = poisson()$variance(x)
  c.variance = drop(sgdGMF::cpp.poisson.variance(x))
  print(all.equal(r.variance, c.variance))

  r.devresid = poisson()$dev.resid(z, x, 1)
  c.devresid = drop(sgdGMF::cpp.poisson.devresid(z, x))
  print(all.equal(r.devresid, c.devresid))
}

## Test: gamma ----
{
  n = 100
  x = seq(from = 0.1, to = 5, length = n)
  y = seq(from = 0.1, to = 5, length = n)
  z = rep(1, length = n)
  plot.link(x, sgdGMF::cpp.gamma.variance(x), main = "Gamma \n variance")
  plot.link(x, sgdGMF::cpp.gamma.initialize(y), main = "Gamma \n initialize")
  plot.link(x, sgdGMF::cpp.gamma.devresid(z, x), main = "Gamma \n devresid")

  r.variance = Gamma()$variance(x)
  c.variance = drop(sgdGMF::cpp.gamma.variance(x))
  print(all.equal(r.variance, c.variance))

  r.devresid = Gamma()$dev.resid(z, x, 1)
  c.devresid = drop(sgdGMF::cpp.gamma.devresid(z, x))
  print(all.equal(r.devresid, c.devresid))
}

## Test: negative binomial ----
{
  n = 100
  x = seq(from = 0.1, to = 5, length = n)
  y = seq(from = 0.1, to = 5, length = n)
  z = rep(1, length = n)
  plot.link(x, sgdGMF::cpp.negbinom.variance(x), main = "Negative Binomial \n variance")
  plot.link(x, sgdGMF::cpp.negbinom.initialize(y), main = "Negative Binomial \n initialize")
  plot.link(x, sgdGMF::cpp.negbinom.devresid(z, x), main = "Negative Binomial \n devresid")

  r.variance = MASS::neg.bin(10)$variance(x)
  c.variance = drop(sgdGMF::cpp.negbinom.variance(x))
  print(all.equal(r.variance, c.variance))

  r.devresid = MASS::neg.bin(10)$dev.resid(z, x, 1)
  c.devresid = drop(sgdGMF::cpp.negbinom.devresid(z, x))
  print(all.equal(r.devresid, c.devresid))
}


