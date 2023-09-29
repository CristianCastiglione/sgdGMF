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

## Test data ----
par(mfrow = c(1, 3))


## Test: gaussian ----
{
  n = 100
  x = seq(from = -3, to = +3, length = n)
  y = seq(from = -3, to = +3, length = n)
  z = rep(0, length = n)
  plot.link(x, sgdGMF::c_gaussian_variance(x), main = "Gaussian \n variance")
  plot.link(x, sgdGMF::c_gaussian_initialize(y), main = "Gaussian \n initialize")
  plot.link(x, sgdGMF::c_gaussian_devresid(z, x), main = "Gaussian \n devresid")
}

## Test: binomial ----
{
  n = 100
  x = seq(from = -0.001, to = +0.999, length = n)
  y = c(rep(0, length = n/2), rep(1, length = n/2))
  z = rep(0, length = n)
  plot.link(x, sgdGMF::c_binomial_variance(x), main = "Binomial \n variance")
  plot.link(x, sgdGMF::c_binomial_initialize(y), main = "Binomial \n initialize")
  plot.link(x, sgdGMF::c_binomial_devresid(z, x), main = "Binomial \n devresid")
}

## Test: poisson ----
{
  n = 100
  x = seq(from = 1, to = 10, length = n)
  y = seq(from = 1, to = 20, by = 1)
  z = rep(3, length = n)
  plot.link(x, sgdGMF::c_poisson_variance(x), main = "Poisson \n variance")
  plot.link(y, sgdGMF::c_poisson_initialize(y), main = "Poisson \n initialize")
  plot.link(x, sgdGMF::c_poisson_devresid(z, x), main = "Poisson \n devresid")
}

## Test: gamma ----
{
  x = seq(from = 0.1, to = 5, length = 100)
  y = seq(from = 0.1, to = 5, length = 100)
  z = rep(1, length = 100)
  plot.link(x, sgdGMF::c_gamma_variance(x), main = "Gamma \n variance")
  plot.link(x, sgdGMF::c_gamma_initialize(y), main = "Gamma \n initialize")
  plot.link(x, sgdGMF::c_gamma_devresid(z, x), main = "Gamma \n devresid")
}

