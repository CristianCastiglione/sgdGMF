# test-link.R
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
x = seq(from = -3, to = +3, length = 201)
y = seq(from = 0.1, to = +5, length = 201)
z = seq(from = 0.001, to = 0.999, length = 201)

par(mfrow = c(1, 3))

## Test: identity ----
{
  plot.link(x, sgdGMF::c_link_identity_linkfun(x), main = "linkfun")
  plot.link(x, sgdGMF::c_link_identity_linkinv(x), main = "linkinv")
  plot.link(x, sgdGMF::c_link_identity_mueta(x), main = "linkmueta")
}

## Test: logit ----
{
  plot.link(x, sgdGMF::c_link_logit_linkfun(z), main = "linkfun")
  plot.link(x, sgdGMF::c_link_logit_linkinv(x), main = "linkinv")
  plot.link(x, sgdGMF::c_link_logit_mueta(x), main = "linkmueta")
}

## Test: probit ----
{
  plot.link(x, sgdGMF::c_link_probit_linkfun(z), main = "linkfun")
  plot.link(x, sgdGMF::c_link_probit_linkinv(x), main = "linkinv")
  plot.link(x, sgdGMF::c_link_probit_mueta(x), main = "linkmueta")
}

## Test: cauchy ----
{
  plot.link(x, sgdGMF::c_link_cauchy_linkfun(z), main = "linkfun")
  plot.link(x, sgdGMF::c_link_cauchy_linkinv(x), main = "linkinv")
  plot.link(x, sgdGMF::c_link_cauchy_mueta(x), main = "linkmueta")
}

## Test: cloglog ----
{
  plot.link(x, sgdGMF::c_link_cloglog_linkfun(z), main = "linkfun")
  plot.link(x, sgdGMF::c_link_cloglog_linkinv(x), main = "linkinv")
  plot.link(x, sgdGMF::c_link_cloglog_mueta(x), main = "linkmueta")
}

## Test: log ----
{
  plot.link(x, sgdGMF::c_link_log_linkfun(y), main = "linkfun")
  plot.link(x, sgdGMF::c_link_log_linkinv(x), main = "linkinv")
  plot.link(x, sgdGMF::c_link_log_mueta(x), main = "linkmueta")
}

## Test: inverse ----
{
  plot.link(x, sgdGMF::c_link_inverse_linkfun(z), main = "linkfun")
  plot.link(x, sgdGMF::c_link_inverse_linkinv(z), main = "linkinv")
  plot.link(x, sgdGMF::c_link_inverse_mueta(z), main = "linkmueta")
}

## Test: sqrt ----
{
  plot.link(x, sgdGMF::c_link_sqrt_linkfun(z), main = "linkfun")
  plot.link(x, sgdGMF::c_link_sqrt_linkinv(z), main = "linkinv")
  plot.link(x, sgdGMF::c_link_sqrt_mueta(z), main = "linkmueta")
}
