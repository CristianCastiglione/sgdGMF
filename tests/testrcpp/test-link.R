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
  plot.link(x, sgdGMF::cpp.link.identity.linkfun(x), main = "linkfun")
  plot.link(x, sgdGMF::cpp.link.identity.linkinv(x), main = "linkinv")
  plot.link(x, sgdGMF::cpp.link.identity.mueta(x), main = "linkmueta")
}

## Test: logit ----
{
  plot.link(x, sgdGMF::cpp.link.logit.linkfun(z), main = "linkfun")
  plot.link(x, sgdGMF::cpp.link.logit.linkinv(x), main = "linkinv")
  plot.link(x, sgdGMF::cpp.link.logit.mueta(x), main = "linkmueta")
}

## Test: probit ----
{
  plot.link(x, sgdGMF::cpp.link.probit.linkfun(z), main = "linkfun")
  plot.link(x, sgdGMF::cpp.link.probit.linkinv(x), main = "linkinv")
  plot.link(x, sgdGMF::cpp.link.probit.mueta(x), main = "linkmueta")
}

## Test: cauchy ----
{
  plot.link(x, sgdGMF::cpp.link.cauchy.linkfun(z), main = "linkfun")
  plot.link(x, sgdGMF::cpp.link.cauchy.linkinv(x), main = "linkinv")
  plot.link(x, sgdGMF::cpp.link.cauchy.mueta(x), main = "linkmueta")
}

## Test: cloglog ----
{
  plot.link(x, sgdGMF::cpp.link.cloglog.linkfun(z), main = "linkfun")
  plot.link(x, sgdGMF::cpp.link.cloglog.linkinv(x), main = "linkinv")
  plot.link(x, sgdGMF::cpp.link.cloglog.mueta(x), main = "linkmueta")
}

## Test: log ----
{
  plot.link(x, sgdGMF::cpp.link.log.linkfun(y), main = "linkfun")
  plot.link(x, sgdGMF::cpp.link.log.linkinv(x), main = "linkinv")
  plot.link(x, sgdGMF::cpp.link.log.mueta(x), main = "linkmueta")
}

## Test: inverse ----
{
  plot.link(x, sgdGMF::cpp.link.inverse.linkfun(z), main = "linkfun")
  plot.link(x, sgdGMF::cpp.link.inverse.linkinv(z), main = "linkinv")
  plot.link(x, sgdGMF::cpp.link.inverse.mueta(z), main = "linkmueta")
}

## Test: sqrt ----
{
  plot.link(x, sgdGMF::cpp.link.sqrt.linkfun(z), main = "linkfun")
  plot.link(x, sgdGMF::cpp.link.sqrt.linkinv(z), main = "linkinv")
  plot.link(x, sgdGMF::cpp.link.sqrt.mueta(z), main = "linkmueta")
}

