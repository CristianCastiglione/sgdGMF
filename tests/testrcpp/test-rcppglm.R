# test-rcppglm.R
# author: Cristian Castiglione
# creation: 01/10/2023
# last change: 01/10/2023

## Workspace setup ----
rm(list = ls())
graphics.off()

# Package compilation and import
devtools::load_all()

## Test: link-family constructors ----
sgdGMF::test_make_gaussian("identity")
sgdGMF::test_make_binomial("logit")
sgdGMF::test_make_binomial("probit")
sgdGMF::test_make_binomial("cauchy")
sgdGMF::test_make_binomial("cloglog")
sgdGMF::test_make_poisson("log")
sgdGMF::test_make_poisson("sqrt")
sgdGMF::test_make_gamma("log")
sgdGMF::test_make_gamma("inverse")
sgdGMF::test_make_gamma("sqrt")

## End of file ----
