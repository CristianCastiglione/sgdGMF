# test-newton.R
# author: Cristian Castiglione
# creation: 02/10/2023
# last change: 04/10/2023

## Workspace setup ----
rm(list = ls())
graphics.off()

# Package compilation and import
devtools::load_all()

## Test: synthetic data ----
n = 100
m = 20
d = 3
p = 3
q = 4



family = poisson()

X = matrix(rnorm(n*p), nrow = n, ncol = p) / sqrt(3)
B = matrix(rnorm(m*p), nrow = m, ncol = p) / sqrt(3)
A = matrix(rnorm(n*q), nrow = n, ncol = q) / sqrt(3)
Z = matrix(rnorm(m*q), nrow = m, ncol = q) / sqrt(3)
U = matrix(rnorm(n*d), nrow = n, ncol = d) / sqrt(3)
V = matrix(rnorm(m*d), nrow = m, ncol = d) / sqrt(3)

eta = tcrossprod(cbind(X, A, U), cbind(B, Z, V))
mu = family$linkinv(eta)

Y = matrix(rpois(n*m, mu), nrow = n, ncol = m)

plot3D::image2D(log1p(Y))

logY = log(Y + 0.1)
B0 = t(solve(crossprod(X), crossprod(X, logY)))
A0 = t(solve(crossprod(Z), crossprod(Z, t(logY - tcrossprod(X, B0)))))
UV = svd::propack.svd(logY - tcrossprod(cbind(X, A0), cbind(B0, Z)), neig = d)
U0 = UV$u %*% diag(sqrt(UV$d))
V0 = UV$v %*% diag(sqrt(UV$d))

cfit = sgdGMF::c_fit_newton(
  Y = Y, X = X, B = B0, A = A0, Z = Z, U = U0, V = V0,
  familyname = "poisson", linkname = "log", ncomp = d,
  lambda = c(0, 0, 1, 0), maxiter = 500, stepsize = 0.1,
  tol = 1e-05, frequency = 50)

rfit = sgdGMF::sgdgmf(Y, X, Z, family = poisson(), ncomp = d,
                      init = list(niter = 0),
                      control = list(maxiter = 500, stepsize = 0.1))

fit$mu
fit$eta
fit$U
fit$V
tcrossprod(fit$U, fit$V)

plot(c(rfit$pred$mu), c(cfit$mu))
plot(c(Y), c(cfit$mu))
cor(c(rfit$pred$mu), c(cfit$mu))

plot3D::image2D(rfit$pred$mu)
plot3D::image2D(cfit$mu)
plot3D::image2D(Y)

all.equal(c(rfit$pred$mu), c(cfit$mu))



