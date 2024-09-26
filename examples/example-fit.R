library(sgdGMF)

# Set the data dimensions
n = 100; m = 20; d = 5

# Set the time range, phases, frequences and amplitudes of the underlying signals
time = seq(from = 0, to = 1, length = n)
phase = sort(2*runif(d), decreasing = FALSE)
freq = sort(2*runif(d), decreasing = FALSE)
amp = sort(runif(d), decreasing = TRUE)

Time = tcrossprod(time, rep(1, d))
Phase = tcrossprod(rep(1, n), phase)
Freq = tcrossprod(rep(1, n), freq)
Amp = diag(amp)

# Combine the latent signals using independent Gaussian coefficients
U = sin(2 * pi * Freq * (Time + Phase)) %*% Amp
V = matrix(rnorm(m*d), nrow = m, ncol = d)
eta = tcrossprod(U, V)

# Generate data using Poisson, Binomial and Gamma models
mu_pois = exp(eta)
mu_bin = plogis(2*eta)
mu_gam = exp(eta)

Y_pois = matrix(rpois(n*m, lambda = mu_pois), nrow = n, ncol = m)
Y_bin = matrix(rbinom(n*m, size = 1, prob = mu_bin), nrow = n, ncol = m)
Y_gam = matrix(rgamma(n*m, shape = 2, scale = mu_gam), nrow = n, ncol = m)

# Initialize the GMF parameters assuming 3 latent factors
gmf_pois = sgdgmf.fit(Y_pois, ncomp = 3, family = poisson(), method = "airwls")
gmf_bin = sgdgmf.fit(Y_bin, ncomp = 3, family = binomial(), method = "airwls")
gmf_gam = sgdgmf.fit(Y_gam, ncomp = 3, family = Gamma(link = "log"), method = "airwls")

# Get the fitted values in the link and response scales
mu_hat_pois = fitted(gmf_pois, type = "response")
mu_hat_bin = fitted(gmf_bin, type = "response")
mu_hat_gam = fitted(gmf_gam, type = "response")

# Compare the results
par(mfrow = c(3,3), mar = c(1,1,3,1))
image(Y_pois, axes = FALSE, main = expression(Y[Pois]))
image(mu_pois, axes = FALSE, main = expression(mu[Pois]))
image(mu_hat_pois, axes = FALSE, main = expression(hat(mu)[Pois]))
image(Y_bin, axes = FALSE, main = expression(Y[Bin]))
image(mu_bin, axes = FALSE, main = expression(mu[Bin]))
image(mu_hat_bin, axes = FALSE, main = expression(hat(mu)[Bin]))
image(Y_gam, axes = FALSE, main = expression(Y[Gam]))
image(mu_gam, axes = FALSE, main = expression(mu[Gam]))
image(mu_hat_gam, axes = FALSE, main = expression(hat(mu)[Gam]))

