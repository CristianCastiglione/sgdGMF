library(sgdGMF)

# Set the data dimensions
n = 100; m = 20; d = 5

# Generate data using Poisson, Binomial and Gamma models
data_pois = sim.gmf.data(n = n, m = m, ncomp = d, family = poisson())
data_bin = sim.gmf.data(n = n, m = m, ncomp = d, family = binomial())
data_gam = sim.gmf.data(n = n, m = m, ncomp = d, family = Gamma(link = "log"), dispersion = 0.25)

# Estimate the GMF parameters assuming 3 latent factors
gmf_pois = sgdgmf.fit(data_pois$Y, ncomp = 3, family = poisson(), method = "airwls")
gmf_bin = sgdgmf.fit(data_bin$Y, ncomp = 3, family = binomial(), method = "airwls")
gmf_gam = sgdgmf.fit(data_gam$Y, ncomp = 3, family = Gamma(link = "log"), method = "airwls")

# Get the fitted values in the link and response scales
mu_hat_pois = fitted(gmf_pois, type = "response")
mu_hat_bin = fitted(gmf_bin, type = "response")
mu_hat_gam = fitted(gmf_gam, type = "response")

# Compare the results
par(mfrow = c(3,3), mar = c(1,1,3,1))
image(data_pois$Y, axes = FALSE, main = expression(Y[Pois]))
image(data_pois$mu, axes = FALSE, main = expression(mu[Pois]))
image(mu_hat_pois, axes = FALSE, main = expression(hat(mu)[Pois]))
image(data_bin$Y, axes = FALSE, main = expression(Y[Bin]))
image(data_bin$mu, axes = FALSE, main = expression(mu[Bin]))
image(mu_hat_bin, axes = FALSE, main = expression(hat(mu)[Bin]))
image(data_gam$Y, axes = FALSE, main = expression(Y[Gam]))
image(data_gam$mu, axes = FALSE, main = expression(mu[Gam]))
image(mu_hat_gam, axes = FALSE, main = expression(hat(mu)[Gam]))

