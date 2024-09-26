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
ncomp_pois = sgdgmf.rank(Y_pois, family = poisson(), normalize = TRUE)
ncomp_bin = sgdgmf.rank(Y_bin, family = binomial(), normalize = TRUE)
ncomp_gam = sgdgmf.rank(Y_gam, family = Gamma(link = "log"), normalize = TRUE)

# Get the selected number of components
print(paste("Poisson:", ncomp_pois$ncomp))
print(paste("Binomial:", ncomp_bin$ncomp))
print(paste("Gamma:", ncomp_gam$ncomp))

# Plot the screeplot used for the component determination
par(mfrow = c(3,1))
barplot(ncomp_pois$lambdas, main = "Poisson screeplot")
barplot(ncomp_bin$lambdas, main = "Binomial screeplot")
barplot(ncomp_gam$lambdas, main = "Gamma screeplot")
par(mfrow = c(1,1))
