
## WORKSPACE SETUP ----

# Clean the workspace
# rm(list = ls())
# graphics.off()

# Load the GMF functions
# devtools::load_all()

# Import the needed packages
suppressPackageStartupMessages({
  # Latent variable modelling
  require(MASS)
  library(gllvm)
  library(glmpca)
  library(NewWave)

  # Matrix factorization
  library(NMF)
  library(NNLM)
  library(cmfrec)
  library(gmf)

  # Low-dimensional embedding
  library(Rtsne)
  library(umap)

  # Omics classes and methods
  library(scry)
  library(splatter)
  library(scater)
  library(zinbwave)

  # Clustering
  library(cluster)

  # Visualization
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(GGally)
  library(factoextra)
})

## GLOBAL VARIABLES ----
NCORES = 8

## UTILITIES ----

## Join several paths
join.path <- function (...) {
  paste(..., sep = "/")
}

## Join several strings
join.string <- function (...) {
  paste(..., sep = "")
}

## TRANSFORMATIONS ----

minmax <- function (x, a = -1, b = +1) {
  minx <- min(x); maxx <- max(x)
  (b - a) * (x - minx) / (maxx - minx) + a
}

boxcox <- function (x, lambda) {
  (x^lambda - 1) / lambda
}

log1p <- function (x) {
  log(1 + x)
}

log10p1 <- function (x) {
  log10(1 + x)
}

slog10p1 <- function (x) {
  sign(x) * log10(1 + abs(x))
}

ssqrt <- function (x) {
  sign(x) * sqrt(abs(x) + 1e-03)
}

rdirichlet = function (n, m = 2, alpha = 1) {
  x = rgamma(n*m, shape = alpha, rate = 1)
  x = matrix(x, nrow = m, ncol = n)
  y = t(apply(x, 2, function(x) x / sum(x)))
  return (y)
}

## ERROR MEASURES ----

# Residual sum of squares
rss <- function (y, x, f = log1p) {
  m = mean(y, na.rm = TRUE)
  err0 = mean(c(f(y) - f(m))^2, na.rm = TRUE)
  errf = mean(c(f(y) - f(x))^2, na.rm = TRUE)
  return (errf / err0)
}

# Root mean squared error
rmse <- function (y, x, f = log1p) {
  sqrt(mean((f(y) - f(x))^2, na.rm = TRUE))
}

# Cosine distance
cosdist <- function (y, x, f = log1p) {
  xy <- sum(f(y) * f(x), na.rm = TRUE)
  xx <- sqrt(sum(f(y)^2, na.rm = TRUE))
  yy <- sqrt(sum(f(x)^2, na.rm = TRUE))
  return(1 - abs(xy) / (xx * yy))
}

# Explained deviance
expdev <- function (y, x, family) {
  m <- mean(y, na.rm = TRUE)
  # nulldev <- mean(family$dev.resids(y, m, 1), na.rm = TRUE)
  # fitdev <- mean(family$dev.resids(y, x, 1), na.rm = TRUE)
  dev0 <- matrix.deviance(m, y, family = family)
  devf <- matrix.deviance(x, y, family = family)
  return (devf / dev0)
}

avgsil <- function (y, g) {
  mean(cluster::silhouette(as.numeric(g), dist(y))[,3])
}

# Summary error matrix
error.matrix = function (y, ...) {
  object.list = list(...)
  error = data.frame(Model = c(), Time = c(), RSS = c(), Cos = c(), Dev = c())
  for (k in 1:length(object.list)) {
    .object = object.list[[k]]
    .mu = .object$mu
    .model = .object$model
    .time = round(.object$time[3], 4)
    .rss = round(rss(y, .mu), 4)
    .cos = round(cosdist(y, .mu), 4)
    .dev = round(expdev(y, .mu, family = poisson()), 4)
    .error = data.frame(Model = .model, Time = .time, RSS = .rss, Cos = .cos, Dev = .dev)
    error = rbind(error, .error)
    rownames(error)[k] = k
  }
  return (error)
}

# Summary error matrix + silhouette
error.summary = function (y, g, ...) {
  object.list = list(...)
  error = data.frame(Model = c(), Time = c(), RSS = c(), Cos = c(), Dev = c(), Sil = c())
  for (k in 1:length(object.list)) {
    .object = object.list[[k]]
    .mu = .object$mu
    .model = .object$model
    .time = ifelse(is.null(.object$time), NA, round(.object$time[3], 4))
    .rss = ifelse(is.null(.mu), NA, round(rss(y, .mu), 4))
    .cos = ifelse(is.null(.mu), NA, round(cosdist(y, .mu), 4))
    .dev = ifelse(is.null(.mu), NA, round(expdev(y, .mu, family = poisson()), 4))
    .sil = ifelse(is.null(.object$tsne), NA, round(avgsil(.object$tsne, g), 4))
    .error = data.frame(Model = .model, Time = .time, RSS = .rss, Cos = .cos, Dev = .dev, Sil = .sil)
    error = rbind(error, .error)
    rownames(error)[k] = k
  }
  return (error)
}

## POSTPROCESSING AND PLOTTING ----

## Plotting function
plot.coeff.matrix <- function (scores, limits = NULL, colours = NULL,
                               transpose = TRUE, symmetric = TRUE) {

  if (is.null(colours)) colours = c("navyblue", "grey95", "darkred")

  if (transpose) scores <- t(scores)

  scores[scores > quantile(scores, .99, na.rm = TRUE)] <- quantile(scores, .99, na.rm = TRUE)
  scores[scores < quantile(scores, .01, na.rm = TRUE)] <- quantile(scores, .01, na.rm = TRUE)

  df <- expand.grid(x = 1:nrow(scores), y = 1:ncol(scores))
  df$z <- as.vector(scores)

  if (is.null(limits)) {
    if (symmetric) {
      limits <- c(-1,+1) * max(abs(scores))
    } else {
      limits <- range(scores)
    }
  }

  plt <- ggplot(data = df, mapping = aes(x = x, y = y, fill = z)) +
    geom_raster() + labs(x = "", y = "", fill = "") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(colours = colours, limits = limits) +
    theme_bw() + theme(axis.text = element_blank(), axis.ticks = element_blank())
  plt
}

# Compare the latent representations
plot.eigvectors = function (u, idx = 1:2, title = "Biplot") {
  ggplot(data = data.frame(x = u[,idx[1]], y = u[,idx[2]], z = 1:nrow(u))) +
    geom_hline(yintercept = 0, lty = 2, color = 1) +
    geom_vline(xintercept = 0, lty = 2, color = 1) +
    geom_point(mapping = aes(x = x, y = y, color = z), size = 2) +
    geom_text(mapping = aes(x = x, y = y - 0.025, color = z, label = z), size = 3) +
    scale_color_gradientn(colors = c("navyblue", "gold"), name = "") +
    labs(x = "PC1", y = "PC2", title = title) + theme_bw() +
    theme(axis.title = element_blank())
}

## TRAIN-TEST SPLIT ----

## Train-test split
train.test.split <- function (y, test = 0.3) {

  n = nrow(y)
  m = ncol(y)

  mask <- cbind(x = sample.int(n = n, size = floor(test * n * m), replace = TRUE),
                y = sample.int(n = m, size = floor(test * n * m), replace = TRUE))

  ycc <- y
  yna <- y
  ynn <- y

  ynn[mask] <- NA
  isna <- is.na(ynn)

  yna[!isna] <- NA

  list(mask = mask, isna = isna, ycc = ycc, yna = yna, ynn = ynn)
}

train.test.split <- function (y, test = 0.3) {

  n = nrow(y)
  m = ncol(y)

  idx = cbind(x = sample.int(n = n, size = floor(test * n * m), replace = TRUE),
              y = sample.int(n = m, size = floor(test * n * m), replace = TRUE))

  mask = matrix(0, nrow = n, ncol = m)
  mask[idx] = 1

  train = (mask == 1)
  test = (mask != 1)

  list(train = train, test = test)
}

## Matrix completion for gllvm and glmpca
matrix.completion <- function (y, x = NULL, z = NULL, ncomp = 2,
                               family = poisson(), niter = 10) {

  # data dimensions
  n = nrow(y)
  m = ncol(y)
  d = ncomp
  p = if (is.null(x)) 0 else ncol(x)
  q = if (is.null(z)) 0 else ncol(z)

  # Initialization via iterated least squares
  init = sgdGMF:::ls.svd.init(Y = y, X = x, Z = z, d = d, family = family)

  # create proxy model matrices if they are null
  if (p == 0) x = matrix(0, nrow = n, ncol = p)
  if (q == 0) z = matrix(0, nrow = q, ncol = m)

  # compute the linear predictor
  eta = matrix(0, nrow = n, ncol = m)
  eta[] = eta + tcrossprod(init$u, init$v)
  eta[] = eta + tcrossprod(x, init$bx)
  eta[] = eta + tcrossprod(init$bz, z)

  # compute the estimated mean, variance and deviance residuals
  mu = family$linkinv(eta)
  var = family$variance(mu)

  # fill the missing values
  isna = is.na(y)
  y[isna] = round(mu[isna])

  # compute the deviance and Pearson residuals
  dr = family$dev.resids(y, mu, 1)
  pr = (y - mu) / sqrt(var)

  # return the output
  # list(n = n, m = m, p = p, q = q,
  #      u = init$u, v = init$v,
  #      beta.x = init$bx, beta.z = init$bz,
  #      y = y, eta = eta, mu = mu, var = var,
  #      dr = dr, pr = pr, isna = isna)

  # return the output
  return (y)
}

naive.completion = function (y) {
  apply(y, 2, function (x) {
    m = mean(x, na.rm = TRUE)
    x[is.na(x)] = round(m)
    return (x)
  })
}

## PREPROCESSING ----

null.residuals = function (
    y, x = NULL, z = NULL,
    family = poisson(), type = "deviance"
) {
  n = nrow(y); m = ncol(y)
  y = as.matrix(y)
  yr = rowMeans(y) # row-average
  yc = colMeans(y) # col-average
  yt = mean(yr) # total average
  mu = tcrossprod(yr, yc) / yt

  res = matrix(NA, nrow = n, ncol = m)
  if (type == "deviance") {
    res[] = sign(y - mu) * sqrt(family$dev.resids(y, mu, 1))
  }
  if (type == "pearson") {
    res[] = (y - mu) / sqrt(family$variance(mu))
  }
  if (!is.null(x)) {
    beta = matrix(NA, nrow = m, ncol = ncol(x))
    beta[] = t(solve(crossprod(x, x), crossprod(x, res)))
    res[] = res - tcrossprod(x, beta)
  }
  if (!is.null(z)) {
    gamma = matrix(NA, nrow = n, ncol = ncol(z))
    gamma[] = t(solve(crossprod(z, z), crossprod(z, t(res))))
    res[] = res - tcrossprod(gamma, z)
  }
  return (res)
}

## MODEL FIT ----

# Latent feature extraction via Pearson residuals
fit.pearson = function (y, x = NULL, z = NULL, ncomp = 2,
                        family = poisson(), verbose = FALSE) {

  # model fitting
  time0 = proc.time()
  res = scry::nullResiduals(object = as.matrix(y), fam = "poisson", type = "pearson")

  if (is.null(x)) {
    beta.x = 0
  } else {
    beta.x = t(solve(crossprod(x), crossprod(x, res)))
    res = res - res - crossprod(x, beta.x)
  }

  if (is.null(z)) {
    beta.z = 0
  } else {
    beta.z = t(solve(crossprod(z), crossprod(z, t(res))))
    res = res - crossprod(beta.z, z)
  }

  # if (!is.null(x)) res = res - x %*% solve(crossprod(x), crossprod(x, res))
  SVD = svd::propack.svd(res, neig = ncomp)
  timef = proc.time()

  # if (is.null(X)) beta.x = matrix(0, nrow = ncol(y), ncol = 0)
  # if (is.null(z)) beta.z = matrix(0, nrow = nrow(y), ncol = 0)

  eta = tcrossprod(SVD$u, SVD$v %*% diag(SVD$d))
  mu = family$linkinv(eta)
  tsne = Rtsne::Rtsne(SVD$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "Pearson",
    u = SVD$u %*% diag(sqrt(SVD$d)),
    v = SVD$v %*% diag(sqrt(SVD$d)),
    d = SVD$d,
    bx = beta.x,
    bz = beta.z,
    eta = eta,
    mu = mu,
    tsne = tsne$Y,
    dev = -1,
    error = -1,
    time = timef - time0)
}

# Latent feature extraction via deviance residuals
fit.deviance = function (y, x = NULL, z = NULL, ncomp = 2,
                         family = poisson(), verbose = FALSE) {

  # model fitting
  time0 = proc.time()
  res = scry::nullResiduals(object = as.matrix(y), fam = "poisson", type = "deviance")

  if (is.null(x)) {
    beta.x = 0
  } else {
    beta.x = t(solve(crossprod(x), crossprod(x, res)))
    res = res - res - crossprod(x, beta.x)
  }

  if (is.null(z)) {
    beta.z = 0
  } else {
    beta.z = t(solve(crossprod(z), crossprod(z, t(res))))
    res = res - crossprod(beta.z, z)
  }

  # if (!is.null(x)) res = res - x %*% solve(crossprod(x), crossprod(x, res))
  SVD = svd::propack.svd(res, neig = ncomp)
  timef = proc.time()

  if (is.null(X)) beta.x = matrix(0, nrow = ncol(y), ncol = 0)
  if (is.null(z)) beta.z = matrix(0, nrow = nrow(y), ncol = 0)

  eta = tcrossprod(SVD$u, SVD$v %*% diag(SVD$d))
  mu = family$linkinv(eta)
  tsne = Rtsne::Rtsne(SVD$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "Deviance",
    u = SVD$u %*% diag(sqrt(SVD$d)),
    v = SVD$v %*% diag(sqrt(SVD$d)),
    d = SVD$d,
    bx = beta.x,
    bz = beta.z,
    eta = eta,
    mu = mu,
    tsne = tsne$Y,
    dev = -1,
    error = -1,
    time = timef - time0)
}

# GLLVM via variational approximations
fit.gllvm = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    maxiter = 1000, stepsize = 0.01, tol = 1e-03,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  # model fitting
  time0 = proc.time()
  fit = gllvm::gllvm(
    y = y, X = x, Z = z,
    formula = ~ .,
    num.lv = ncomp,
    family = family,
    method = "EVA",
    reltol = tol) %>%
    suppressWarnings()
  timef = proc.time()

  # Latent factor normalization
  uv = tcrossprod(fit$lvs, fit$params$theta)
  uv = svd::propack.svd(uv, neig = ncomp)

  # Fitted values
  eta = residuals(fit)$linpred
  mu = family$linkinv(eta)

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Explained deviance metrics
  # rss  =  RSS(y, mu)
  # rmse = RMSE(y, mu)
  # cosd = COSD(y, mu)
  # dev  = RDEV(y, mu, family = family)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "gllvm",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = fit$coefZ,
    bz = fit$coefX,
    eta = eta,
    mu = mu,
    tsne = tsne$Y,
    error = error,
    time = timef - time0)
}

# GLMPCA via AVAGRAD
fit.glmpca = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    maxiter = 1000, stepsize = 0.01, tol = 1e-06,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  # glmPCA model fitting
  time0 = proc.time()
  fit = glmpca::glmpca(
    Y = y, X = z, Z = x,
    L = ncomp,
    fam = "poi",
    minibatch = "none",
    optimizer = "fisher",
    ctl = list(verbose = verbose,
               maxIter = maxiter,
               tol = tol)) %>%
    suppressWarnings()
  timef = proc.time()

  # Latent factor normalization
  uv = tcrossprod(as.matrix(fit$loadings), as.matrix(fit$factors))
  uv = svd::propack.svd(uv, neig = ncomp)

  uv = list(d = rep(1, ncomp),
            u = as.matrix(fit$loadings),
            v = as.matrix(fit$factors))

  # Fitted values
  mu = predict(fit)
  eta = family$linkfun(mu)

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Explained deviance metrics
  # rss  =  RSS(y, mu)
  # rmse = RMSE(y, mu)
  # cosd = COSD(y, mu)
  # dev  = RDEV(y, mu, family = family)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "glmPCA",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = fit$coefZ,
    bz = fit$coefX,
    eta = eta,
    mu = mu,
    tsne = tsne$Y,
    dev = fit$dev,
    error = error,
    time = timef - time0)
}

## NBWaVE
fit.nbwave = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    maxiter = 100, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  idx = which(apply(y, 2, sd) == 0)
  if (length(idx) > 0) {
    for (j in idx) {
      y[sample(1:n, 1),j] = 1
    }
  }

  # NBWaVE model fitting
  time0 = proc.time()
  fit = NewWave::newFit(
    Y = t(y),
    X = x,
    V = z,
    K = ncomp,
    commondispersion = TRUE,
    verbose = verbose,
    maxiter_optimize = maxiter,
    stop_epsilon = tol,
    children = NCORES) %>%
    suppressWarnings()
  timef = proc.time()

  # Latent factor normalization
  b = t(fit@beta)
  xb = tcrossprod(x, b)

  a = t(fit@gamma)
  az = tcrossprod(a, z)

  uv = fit@W %*% fit@alpha

  # Fitted values
  eta = xb + az + uv
  mu = family$linkinv(eta)

  # Orthogonalization
  uv = svd::propack.svd(uv, neig = ncomp)

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Explained deviance metrics
  # rss  =  RSS(y, mu)
  # rmse = RMSE(y, mu)
  # cosd = COSD(y, mu)
  # dev  = RDEV(y, mu, family = family)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "NBWaVE",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = b,
    bz = a,
    eta = eta,
    mu = mu,
    tsne = tsne$Y,
    dev = NULL,
    error = error,
    time = timef - time0)
}

## NMF
fit.nmf = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    verbose = FALSE, train = NULL, test = NULL) {

  suppressPackageStartupMessages(require("NMF"))

  time0 = proc.time()
  if (verbose) {
    fit = NMF::nmf(x = y, rank = ncomp, method = "brunet", seed = "nndsvd", nrun = 1, .options = "v")
  } else {
    fit = NMF::nmf(x = y, rank = ncomp, method = "brunet", seed = "nndsvd", nrun = 1)
  }
  timef = proc.time()

  mu = fit@fit@W %*% fit@fit@H

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(fit@fit@W, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "NMF",
    u = fit@fit@W,
    v = t(fit@fit@H),
    d = NULL,
    bx = NULL,
    bz = NULL,
    eta = NULL,
    mu = mu,
    tsne = tsne$Y,
    dev = NULL,
    error = error,
    time = timef - time0)
}


## NNLM
fit.nnlm = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 1000, verbose = FALSE, train = NULL, test = NULL) {

  suppressPackageStartupMessages(require("NNLM"))

  time0 = proc.time()
  fit = NNLM::nnmf(
    A = y, k = ncomp, alpha = penalty, beta = penalty, n.threads = NCORES,
    method = "lee", loss = "mkl", max.iter = maxiter, verbose = verbose)
  timef = proc.time()

  mu = fit$W %*% fit$H

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(fit$W, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "NNLM",
    u = fit$W,
    v = t(fit$H),
    d = NULL,
    bx = NULL,
    bz = NULL,
    eta = NULL,
    mu = mu,
    tsne = tsne$Y,
    dev = fit$nkl,
    error = error,
    time = timef - time0)
}


## NNLM
fit.cmf = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 1000, verbose = FALSE, train = NULL, test = NULL) {

  suppressPackageStartupMessages(require("cmfrec"))

  if (!is.null(x)) if (sd(x[,1]) == 0) x = x[, -1, drop = FALSE]
  if (!is.null(z)) if (sd(z[,1]) == 0) z = z[, -1, drop = FALSE]

  time0 = proc.time()
  fit = cmfrec::CMF(
    X = y, U = x, I = z, k = ncomp, nonneg = TRUE,
    user_bias = FALSE, item_bias = FALSE, center = FALSE,
    nthreads = NCORES, niter = maxiter, verbose = verbose)
  timef = proc.time()

  eta = crossprod(fit$matrices$A, fit$matrices$B) + fit$matrices$glob_mean
  mu = eta

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(t(fit$matrices$A), dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "CMF",
    u = t(fit$matrices$A),
    v = fit$matrices$B,
    d = NULL,
    bx = NULL,
    bz = NULL,
    eta = NULL,
    mu = mu,
    tsne = tsne$Y,
    dev = NULL,
    error = error,
    time = timef - time0)
}

## MODEL FIT (R) ----

## GMF via SVDReg
fit.svdreg = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 200, stepsize = 0.01, tol = 1e-04,
    verbose = FALSE, train = NULL, test = NULL) {

  time0 = proc.time()
  fit = sgdGMF::init.param(
    y, x, z, ncomp = ncomp, family = family,
    method = "svd", niter = maxiter, verbose = verbose)
  timef = proc.time()

  # Latent factor normalization
  uv = tcrossprod(fit$U, fit$V)
  uv = svd::propack.svd(uv, neig = ncomp)

  # Fitted values
  eta = tcrossprod(cbind(x, fit$A, fit$U), cbind(fit$B, z, fit$V))
  mu = family$linkinv(eta)

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Error matrix
  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c( RSS(train, mu),  RSS(test, mu)),
      cosd = c(COSD(train, mu), COSD(test, mu)),
      dev  = c(RDEV(train, mu, family = family),
               RDEV(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "SVDReg",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = fit$bx,
    bz = fit$bz,
    eta = eta,
    mu = mu,
    tsne = tsne$Y,
    dev = NULL,
    error = error,
    time = timef - time0)
}

## GMF via AIRWLS
fit.R.airwls = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 200, stepsize = 0.01, tol = 1e-04,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  # model fitting
  time0 = proc.time()
  fit = sgdgmf(
    Y = y, X = x, Z = z,
    ncomp = ncomp,
    family = family,
    method = "airwls",
    penalty = list(u = penalty, v = penalty, b = 0),
    init = list(init = "svd", niter = 10),
    control = list(maxiter = maxiter,
                   stepsize = stepsize,
                   verbose = verbose,
                   tol = tol)) %>%
    suppressWarnings()
  timef = proc.time()

  # Latent factor normalization
  uv = tcrossprod(fit$coef$U, fit$coef$V)
  uv = svd::propack.svd(uv, neig = ncomp)

  # Fitted values
  eta = fit$pred$eta
  mu = fit$pred$mu

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c( RSS(train, mu),  RSS(test, mu)),
      cosd = c(COSD(train, mu), COSD(test, mu)),
      dev  = c(RDEV(train, mu, family = family),
               RDEV(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "AIRWLS",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = fit$coef$betaX,
    bz = fit$coef$betaZ,
    eta = fit$pred$eta,
    mu = fit$pred$mu,
    dev = fit$trace$deviance,
    error = error,
    time = timef - time0)
}

## GMF via quasi-Newton algorithm
fit.R.newton = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 1000, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  # model fitting
  time0 = proc.time()
  fit = sgdgmf(
    Y = y, X = x, Z = z,
    ncomp = ncomp,
    family = family,
    method = "newton",
    penalty = list(u = penalty, v = penalty, b = 0),
    init = list(init = "svd", niter = 10),
    control = list(maxiter = maxiter,
                   stepsize = stepsize,
                   verbose = verbose,
                   frequency = 25,
                   tol = tol)) %>%
    suppressWarnings()
  timef = proc.time()

  # Latent factor normalization
  uv = tcrossprod(fit$coef$U, fit$coef$V)
  uv = svd::propack.svd(uv, neig = ncomp)

  # Fitted values
  eta = fit$pred$eta
  mu = fit$pred$mu

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c( RSS(train, mu),  RSS(test, mu)),
      cosd = c(COSD(train, mu), COSD(test, mu)),
      dev  = c(RDEV(train, mu, family = family),
               RDEV(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "Newton",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = fit$coef$betaX,
    bz = fit$coef$betaZ,
    eta = fit$pred$eta,
    mu = fit$pred$mu,
    dev = fit$trace$deviance,
    error = error,
    time = timef - time0)
}

## Coordinate stochastic gradient descent
fit.R.csgd = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 1000, stepsize = 0.001, tol = 1e-04,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  # model fitting
  time0 = proc.time()
  fit = sgdgmf(
    Y = y, X = x, Z = z,
    ncomp = ncomp,
    family = family,
    method = "c-sgd",
    penalty = list(u = penalty, v = penalty, b = 0),
    init = list(init = "svd", niter = 10),
    control = list(maxiter = maxiter,
                   size = c(5, 5),
                   rate0 = stepsize,
                   verbose = verbose,
                   frequency = 100)) %>%
    suppressWarnings()
  timef = proc.time()

  # Latent factor normalization
  uv = tcrossprod(fit$coef$U, fit$coef$V)
  uv = svd::propack.svd(uv, neig = ncomp)

  # Fitted values
  eta = fit$pred$eta
  mu = fit$pred$mu

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c( RSS(train, mu),  RSS(test, mu)),
      cosd = c(COSD(train, mu), COSD(test, mu)),
      dev  = c(RDEV(train, mu, family = family),
               RDEV(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "C-SGD",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = fit$coef$betaX,
    bz = fit$coef$betaZ,
    eta = fit$pred$eta,
    mu = fit$pred$mu,
    dev = fit$trace$deviance,
    error = error,
    time = timef - time0)
}

## Averaged stochastic gradient descent
fit.R.bsgd = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 5000, stepsize = 0.001, tol = 1e-04,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  # model fitting
  time0 = proc.time()
  fit = sgdgmf(
    Y = y, X = x, Z = z,
    ncomp = ncomp,
    family = family,
    method = "b-sgd",
    penalty = list(u = penalty, v = penalty, b = 0),
    init = list(init = "svd", niter = 10),
    control = list(maxiter = maxiter,
                   size = c(100, 20),
                   rate0 = stepsize,
                   verbose = verbose,
                   frequency = 250)) %>%
    suppressWarnings()
  timef = proc.time()

  # Latent factor normalization
  uv = tcrossprod(fit$coef$U, fit$coef$V)
  uv = svd::propack.svd(uv, neig = ncomp)

  # Fitted values
  eta = fit$pred$eta
  mu = fit$pred$mu

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c( RSS(train, mu),  RSS(test, mu)),
      cosd = c(COSD(train, mu), COSD(test, mu)),
      dev  = c(RDEV(train, mu, family = family),
               RDEV(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "B-SGD",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = fit$coef$betaX,
    bz = fit$coef$betaZ,
    eta = fit$pred$eta,
    mu = fit$pred$mu,
    dev = fit$trace$deviance,
    error = error,
    time = timef - time0)
}


## Averaged stochastic gradient descent
fit.R.msgd = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 100, stepsize = 0.001, tol = 1e-04,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  # model fitting
  time0 = proc.time()
  fit = sgdgmf(
    Y = y, X = x, Z = z,
    ncomp = ncomp,
    family = family,
    method = "m-sgd",
    penalty = list(u = penalty, v = penalty, b = 0),
    init = list(init = "svd", niter = 10),
    control = list(maxiter = maxiter,
                   size = 100,
                   rate0 = stepsize,
                   verbose = verbose,
                   frequency = 10)) %>%
    suppressWarnings()
  timef = proc.time()

  # Latent factor normalization
  uv = tcrossprod(fit$coef$U, fit$coef$V)
  uv = svd::propack.svd(uv, neig = ncomp)

  # Fitted values
  eta = fit$pred$eta
  mu = fit$pred$mu

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c( RSS(train, mu),  RSS(test, mu)),
      cosd = c(COSD(train, mu), COSD(test, mu)),
      dev  = c(RDEV(train, mu, family = family),
               RDEV(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "M-SGD",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = fit$coef$betaX,
    bz = fit$coef$betaZ,
    eta = fit$pred$eta,
    mu = fit$pred$mu,
    dev = fit$trace$deviance,
    error = error,
    time = timef - time0)
}


## MODEL FIT (C++) ----


## GMF via AIRWLS
fit.C.airwls = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 200, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  p = ncol(x)
  q = ncol(z)

  familyname = family$family
  linkname = family$link

  if (familyname == "Negative Binomial") {
    familyname = "negbinom"
  }

##   isna = is.na(y)
##
##   yc = apply(y, 2, function (x) {
##     m = mean(x, na.rm = TRUE)
##     x[is.na(x)] = round(m)
##     return (x)
##   })
##
##   r = log(yc + 0.1)
##
##   # r = matrix(NA, nrow = nrow(y), ncol = ncol(y))
##   # r[!isna] = log(y[!isna] + 0.1)
##   # r[isna] = mean(r[!isna])
##
##   # Compute the initial column-specific regression parameters (if any)
##   B = t(solve(crossprod(x), crossprod(x, r)))
##   XB = tcrossprod(x, B)
##
##   # Compute the initial row-specific regression parameter (if any)
##   A = t(solve(crossprod(z), crossprod(z, t(r - XB))))
##   AZ = tcrossprod(A, z)
##
##   # Compute the initial latent factors via incomplete SVD
##   s = svd::propack.svd(r - XB - AZ, neig = ncomp)
##   U = s$u %*% diag(sqrt(s$d))
##   V = s$v %*% diag(sqrt(s$d))
##   UV = tcrossprod(U, V)
##
##   niter = 10
##   for (iter in 1:niter) {
##     r[isna] = (XB + AZ + UV)[isna]
##
##     XtX = crossprod(x)
##     Xty = crossprod(x, r - AZ - UV)
##     B = t(solve(XtX, Xty))
##     XB = tcrossprod(x, B)
##
##     ZtZ = crossprod(z)
##     Zty = crossprod(z, t(r - XB - UV))
##     A = t(solve(ZtZ, Zty))
##     AZ = tcrossprod(A, z)
##
##     s = svd::propack.svd(r - XB - AZ, neig = ncomp)
##     U = s$u %*% diag(sqrt(s$d))
##     V = s$v %*% diag(sqrt(s$d))
##     UV = tcrossprod(U, V)
##   }

  init = sgdGMF::init.param.svd(
    Y = y, X = x, Z = z, ncomp = ncomp,
    family = family, niter = 10, FALSE)

  # model fitting
  time0 = proc.time()
  lambda = (1 - 1e-04) * c(0, 0, penalty, 0) + 1e-04
  fit = sgdGMF::cpp.fit.airwls(
    Y = y, X = x, B = init$B, A = init$A, Z = z, U = init$U, V = init$V,
    ncomp = ncomp, familyname = familyname, linkname = linkname,
    lambda = lambda, maxiter = maxiter, nsteps = 1, stepsize = stepsize,
    eps = 1e-08, nafill = 1, tol = tol, damping = 1e-03,
    verbose = verbose, frequency = 25, parallel = TRUE, nthreads = NCORES)
  timef = proc.time()

  bz = fit$U[, (p+1):(p+q)]
  bx = fit$V[, 1:p]
  u = fit$U[, (p+q+1):(p+q+ncomp)]
  v = fit$V[, (p+q+1):(p+q+ncomp)]

  # Latent factor normalization
  uv = tcrossprod(u, v)
  uv = svd::propack.svd(uv, neig = ncomp)

  mu = fit$mu

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "AIRWLS",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = bx,
    bz = bz,
    eta = fit$eta,
    mu = fit$mu,
    phi = fit$phi,
    tsne = tsne$Y,
    dev = fit$trace[,2],
    error = error,
    time = timef - time0)
}


## GMF via qiasi-Newton
fit.C.newton = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 200, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  p = ncol(x)
  q = ncol(z)

  # init = gmf.init(y, x, z, d = ncomp, method = "svd",
  #                 niter = 0, verbose = FALSE)

  familyname = family$family
  linkname = family$link

  if (familyname == "Negative Binomial") {
    familyname = "negbinom"
  }

##  isna = is.na(y)
##
##  yc = apply(y, 2, function (x) {
##    m = mean(x, na.rm = TRUE)
##    x[is.na(x)] = round(m)
##    return (x)
##  })
##
##  r = log(yc + 0.1)
##
##  # r = matrix(NA, nrow = nrow(y), ncol = ncol(y))
##  # r[!isna] = log(y[!isna] + 0.1)
##  # r[isna] = mean(r[!isna])
##
##  # Compute the initial column-specific regression parameters (if any)
##  B = t(solve(crossprod(x), crossprod(x, r)))
##  XB = tcrossprod(x, B)
##
##  # Compute the initial row-specific regression parameter (if any)
##  A = t(solve(crossprod(z), crossprod(z, t(r - XB))))
##  AZ = tcrossprod(A, z)
##
##  # Compute the initial latent factors via incomplete SVD
##  s = svd::propack.svd(r - XB - AZ, neig = ncomp)
##  U = s$u %*% diag(sqrt(s$d))
##  V = s$v %*% diag(sqrt(s$d))
##  UV = tcrossprod(U, V)
##
##  niter = 10
##  for (iter in 1:niter) {
##    r[isna] = (XB + AZ + UV)[isna]
##
##    XtX = crossprod(x)
##    Xty = crossprod(x, r - AZ - UV)
##    B = t(solve(XtX, Xty))
##    XB = tcrossprod(x, B)
##
##    ZtZ = crossprod(z)
##    Zty = crossprod(z, t(r - XB - UV))
##    A = t(solve(ZtZ, Zty))
##    AZ = tcrossprod(A, z)
##
##    s = svd::propack.svd(r - XB - AZ, neig = ncomp)
##    U = s$u %*% diag(sqrt(s$d))
##    V = s$v %*% diag(sqrt(s$d))
##    UV = tcrossprod(U, V)
##  }

  init = sgdGMF::init.param.svd(
    Y = y, X = x, Z = z, ncomp = ncomp,
    family = family, niter = 10, FALSE)

  # model fitting
  time0 = proc.time()
  fit = sgdGMF::cpp.fit.newton(
    Y = y, X = x, B = init$B, A = init$A, Z = z, U = init$U, V = init$V,
    ncomp = ncomp, familyname = familyname, linkname = linkname,
    lambda = c(0,0,1,0), maxiter = maxiter, stepsize = stepsize,
    eps = 1e-08, nafill = 1, tol = tol, damping = 1e-03,
    verbose = verbose, frequency = 25, parallel = TRUE, nthreads = NCORES)
  timef = proc.time()

  bz = fit$U[, (p+1):(p+q)]
  bx = fit$V[, 1:p]
  u = fit$U[, (p+q+1):(p+q+ncomp)]
  v = fit$V[, (p+q+1):(p+q+ncomp)]

  # Latent factor normalization
  uv = tcrossprod(u, v)
  uv = svd::propack.svd(uv, neig = ncomp)

  mu = fit$mu

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "Newton",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = bx,
    bz = bz,
    eta = fit$eta,
    mu = fit$mu,
    phi = fit$phi,
    tsne = tsne$Y,
    dev = fit$trace[,2],
    error = error,
    time = timef - time0)
}


## GMF via memoized SGD
fit.C.msgd = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 200, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  p = ncol(x)
  q = ncol(z)

  # init = gmf.init(y, x, z, d = ncomp, method = "svd",
  #                 niter = 0, verbose = FALSE)

  familyname = family$family
  linkname = family$link

  if (familyname == "Negative Binomial") {
    familyname = "negbinom"
  }

##  isna = is.na(y)
##
##  yc = apply(y, 2, function (x) {
##    m = mean(x, na.rm = TRUE)
##    x[is.na(x)] = round(m)
##    return (x)
##  })
##
##  r = log(yc + 0.1)
##
##  # r = matrix(NA, nrow = nrow(y), ncol = ncol(y))
##  # r[!isna] = log(y[!isna] + 0.1)
##  # r[isna] = mean(r[!isna])
##
##  # Compute the initial column-specific regression parameters (if any)
##  B = t(solve(crossprod(x), crossprod(x, r)))
##  XB = tcrossprod(x, B)
##
##  # Compute the initial row-specific regression parameter (if any)
##  A = t(solve(crossprod(z), crossprod(z, t(r - XB))))
##  AZ = tcrossprod(A, z)
##
##  # Compute the initial latent factors via incomplete SVD
##  s = svd::propack.svd(r - XB - AZ, neig = ncomp)
##  U = s$u %*% diag(sqrt(s$d))
##  V = s$v %*% diag(sqrt(s$d))
##  UV = tcrossprod(U, V)
##
##  niter = 10
##  for (iter in 1:niter) {
##    r[isna] = (XB + AZ + UV)[isna]
##
##    XtX = crossprod(x)
##    Xty = crossprod(x, r - AZ - UV)
##    B = t(solve(XtX, Xty))
##    XB = tcrossprod(x, B)
##
##    ZtZ = crossprod(z)
##    Zty = crossprod(z, t(r - XB - UV))
##    A = t(solve(ZtZ, Zty))
##    AZ = tcrossprod(A, z)
##
##    s = svd::propack.svd(r - XB - AZ, neig = ncomp)
##    U = s$u %*% diag(sqrt(s$d))
##    V = s$v %*% diag(sqrt(s$d))
##    UV = tcrossprod(U, V)
##  }

  init = sgdGMF::init.param.svd(
    Y = y, X = x, Z = z, ncomp = ncomp,
    family = family, niter = 10, FALSE)

  # model fitting
  time0 = proc.time()
  fit = sgdGMF::cpp.fit.msgd(
    Y = y, X = x, B = init$B, A = init$A, Z = z, U = init$U, V = init$V,
    ncomp = ncomp, familyname = familyname, linkname = linkname,
    lambda = c(0,0,1,0), maxiter = maxiter, rate0 = stepsize,
    size = 100, burn = .9,
    eps = 1e-08, nafill = 1, tol = tol, damping = 1e-03,
    verbose = verbose, frequency = 10)
  timef = proc.time()

  bz = fit$U[, (p+1):(p+q)]
  bx = fit$V[, 1:p]
  u = fit$U[, (p+q+1):(p+q+ncomp)]
  v = fit$V[, (p+q+1):(p+q+ncomp)]

  # Latent factor normalization
  uv = tcrossprod(u, v)
  uv = svd::propack.svd(uv, neig = ncomp)

  mu = fit$mu

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "M-SGD",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = bx,
    bz = bz,
    eta = fit$eta,
    mu = fit$mu,
    phi = fit$phi,
    tsne = tsne$Y,
    dev = fit$trace[,2],
    error = error,
    time = timef - time0)
}

## GMF via coordinate SGD
fit.C.csgd = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 200, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  p = ncol(x)
  q = ncol(z)

  familyname = family$family
  linkname = family$link

  if (familyname == "Negative Binomial") {
    familyname = "negbinom"
  }

  init = sgdGMF::init.param.svd(
    Y = y, X = x, Z = z, ncomp = ncomp,
    family = family, niter = 10, FALSE)

  # model fitting
  time0 = proc.time()
  fit = sgdGMF::cpp.fit.csgd(
    Y = y, X = x, B = init$B, A = init$A, Z = z, U = init$U, V = init$V,
    ncomp = ncomp, familyname = familyname, linkname = linkname,
    lambda = c(0,0,1,0), maxiter = maxiter, rate0 = stepsize,
    size1 = 100, size2 = 10,
    eps = 1e-08, nafill = 1, tol = tol, damping = 1e-03,
    verbose = verbose, frequency = 100, parallel = FALSE)
  timef = proc.time()

  bz = fit$U[, (p+1):(p+q)]
  bx = fit$V[, 1:p]
  u = fit$U[, (p+q+1):(p+q+ncomp)]
  v = fit$V[, (p+q+1):(p+q+ncomp)]

  # Latent factor normalization
  uv = tcrossprod(u, v)
  uv = svd::propack.svd(uv, neig = ncomp)

  mu = fit$mu

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "C-SGD",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = bx,
    bz = bz,
    eta = fit$eta,
    mu = fit$mu,
    phi = fit$phi,
    tsne = tsne$Y,
    dev = fit$trace[,2],
    error = error,
    time = timef - time0)
}


## GMF via row-wise SGD
fit.C.rsgd = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 200, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  p = ncol(x)
  q = ncol(z)

  # init = gmf.init(y, x, z, d = ncomp, method = "svd",
  #                 niter = 0, verbose = FALSE)

  familyname = family$family
  linkname = family$link

  if (familyname == "Negative Binomial") {
    familyname = "negbinom"
  }

##  isna = is.na(y)
##
##  yc = apply(y, 2, function (x) {
##    m = mean(x, na.rm = TRUE)
##    x[is.na(x)] = round(m)
##    return (x)
##  })
##
##  r = log(yc + 0.1)
##
##  # r = matrix(NA, nrow = nrow(y), ncol = ncol(y))
##  # r[!isna] = log(y[!isna] + 0.1)
##  # r[isna] = mean(r[!isna])
##
##  # Compute the initial column-specific regression parameters (if any)
##  B = t(solve(crossprod(x), crossprod(x, r)))
##  XB = tcrossprod(x, B)
##
##  # Compute the initial row-specific regression parameter (if any)
##  A = t(solve(crossprod(z), crossprod(z, t(r - XB))))
##  AZ = tcrossprod(A, z)
##
##  # Compute the initial latent factors via incomplete SVD
##  s = svd::propack.svd(r - XB - AZ, neig = ncomp)
##  U = s$u %*% diag(sqrt(s$d))
##  V = s$v %*% diag(sqrt(s$d))
##  UV = tcrossprod(U, V)
##
##  niter = 10
##  for (iter in 1:niter) {
##    r[isna] = (XB + AZ + UV)[isna]
##
##    XtX = crossprod(x)
##    Xty = crossprod(x, r - AZ - UV)
##    B = t(solve(XtX, Xty))
##    XB = tcrossprod(x, B)
##
##    ZtZ = crossprod(z)
##    Zty = crossprod(z, t(r - XB - UV))
##    A = t(solve(ZtZ, Zty))
##    AZ = tcrossprod(A, z)
##
##    s = svd::propack.svd(r - XB - AZ, neig = ncomp)
##    U = s$u %*% diag(sqrt(s$d))
##    V = s$v %*% diag(sqrt(s$d))
##    UV = tcrossprod(U, V)
##  }

  init = sgdGMF::init.param.svd(
    Y = y, X = x, Z = z, ncomp = ncomp,
    family = family, niter = 10, FALSE)

  # model fitting
  time0 = proc.time()
  fit = sgdGMF::cpp.fit.rsgd(
    Y = y, X = x, B = init$B, A = init$A, Z = z, U = init$U, V = init$V,
    ncomp = ncomp, familyname = familyname, linkname = linkname,
    lambda = c(0,0,1,0), maxiter = maxiter, burn = 1, rate0 = stepsize,
    size1 = 100, size2 = m, eps = 1e-08, nafill = 1, tol = tol,
    damping = 1e-03, verbose = verbose, frequency = 250)
  timef = proc.time()

  bz = fit$U[, (p+1):(p+q)]
  bx = fit$V[, 1:p]
  u = fit$U[, (p+q+1):(p+q+ncomp)]
  v = fit$V[, (p+q+1):(p+q+ncomp)]

  # Latent factor normalization
  uv = tcrossprod(u, v)
  uv = svd::propack.svd(uv, neig = ncomp)

  mu = fit$mu

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "R-SGD",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = bx,
    bz = bz,
    eta = fit$eta,
    mu = fit$mu,
    phi = fit$phi,
    tsne = tsne$Y,
    dev = fit$trace[,2],
    error = error,
    time = timef - time0)
}


## GMF via block-wise SGD
fit.C.bsgd = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 200, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  p = ncol(x)
  q = ncol(z)

  # init = gmf.init(y, x, z, d = ncomp, method = "svd",
  #                 niter = 0, verbose = FALSE)

  familyname = family$family
  linkname = family$link

  if (familyname == "Negative Binomial") {
    familyname = "negbinom"
  }

  init = sgdGMF::init.param.svd(
    Y = y, X = x, Z = z, ncomp = ncomp,
    family = family, niter = 10, FALSE)

  # model fitting
  time0 = proc.time()
  fit = sgdGMF::cpp.fit.bsgd(
    Y = y, X = x, B = init$B, A = init$A, Z = z, U = init$U, V = init$V,
    ncomp = ncomp, familyname = familyname, linkname = linkname,
    lambda = c(0,0,1,0), maxiter = maxiter, rate0 = stepsize,
    size1 = 100, size2 = 20, eps = 1e-08, nafill = 1, tol = tol,
    damping = 1e-03, verbose = verbose, frequency = 100)
  timef = proc.time()

  bz = fit$U[, (p+1):(p+q)]
  bx = fit$V[, 1:p]
  u = fit$U[, (p+q+1):(p+q+ncomp)]
  v = fit$V[, (p+q+1):(p+q+ncomp)]

  # Latent factor normalization
  uv = tcrossprod(u, v)
  uv = svd::propack.svd(uv, neig = ncomp)

  mu = fit$mu

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "B-SGD",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = bx,
    bz = bz,
    eta = fit$eta,
    mu = fit$mu,
    phi = fit$phi,
    tsne = tsne$Y,
    dev = fit$trace[,2],
    error = error,
    time = timef - time0)
}




