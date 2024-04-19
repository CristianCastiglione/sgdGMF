
#' @title S3 class: initgmf
#'
#' @description
#' A short description...
#'
#' @slot method initialization method
#' @slot family model family
#' @slot ncomp number of latent components
#' @slot type residual type
#' @slot verbose if \code{TRUE}, ...
#' @slot parallel if \code{TRUE}, ...
#' @slot nthreads ...
#' @slot A matrix of row-specific regression effects (\eqn{n \times q})
#' @slot B matrix of column-specific regression effects (\eqn{m \times p})
#' @slot U matrix of latent scores (\eqn{n \times d})
#' @slot V matrix of factor loadings (\eqn{m \times d})
#' @slot phi scalar dispersion parameter
#'
#' @export
setClass("initgmf",
  slots = list(
    method = "character",
    family = "list",
    ncomp = "numeric",
    type = "character",
    verbose = "logical",
    parallel = "logical",
    nthreads = "numeric",
    A = "matrix",
    B = "matrix",
    U = "matrix",
    V = "matrix",
    phi = "vector"
))

#' @method deviance initgmf
#' @export
deviance.initgmf = function (object, Y, X, Z, normalize = FALSE) {
  U = cbind(X, object$A, object$U)
  V = cbind(object$B, Z, object$V)
  eta = tcrossprod(U, V)
  mu = object$family$linkinv(eta)
  dev = matrix.deviance(mu, Y, object$family)
  if (normalize) {
    mu0 = matrix(mean(Y, na.rm = TRUE), nrow(Y), ncol(Y))
    dev0 = matrix.deviance(mu0, Y, object$family)
    dev = dev / dev0
  }
  return (dev)
}

#' @method coefficients initgmf
#' @export
coefficients.initgmf = function (
    object, type = c("all", "colreg", "rowreg", "scores", "loadings")
) {
  type = match.arg(type)
  switch(type,
    "colreg" = object$B,
    "rowreg" = object$A,
    "scores" = object$U,
    "loadings" = object$V,
    "all" = list(colcoef = object$B, rowcoef = object$A,
                 scores = object$U, loadings = object$V))
}

#' @method residuals initgmf
#' @export
residuals.initgmf = function (
    object, Y, X = NULL, Z = NULL,
    type = c("deviance", "pearson", "working", "response", "link"),
    partial = FALSE, normalize = FALSE, fillna = FALSE, spectrum = FALSE, ncomp = 50
) {
  # Set the residual type
  type = match.arg(type)

  # Set the covariate matrices
  if (is.null(X)) X = matrix(1, nrow = nrow(Y), ncol = 1)
  if (is.null(Z)) Z = matrix(1, nrow = ncol(Y), ncol = 1)

  # Safety checks
  if (!is.numeric(X) | !is.matrix(X)) stop("`X` is not a numeric matrix.")
  if (!is.numeric(Z) | !is.matrix(Z)) stop("`Z` is not a numeric matrix.")

  if (nrow(X) != nrow(object$U)) stop("Incompatible dimensions.")
  if (ncol(X) != ncol(object$B)) stop("Incompatible dimensions.")

  if (nrow(Z) != nrow(object$V)) stop("Incompatible dimensions.")
  if (ncol(Z) != ncol(object$A)) stop("Incompatible dimensions.")

  # Compute the predicted values
  family = object$family
  if (partial) {
    U = cbind(X, object$A)
    V = cbind(object$B, Z)
  } else {
    U = cbind(X, object$A, object$U)
    V = cbind(object$B, Z, object$V)
  }
  eta = tcrossprod(U, V)
  mu = family$linkinv(eta)

  # Compute the residuals
  res = switch(type,
    "deviance" = sign(Y - mu) * sqrt(abs(family$dev.resids(Y, mu, 1))),
    "pearson" = (Y - mu) / sqrt(abs(family$variance(mu))),
    "working" = (Y - mu) * family$mu.eta(mu) / abs(family$variance(mu)),
    "response" = (Y - mu),
    "link" = (family$transform(Y) - eta))

  # Fill the missing values using Gaussian random values
  if (anyNA(res) & (fillna | spectrum)) {
    res = apply(res, 2, function (x) {
      if (anyNA(x)) {
        na = which(is.na(x) | is.nan(x))
        r = length(na)
        m = mean(x, na.rm = TRUE)
        s = sd(x, na.rm = TRUE)
        x[na] = rnorm(r, mean = m, sd = s)
      }
      return (x)
    })
  }

  # Standardize the residuals column-by-column
  if (normalize) {
    res = scale(res, center = TRUE, scale = TRUE)
  }

  # Decompose the residuals using incomplete SVD
  if (spectrum) {
    rcov = cov(res)
    ncomp = max(1, min(ncomp, ncol(res)))
    pca = RSpectra::eigs_sym(rcov, ncomp)

    # Estimate the explained and residual variance
    var.eig = pca$values
    var.tot = sum(diag(rcov))
    var.exp = sum(var.eig)
    var.res = var.tot - var.exp
  }

  # Return the residuals and the corresponding spectrum
  if (!spectrum) {
    return (res)
  } else {
    return (
      list(residuals = res,
           lambdas = var.eig,
           explained.var = var.exp,
           reminder.var = var.res,
           total.var = var.tot))
  }
}

#' @method fitted initgmf
#' @export
fitted.initgmf = function (
    object, X = NULL, Z = NULL,
    type = c("link", "response", "terms"), partial = FALSE
) {
  # Set the fitted value type
  type = match.arg(type)

  # Set the covariate matrices
  if (is.null(X)) X = matrix(1, nrow = nrow(Y), ncol = 1)
  if (is.null(Z)) Z = matrix(1, nrow = ncol(Y), ncol = 1)

  # Safety checks
  if (!is.numeric(X) | !is.matrix(X)) stop("`X` is not a numeric matrix.")
  if (!is.numeric(Z) | !is.matrix(Z)) stop("`Z` is not a numeric matrix.")

  if (nrow(X) != nrow(object$U)) stop("Incompatible dimensions.")
  if (ncol(X) != ncol(object$B)) stop("Incompatible dimensions.")

  if (nrow(Z) != nrow(object$V)) stop("Incompatible dimensions.")
  if (ncol(Z) != ncol(object$A)) stop("Incompatible dimensions.")

  # Return the fitted values depending on the prediction type
  XB = tcrossprod(X, object$B)
  AZ = tcrossprod(object$A, Z)
  if (!partial) {
    UV = tcrossprod(object$U, object$V)
    switch(type,
      "link" = XB + AZ + UV,
      "response" = object$family$linkinv(XB + AZ + UV),
      "terms" = list(XB = XB, AZ = AZ, UV = UV))
  } else {
    switch(type,
      "link" = XB + AZ,
      "response" = object$family$linkinv(XB + AZ),
      "terms" = list(XB = XB, AZ = AZ, UV = NULL))
  }
}

#' @method plot initgmf
#' @export
plot.initgmf = function (
    object, Y, X = NULL, Z = NULL,
    type = c("1", "2", "3", "4", "5", "6", "idx", "fit", "std", "hist", "qq", "ecdf"),
    resid = c("deviance", "pearson", "working", "response", "link"),
    subsample = FALSE, sample.size = 500, partial = FALSE,
    normalize = FALSE, fillna = FALSE, bycol = FALSE
) {
  type = match.arg(type)
  resid = match.arg(resid)

  fit = switch(resid,
    "deviance" = fitted(object, X, Z, type = "response", partial = partial),
    "pearson" = fitted(object, X, Z, type = "response", partial = partial),
    "working" = fitted(object, X, Z, type = "response", partial = partial),
    "response" = fitted(object, X, Z, type = "response", partial = partial),
    "link" = fitted(object, X, Z, type = "link", partial = partial))

  res = residuals(
    object, Y, X, Z, type = resid, partial = partial,
    normalize = normalize, fillna = fillna, spectrum = FALSE)

  if (subsample) {
    n = nrow(object$U)
    m = nrow(object$V)
    if (sample.size < n*m) {
      idx = cbind(
        row = sample.int(n = n, size = sample.size, replace = TRUE),
        col = sample.int(n = m, size = sample.size, replace = TRUE))
      fit = fit[idx]
      res = res[idx]
      col = idx$col
    }
  } else {
    col = expand.grid(row = 1:n, col = 1:m)$col
  }

  if (type %in% c("1", "idx")) {
    df = data.frame(residuals = c(res), index = c(1:prod(dim(res))), column = as.factor(col))
    if (!bycol) plt = ggplot(data = df, map = aes(x = index, y = residuals))
    if ( bycol) plt = ggplot(data = df, map = aes(x = index, y = residuals, color = column))
    plt = plt + geom_point(alpha = 0.5) + geom_hline(yintercept = 0, col = 2, lty = 2) +
      labs(x = "Index", y = "Residuals", title = "Residuals vs Fitted values")
  }
  if (type %in% c("2", "fit")) {
    df = data.frame(residuals = c(res), fitted = c(fit), column = as.factor(col))
    if (!bycol) plt = ggplot(data = df, map = aes(x = fitted, y = residuals))
    if ( bycol) plt = ggplot(data = df, map = aes(x = fitted, y = residuals, color = column))
    plt = plt + geom_point(alpha = 0.5) + geom_hline(yintercept = 0, col = 2, lty = 2) +
      labs(x = "Fitted values", y = "Residuals", title = "Residuals vs Fitted values")
  }
  if (type %in% c("3", "std")) {
    df = data.frame(residuals = c(sqrt(abs(res))), fitted = c(fit), column = as.factor(col))
    if (!bycol) plt = ggplot(data = df, map = aes(x = fitted, y = residuals))
    if ( bycol) plt = ggplot(data = df, map = aes(x = fitted, y = residuals, color = column))
    plt = plt + geom_point(alpha = 0.5) + geom_hline(yintercept = 0, col = 2, lty = 2) +
      labs(x = "Fitted values", y = "Abs. Residuals", title = "Residuals vs Fitted values")
  }
  if (type %in% c("4", "hist")) {
    df = data.frame(residuals = c(res), column = as.factor(col))
    if (!bycol) plt = ggplot(data = df, map = aes(x = residuals, y = after_stat(density)))
    if ( bycol) plt = ggplot(data = df, map = aes(x = residuals, y = after_stat(density), color = column, fill = column))
    plt = plt + geom_histogram(bins = 30) + geom_vline(xintercept = 0, col = 2, lty = 2) +
      labs(x = "Residuals", y = "Frequency", title = "Histogram of the residuals")
  }
  if (type %in% c("5", "qq")) {
    df = list2DF(qqnorm(scale(c(res)), plot.it = FALSE))
    plt = ggplot(data = df, map = aes(x = x, y = y)) +
      geom_abline(intercept = 0, slope = 1, color = 2, lty = 2) + geom_point(alpha = 0.5) +
      labs(x = "Theoretical quantiles", y = "Empirical quantiles", title = "Residual QQ-plot")
  }
  if (type %in% c("6", "ecdf")) {
    zn = scale(c(res))
    zz = seq(from = min(zn), to = max(zn), length = 100)
    df1 = data.frame(x = zn, y = ecdf(zn)(zn))
    df2 = data.frame(x = zz, y = pnorm(zz))
    plt = ggplot() +
      geom_line(data = df2, map = aes(x = x, y = y), color = 2) +
      geom_point(data = df1, map = aes(x = x, y = y), alpha = 0.5) +
      labs(x = "Standardized residuals", y = "Empirical CDF", title = "Residual ECDF plot")
  }

  return (plt)
}

#' @method screeplot initgmf
#' @export
screeplot.initgmf = function (
    object, Y, X = NULL, Z = NULL, ncomp = 20,
    type = c("deviance", "pearson", "working", "response", "link"),
    partial = FALSE, normalize = FALSE,
    cumulative = FALSE, proportion = FALSE
) {

  ncomp = max(1, min(ncomp, nrow(object$V)))
  res = residuals(object, Y = Y, X = X, Z = Z,
                  type = type, partial = partial, normalize = normalize,
                  fillna = TRUE, spectrum = TRUE, ncomp = ncomp)

  lambdas = res$lambdas
  if (cumulative) lambdas = cumsum(lambdas)
  if (proportion) lambdas = lambdas / res$total.var

  df = data.frame(components = 1:ncomp, lambdas = lambdas)
  plt = ggplot(data = df, map = aes(x = components, y = lambdas)) + geom_col() +
    labs(x = "Components", y = "Eigenvalues", title = "Residual screeplot")

  return (plt)
}

#' @method biplot initgmf
#' @export
biplot.initgmf = function (
    object, choices = 1:2, normalize = FALSE,
    labels = NULL, palette = c("viridis", "inferno")
) {
  n = nrow(object$U)
  m = nrow(object$V)

  i = max(1, min(choices[1], object$ncomp))
  j = max(1, min(choices[2], object$ncomp))

  if (is.null(labels)) {
    labels = list(scores = c(1:n), loadings = c(1:m))
  }

  if (normalize) {
    pca = RSpectra::svds(tcrossprod(object$U, object$V), object$ncomp)
    u = scale(pca$u[,c(i,j)])
    v = scale(pca$v[,c(i,j)])
  } else {
    u = scale(object$U[,c(i,j)])
    v = scale(object$V[,c(i,j)])
  }

  scores = data.frame(idx = labels$scores, pc1 = c(u[,1]), pc2 = c(u[,2]))
  loadings = data.frame(idx = labels$loadings, pc1 = c(v[,1]), pc2 = c(v[,2]))

  plt.scores =
    ggplot(data = scores, map = aes(x = pc1, y = pc2, color = idx, label = idx)) +
    geom_hline(yintercept = 0, lty = 2, color = "grey40") +
    geom_vline(xintercept = 0, lty = 2, color = "grey40") +
    geom_point() + geom_text(color = 1, size = 2.5, nudge_x = -0.1, nudge_y = +0.1) +
    scale_color_viridis(option = palette[1]) + theme(legend.position = "bottom") +
    labs(x = paste("PC", i), y = paste("PC", j), color = "Index", title = "Scores")

  plt.loadings =
    ggplot(data = loadings, map = aes(x = pc1, y = pc2, color = idx, label = idx)) +
    geom_hline(yintercept = 0, lty = 2, color = "grey40") +
    geom_vline(xintercept = 0, lty = 2, color = "grey40") +
    geom_point() + geom_text(color = 1, size = 2.5, nudge_x = -0.1, nudge_y = +0.1) +
    scale_color_viridis(option = palette[2]) + theme(legend.position = "bottom") +
    labs(x = paste("PC", i), y = paste("PC", j), color = "Index", title = "Loadings")

  list(scores = plt.scores, loadings = plt.loadings)
}

#' @method heatmap initgmf
#' @export
heatmap.initgmf = function (
    object, Y = NULL, X = NULL, Z = NULL,
    type = c("data", "response", "link", "scores", "loadings", "deviance", "pearson", "working"),
    resid = FALSE, symmetric = FALSE, transpose = FALSE, limits = NULL, palette = NULL
) {
  type = match.arg(type)

  if (resid) {
    if (type == "data") stop("type='data' is not allowed with resid=TRUE", call. = FALSE)
    if (type == "scores") stop("type='scores' is not allowed with resid=TRUE", call. = FALSE)
    if (type == "loadings") stop("type='loadings' is not allowed with resid=TRUE", call. = FALSE)

    mat = residuals(object, Y, X, Z, type = type)
    if (transpose) mat = t(mat)

  } else {
    if (type == "deviance") stop("type='deviance' is not allowed with resid=FALSE", call. = FALSE)
    if (type == "pearson") stop("type='pearson' is not allowed with resid=FALSE", call. = FALSE)
    if (type == "working") stop("type='working' is not allowed with resid=FALSE", call. = FALSE)

    mat = switch(type,
      "data" = Y,
      "response" = {
        U = cbind(X, object$A, object$U)
        V = cbind(object$B, Z, object$V)
        eta = tcrossprod(U, V)
        object$family$linkinv(eta)},
      "link" = {
        U = cbind(X, object$A, object$U)
        V = cbind(object$B, Z, object$V)
        eta = tcrossprod(U, V)},
      "scores" = object$U,
      "loadings" = object$V)
    if (transpose) mat = t(mat)
  }

  df = reshape2::melt(mat, varnames = c("sample", "variable"))

  if (is.null(limits)) {
    limits = range(df$value, na.rm = TRUE)
    if (symmetric) limits = c(-1,+1) * max(abs(limits))
  }

  if (is.null(palette)) {
    palette = viridis(100, option = "viridis")
    if (resid) palette = hcl.colors(100, palette = "RdBu")
    if (symmetric) palette = hcl.colors(100, palette = "RdBu")
  }

  plt = ggplot(data = df, map = aes(x = variable, y = sample, fill = value)) +
    geom_raster() + scale_fill_gradientn(colours = palette, limits = limits) +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    theme(legend.position = "bottom", panel.grid = element_blank()) +
    labs(x = "Variables", y = "Samples", fill = "Intensity")

  return (plt)
}
