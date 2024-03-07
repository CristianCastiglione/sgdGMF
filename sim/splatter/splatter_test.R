
## WORKSPACE SETUP ----

## Clean the workspace
rm(list = ls())
graphics.off()

## Load the package
devtools::load_all()

## Load the utility functions
source("sim/utilities.R")

## DATA SIMULATION ----
SETTING = 1
SAVE = FALSE
SHOW = TRUE

a = 500
b = 10
n = a * b
m = a
sim = NULL
if (SETTING == 1) { # ELLIPTIC GROUPS
  nn = floor(n/3)
  params = splatter::newSplatParams()
  params = splatter::setParam(params, "batchCells", c(nn, nn, n - 2 * nn))
  params = splatter::setParam(params, "nGenes", m)
  params = splatter::setParam(params, "group.prob", c(0.1, 0.2, 0.2, 0.2, 0.3))
  params = splatter::setParam(params, "de.prob", c(0.3, 0.1, 0.2, 0.01, 0.1))
  params = splatter::setParam(params, "de.downProb", c(0.1, 0.4, 0.9, 0.6, 0.5))
  params = splatter::setParam(params, "de.facLoc", c(0.6, 0.1, 0.1, 0.01, 0.2))
  params = splatter::setParam(params, "de.facScale", c(0.1, 0.4, 0.2, 0.5, 0.4))
  # params = splatter::setParam(params, "seed", 140275)
  sim = splatter::splatSimulateGroups(params, verbose = FALSE)
}
if (SETTING == 2) { # LINEAR PATHS
  params = splatter::newSplatParams()
  # params = splatter::setParam(params, "batchCells", n)
  params = splatter::setParam(params, "batchCells", c(n/2, n/2))
  params = splatter::setParam(params, "nGenes", m)
  params = splatter::setParam(params, "group.prob", c(0.1, 0.2, 0.2, 0.2, 0.3))
  params = splatter::setParam(params, "de.prob", 0.5)
  params = splatter::setParam(params, "de.facLoc", 0.2)
  params = splatter::setParam(params, "path.from", c(0, 1, 2, 3, 4))
  sim = splatter::splatSimulatePaths(params, verbose = FALSE)
}
if (SETTING == 3) { # BRANCHING PATHS
  params = splatter::newSplatParams()
  # params = splatter::setParam(params, "batchCells", n)
  params = splatter::setParam(params, "batchCells", c(n/2, n/2))
  params = splatter::setParam(params, "nGenes", m)
  params = splatter::setParam(params, "group.prob", c(0.25, 0.25, 0.25, 0.25))
  params = splatter::setParam(params, "de.prob", 0.5)
  params = splatter::setParam(params, "de.facLoc", 0.2)
  params = splatter::setParam(params, "path.from", c(0, 1, 1, 3))
  sim = splatter::splatSimulatePaths(params, verbose = FALSE)

  params = splatter::newSplatParams()
  params = splatter::setParam(params, "batchCells", c(n/4, n/4, n/4, n/4))
  params = splatter::setParam(params, "nGenes", m)
  params = splatter::setParam(params, "group.prob", c(0.1, 0.2, 0.3, 0.3, 0.1))
  params = splatter::setParam(params, "de.prob", 0.5)
  params = splatter::setParam(params, "de.facLoc", 0.2)
  params = splatter::setParam(params, "batch.facLoc", 0.15)
  params = splatter::setParam(params, "path.from", c(0, 1, 1, 3, 3))
  sim = splatter::splatSimulatePaths(params, verbose = FALSE)
}

sim = scater::logNormCounts(sim)



## DATA EXTRACTION ----
logcounts = as.data.frame(logcounts(sim))
counts = as.data.frame(counts(sim))
cells = as.data.frame(colData(sim))
genes = as.data.frame(rowData(sim))
meta = metadata(sim)

# PCA embedding
log.pca = RSpectra::svds(scale(t(as.matrix(logcounts))), k = 2)
plot(log.pca$u, col = cells$Group)
plot(log.pca$u, col = as.factor(cells$Batch))

# t-SNE embedding
log.tsne = Rtsne::Rtsne(t(as.matrix(logcounts)), dims = 2, verbose = TRUE)
plot(log.tsne$Y, col = cells$Group)
plot(log.tsne$Y, col = as.factor(cells$Batch))

# UMAP embedding
log.umap = umap::umap(t(as.matrix(logcounts)), n_components = 2, verbose = TRUE)
plot(log.umap$layout, col = cells$Group)
plot(log.umap$layout, col = as.factor(cells$Batch))

## TRAIN-TEST SPLIT ----
X = model.matrix(~ Batch, data = cells)
Z = matrix(1, nrow = m, ncol = 1)
Y = matrix(NA, nrow = n, ncol = m)
Y[] = t(counts)

ncomp = 5
family = poisson()

# Sparsification and matrix completion
data = train.test.split(Y, test = 0.3)

test = Y
train = Y
ctrain = Y

test[data$test] = NA
train[data$train] = NA
ctrain = naive.completion(train)
# train = matrix.completion(y = train0, x = X, z = Z, ncomp = ncomp, family = family, niter = 10)

if (SAVE) {
  setting = ifelse(SETTING == 1, "bubble", ifelse(SETTING == 2, "linear", "branch"))
  filename = join.string("example_s", setting, "_n", n, "_m", m, "_d", ncomp, "_i25.RData")
  filepath = join.path("data", "splatter", filename)
  save(n, m, ncomp, family, sim, data, file = filename)
}

## MODEL FIT ----
model.pearson = fit.pearson(y = ctrain, x = X, z = Z, ncomp = ncomp, family = family)
model.deviance = fit.deviance(y = ctrain, x = X, z = Z, ncomp = ncomp, family = family)
model.glmpca = fit.glmpca(y = ctrain, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, tol = 1e-05)
model.nbwave = fit.nbwave(y = ctrain, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, tol = 1e-04)
model.cmf = fit.cmf(y = train, x = X, z = NULL, ncomp = ncomp, family = family, verbose = FALSE, maxiter = 500)
model.nmf = fit.nmf(y = ctrain, x = NULL, z = NULL, ncomp = ncomp, family = family, verbose = TRUE)
model.nnlm = fit.nnlm(y = train, x = NULL, z = NULL, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 2000)
model.svdreg = fit.svdreg(y = train, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 10)
model.airwls = fit.C.airwls(y = train, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, stepsize = 0.9)
model.newton = fit.C.newton(y = train, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, stepsize = 0.2)
model.msgd = fit.C.msgd(y = train, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 100, stepsize = 0.1)
model.csgd = fit.C.csgd(y = train, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 1000, stepsize = 0.01)
model.rsgd = fit.C.rsgd(y = train, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 1500, stepsize = 0.1)
model.bsgd = fit.C.bsgd(y = train, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 1000, stepsize = 0.01)

## MODEL CHECK ----

# Execution time
{
cat("glmPCA: ", model.glmpca$time[3], "s \n")
cat("NBWaVE: ", model.nbwave$time[3], "s \n")
cat("NNLM:   ", model.nnlm$time[3], "s \n")
cat("NMF:    ", model.nmf$time[3], "s \n")
cat("CMF:    ", model.cmf$time[3], "s \n")
cat("SVDreg: ", model.svdreg$time[3], "s \n")
cat("AIRWLS: ", model.airwls$time[3], "s \n")
cat("Newton: ", model.newton$time[3], "s \n")
cat("M-SGD:  ", model.msgd$time[3], "s \n")
cat("C-SGD:  ", model.csgd$time[3], "s \n")
cat("R-SGD:  ", model.rsgd$time[3], "s \n")
cat("B-SGD:  ", model.bsgd$time[3], "s \n")
}

# Convergence check
plot(model.glmpca$dev, type = "b", xlab = "iterations", ylab = "deviance")
plot(model.airwls$dev, type = "b", xlab = "iterations", ylab = "deviance")
plot(model.newton$dev, type = "b", xlab = "iterations", ylab = "deviance")
plot(model.msgd$dev, type = "b", xlab = "iterations", ylab = "deviance")
plot(model.csgd$dev, type = "b", xlab = "iterations", ylab = "deviance")
plot(model.rsgd$dev, type = "b", xlab = "iterations", ylab = "deviance")
plot(model.bsgd$dev, type = "b", xlab = "iterations", ylab = "deviance")

## Reconstruction error
full.sample.error = error.summary(
  Y, cells$Group, model.glmpca, model.nbwave, model.nmf, model.nnlm, model.cmf, # model.svdreg,
  model.airwls, model.newton, model.msgd, model.csgd, model.rsgd, model.bsgd)
print(full.sample.error)

# Out-of-sample error
in.sample.error = error.summary(
  train, cells$Group, model.glmpca, model.nbwave, model.nmf, model.nnlm, model.cmf, # model.svdreg,
  model.airwls, model.newton, model.msgd, model.csgd, model.rsgd, model.bsgd)
print(in.sample.error)

# Out-of-sample error
out.sample.error = error.summary(
  test, cells$Group, model.glmpca, model.nbwave, model.nmf, model.nnlm, model.cmf, # model.svdreg,
  model.airwls, model.newton, model.msgd, model.csgd, model.rsgd, model.bsgd)
print(out.sample.error)

df = data.frame(rbind(full.sample.error, in.sample.error, out.sample.error))
df$Sample = rep(c("full", "train", "test"), each = 11)
df$Model = factor(df$Model, levels = c("glmPCA", "NBWaVE", "NMF", "NNLM", "CMF", "SVDReg", "AIRWLS", "Newton", "M-SGD", "C-SGD", "R-SGD", "B-SGD"))
df$Sample = as.factor(df$Sample)
df$Time = as.numeric(df$Time)
df$RSS = as.numeric(df$RSS)
df$Cos = as.numeric(df$Cos)
df$Dev = as.numeric(df$Dev)

ggplot(data = df, map = aes(x = Model, y = Time, color = Sample, pch = Sample)) + geom_point(size = 4) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggplot(data = df, map = aes(x = Model, y = RSS, color = Sample, pch = Sample)) + geom_point(size = 4) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggplot(data = df, map = aes(x = Model, y = Dev, color = Sample, pch = Sample)) + geom_point(size = 4) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggplot(data = df, map = aes(x = Model, y = Sil, color = Sample, pch = Sample)) + geom_point(size = 4) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggplot(data = df, map = aes(x = log10(Time), y = Dev, color = Model, pch = Sample)) + geom_point(size = 4) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


## TSNE PROJECTION ----
plot.tsne.grid = function (tsne, by = 1) {
  models = c()
  df = data.frame(
    x = c(),
    y = c(),
    group = c(),
    batch = c(),
    model = c())

  for (t in tsne) {
    n = nrow(t$tsne)
    dft = data.frame(
      x = t$tsne[,1],
      y = t$tsne[,2],
      group = t$group,
      batch = t$batch,
      model = rep(t$model, n))
    idx = seq(1, n, by = by)
    df = rbind(df, dft[idx, ])
    models = c(models, t$model)
  }

  df$x = as.numeric(df$x)
  df$y = as.numeric(df$y)
  df$group = as.factor(df$group)
  df$batch = as.factor(df$batch)
  df$model = factor(df$model, levels = models)

  levels(df$group) = 1:length(unique(df$group))
  levels(df$batch) = 1:length(unique(df$batch))

  colors = c("#ff7c04", "#387cbc", "#e81c1c", "#50ac4c", "#a04ca4")
  plt = ggplot(data = df, mapping = aes(x = x, y = y, color = group)) +
    geom_point(alpha = 0.3) + facet_wrap(vars(model)) +
    labs(color = "Cell-type") + theme_minimal() +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1.0))) +
    theme(text = element_text(size = 20),
          legend.position = "bottom",
          axis.title = element_blank(),
          axis.text = element_blank())

  return (plt)
}

plt = plot.tsne.grid(list(
  list(model = "Pearson", tsne = model.pearson$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "Deviance", tsne = model.deviance$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "NNLM", tsne = model.nnlm$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "NMF", tsne = model.nmf$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "CMF", tsne = model.cmf$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "glmPCA", tsne = model.glmpca$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "NBWaVE", tsne = model.nbwave$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "SVDreg", tsne = model.svdreg$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "AIRWLS", tsne = model.airwls$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "Newton", tsne = model.newton$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "M-SGD", tsne = model.msgd$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "C-SGD", tsne = model.csgd$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "R-SGD", tsne = model.rsgd$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "B-SGD", tsne = model.bsgd$tsne, group = cells$Group, batch = cells$Batch)
), by = 5)

print(plt)

if (SAVE) {
  filename = paste("tsne-cell-type-", 100*na.freq, ".pdf", sep = "")
  path = "img"
  width = 5
  height = 2
  zoom = 8
  ggsave(filename = filename, plot = plt, path = "img",
         width = zoom * width, height = zoom * height, units = "cm")
}


plot.tsne.grid(list(
  list(model = "Pearson", tsne = model.pearson$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "Deviance", tsne = model.deviance$tsne, group = cells$Group, batch = cells$Batch),
  # list(model = "NNLM", tsne = model.nnlm$tsne, group = cells$Group, batch = cells$Batch),
  # list(model = "NMF", tsne = model.nmf$tsne, group = cells$Group, batch = cells$Batch),
  # list(model = "CMF", tsne = model.cmf$tsne, group = cells$Group, batch = cells$Batch),
  # list(model = "glmPCA", tsne = model.glmpca$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "NBWaVE", tsne = model.nbwave$tsne, group = cells$Group, batch = cells$Batch),
  # list(model = "SVDreg", tsne = model.svdreg$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "AIRWLS", tsne = model.airwls$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "Newton", tsne = model.newton$tsne, group = cells$Group, batch = cells$Batch),
  # list(model = "M-SGD", tsne = model.msgd$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "C-SGD", tsne = model.csgd$tsne, group = cells$Group, batch = cells$Batch),
  # list(model = "R-SGD", tsne = model.rsgd$tsne, group = cells$Group, batch = cells$Batch),
  list(model = "B-SGD", tsne = model.bsgd$tsne, group = cells$Group, batch = cells$Batch)
), by = 5)

## SILHOUETTE SCORES ----
plot.silhouette = function (sil, main = "Silhouette plot", by = 1) {
  n = nrow(sil)
  sil = cluster::sortSilhouette(sil)
  df = data.frame(index = 1:5000, cluster = as.factor(sil[,1]), width = sil[,3])
  mn = mean(df$width)
  df = df[seq(1, n, by = by),]
  plt = ggplot(data = df, map = aes(x = index, y = width, fill = cluster, color = cluster)) +
    geom_area() + ylim(-1,+1) + theme_minimal() +
    geom_hline(yintercept = mn, col = 2, lty = 2) +
    labs(x = "Cell index", y = "Silhouette width") +
    labs(color = "Cell-type", fill = "Cell-type") +
    labs(title = paste(main, " (", round(mn, 3), ") ", sep = ""))
  return(plt)
}

plot.sil.grid = function (sil, by = 1) {

  models = c()
  means = c()

  df = data.frame(
    model = c(),
    index = c(),
    cluster = c(),
    width = c(),
    mean = c())

  for (l in sil) {
    n = nrow(l$sil)
    s = cluster::sortSilhouette(l$sil)
    m = mean(s[,3])
    dft = data.frame(
      model = rep(l$model, n),
      index = 1:n,
      cluster = s[,1],
      width = s[,3],
      mean = rep(m, n))
    idx = seq(1, n, by = by)
    df = rbind(df, dft[idx, ])
    models = c(models, l$model)
    means = c(means, round(m, 3))
  }

  df$model = factor(df$model, levels = models)
  df$index = as.numeric(df$index)
  df$cluster = as.factor(df$cluster)
  df$width = as.numeric(df$width)
  df$mean = as.numeric(df$mean)

  txt = as.character(round(unique(df$mean), 3))

  plt = ggplot(data = df, map = aes(x = index, fill = cluster, color = cluster)) +
    geom_area(map = aes(y = width)) + geom_line(map = aes(y = mean), col = 1, lty = 2) +
    facet_wrap(vars(model), scales = "fixed") +
    ylim(-1,+1) + theme_grey() +
    # geom_hline(yintercept = mn, col = 2, lty = 2) +
    labs(x = "Cell index", y = "Silhouette width") +
    labs(color = "Cell-type", fill = "Cell-type")


  plt = plt + geom_text(
    data = data.frame(
      model = factor(models, levels = models),
      cluster = as.factor(rep(1, length(sil))),
      x = rep(10, length(sil)),
      y = means+0.15),
    map = aes(x = x, y = y),
    label = means,
    col = 1,
    hjust = .0,
    vjust = .5,
    size = 3
  )

  return(plt)
}

{
sil.pearson = cluster::silhouette(as.numeric(cells$Group), dist(model.pearson$tsne))
sil.deviance = cluster::silhouette(as.numeric(cells$Group), dist(model.deviance$tsne))
sil.glmpca = cluster::silhouette(as.numeric(cells$Group), dist(model.glmpca$tsne))
sil.nbwave = cluster::silhouette(as.numeric(cells$Group), dist(model.nbwave$tsne))
sil.nmf = cluster::silhouette(as.numeric(cells$Group), dist(model.nmf$tsne))
sil.cmf = cluster::silhouette(as.numeric(cells$Group), dist(model.cmf$tsne))
sil.nnlm = cluster::silhouette(as.numeric(cells$Group), dist(model.nnlm$tsne))
sil.svdreg = cluster::silhouette(as.numeric(cells$Group), dist(model.svdreg$tsne))
sil.airwls = cluster::silhouette(as.numeric(cells$Group), dist(model.airwls$tsne))
sil.newton = cluster::silhouette(as.numeric(cells$Group), dist(model.newton$tsne))
sil.msgd = cluster::silhouette(as.numeric(cells$Group), dist(model.msgd$tsne))
sil.csgd = cluster::silhouette(as.numeric(cells$Group), dist(model.csgd$tsne))
sil.rsgd = cluster::silhouette(as.numeric(cells$Group), dist(model.rsgd$tsne))
sil.bsgd = cluster::silhouette(as.numeric(cells$Group), dist(model.bsgd$tsne))
}

plt = plot.sil.grid(list(
  list(model = "Pearson", sil = sil.pearson),
  list(model = "Deviance", sil = sil.deviance),
  # list(model = "glmPCA", sil = sil.glmpca),
  list(model = "NBWaVE", sil = sil.nbwave),
  # list(model = "NMF", sil = sil.nmf),
  # list(model = "CMF", sil = sil.cmf),
  # list(model = "NNLM", sil = sil.nnlm),
  # list(model = "SVDreg", sil = sil.svdreg),
  list(model = "AIRWLS", sil = sil.airwls),
  list(model = "Newton", sil = sil.newton),
  # list(model = "M-SGD", sil = sil.msgd),
  list(model = "C-SGD", sil = sil.csgd),
  # list(model = "R-SGD", sil = sil.rsgd),
  list(model = "B-SGD", sil = sil.bsgd)
), by = 5)

print(plt)

# plt.pearson = plot.silhouette(sil.pearson, "Pearson", by = 5)
# plt.deviance = plot.silhouette(sil.deviance, "Deviance", by = 5)
# plt.glmpca = plot.silhouette(sil.glmpca, "glmPCA", by = 5)
# plt.nbwave = plot.silhouette(sil.nbwave, "NBWaVE", by = 5)
# plt.nmf = plot.silhouette(sil.nmf, "NMF", by = 5)
# plt.cmf = plot.silhouette(sil.cmf, "CMF", by = 5)
# plt.nnlm = plot.silhouette(sil.nnlm, "NNLM", by = 5)
# plt.airwls = plot.silhouette(sil.airwls, "AIRWLS", by = 5)
# plt.newton = plot.silhouette(sil.newton, "Newton", by = 5)
# plt.msgd = plot.silhouette(sil.msgd, "M-SGD", by = 5)
# plt.csgd = plot.silhouette(sil.csgd, "C-SGD", by = 5)
# plt.bsgd = plot.silhouette(sil.bsgd, "B-SGD", by = 5)
#
# plt = ggpubr::ggarrange(
#   plt.pearson, plt.deviance, plt.glmpca, plt.nbwave, plt.nmf, plt.cmf,
#   plt.nnlm, plt.airwls, plt.newton, plt.msgd, plt.csgd, plt.bsgd,
#   nrow = 3, ncol = 4, legend = "bottom", common.legend = TRUE
# )
#
# print(plt)




