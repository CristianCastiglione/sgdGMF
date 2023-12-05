
## WORKSPACE SETUP ----

## Clean the workspace
rm(list = ls())
graphics.off()

## Load our functions
devtools::load_all()
source("tests/utilities.R")

## Load the packages
suppressPackageStartupMessages({
  # Matrix factorization
  library(gllvm)
  library(glmpca)
  library(NewWave)
  library(NMF)
  library(NNLM)
  library(cmfrec)

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

## DATA SIMULATION ----
SETTING = 1
SAVE = FALSE
SHOW = TRUE

n = 1000
m = 100
sim = NULL
if (SETTING == 1) { # ELLIPTIC GROUPS
  params = splatter::newSplatParams()
  params = splatter::setParam(params, "batchCells", c(n/2, n/2))
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
}

## DIMENSION REDUCTION ----
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

df = data.frame(
  u = log.tsne$Y[,1],
  v = log.tsne$Y[,2],
  group = as.factor(cells$Group),
  batch = as.factor(cells$Batch))

levels(df$group) = 1:5
levels(df$batch) = 1:2

plt = ggpubr::ggarrange(
  ggplot(data = df, mapping = aes(x = u, y = v, color = group)) +
    geom_point(alpha = 0.3) + labs(color = "Cell-types") + theme_bw() +
    ggtitle("Cell-type groupping stracture") + labs(x = "tSNE 1", y = "tSNE 2") +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1.0))) +
    theme(legend.position = "right", text = element_text(size = 20),
          # axis.title = element_blank(),
          axis.text = element_blank()),
  ggplot(data = df, mapping = aes(x = u, y = v, color = batch)) +
    geom_point(alpha = 0.3) + labs(color = "Batches   ") + theme_bw() +
    ggtitle("Batch groupping stracture") + labs(x = "tSNE 1", y = "tSNE 2") +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1.0))) +
    theme(legend.position = "right", text = element_text(size = 20),
          # axis.title = element_blank(),
          axis.text = element_blank()),
  nrow = 2, ncol = 1)

print(plt)

if (SAVE) {
  filename = "tsne-sim-data.pdf"
  path = "img"
  width = 3
  height = 4
  zoom = 7
  ggsave(filename = filename, plot = plt, path = "img",
         width = zoom * width, height = zoom * height, units = "cm")
}


## LATENT FACTOR ESTIMATION ----
X = matrix(1, nrow = n, ncol = 1)
X = cbind(1, as.numeric(as.factor(cells$Batch))-1)
Z = matrix(1, nrow = m, ncol = 1)
Y = matrix(NA, nrow = n, ncol = m)
Y[] = t(counts)

ncomp = 5
family = poisson()

# Sparsification and matrix completion
na.freq = 0.3
data = train.test.split(Y, test = na.freq)
train0 = data$ynn
test = data$yna

train0 = test = Y
train0[data$train] = NA
test[data$test] = NA
train = matrix.completion(y = train0, x = X, z = Z,
                          ncomp = ncomp, family = family, niter = 10)

## MODEL ESTIMATION ----
model.pearson = fit.pearson(y = train, x = X, z = Z, ncomp = ncomp, family = family)
model.deviance = fit.deviance(y = train, x = X, z = Z, ncomp = ncomp, family = family)
# model.gllvm = fit.gllvm(y = train, x = X, z = Z, ncomp = ncomp, family = family) # TOO SLOW!!
model.glmpca = fit.glmpca(y = Y, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, tol = 1e-05)
model.cmf = fit.cmf(y = train0, x = X, z = NULL, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 500)
model.nmf = fit.nmf(y = Y, x = NULL, z = NULL, ncomp = ncomp, family = family, verbose = TRUE)
model.nnlm = fit.nnlm(y = train0, x = NULL, z = NULL, ncomp = ncomp, family = family, verbose = TRUE)
model.svdreg = fit.svdreg(y = train0, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 10)
model.airwls = fit.C.airwls(y = Y, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, stepsize = 0.1, tol = 1e-05)
model.newton = fit.C.newton(y = train0, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, stepsize = 0.2)
model.msgd = fit.C.msgd(y = train0, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 100, stepsize = 0.1)
model.csgd = fit.C.csgd(y = train0, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 1000, stepsize = 0.01)
model.bsgd = fit.C.bsgd(y = train0, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 1000, stepsize = 0.01)

## MODEL CHECK ----

# Execution time
{
cat("glmPCA: ", model.glmpca$time[3], "s \n")
cat("NNLM:   ", model.nnlm$time[3], "s \n")
cat("NMF:    ", model.nmf$time[3], "s \n")
cat("CMF:    ", model.cmf$time[3], "s \n")
cat("SVDreg: ", model.svdreg$time[3], "s \n")
cat("AIRWLS: ", model.airwls$time[3], "s \n")
cat("Newton: ", model.newton$time[3], "s \n")
cat("M-SGD:  ", model.msgd$time[3], "s \n")
cat("C-SGD:  ", model.csgd$time[3], "s \n")
cat("B-SGD:  ", model.bsgd$time[3], "s \n")
}

# Convergence check
plot(model.glmpca$dev, type = "b", xlab = "iterations", ylab = "deviance")
plot(model.airwls$dev, type = "b", xlab = "iterations", ylab = "deviance")
plot(model.newton$dev, type = "b", xlab = "iterations", ylab = "deviance")
plot(model.msgd$dev, type = "b", xlab = "iterations", ylab = "deviance")
plot(model.csgd$dev, type = "b", xlab = "iterations", ylab = "deviance")
plot(model.bsgd$dev, type = "b", xlab = "iterations", ylab = "deviance")

## Reconstruction error
full.sample.error = error.matrix(
  Y, model.glmpca, model.nmf, model.nnlm, model.cmf, model.svdreg,
  model.airwls, model.newton, model.msgd, model.csgd, model.bsgd)
print(full.sample.error)

# Out-of-sample error
in.sample.error = error.matrix(
  train0, model.glmpca, model.nmf, model.nnlm, model.cmf, model.svdreg,
  model.airwls, model.newton, model.msgd, model.csgd, model.bsgd)
print(in.sample.error)

# Out-of-sample error
out.sample.error = error.matrix(
  test, model.glmpca, model.nmf, model.nnlm, model.cmf, model.svdreg,
  model.airwls, model.newton, model.msgd, model.csgd, model.bsgd)
print(out.sample.error)

## TSNE PROJECTION ----
tsne.pearson = Rtsne::Rtsne(model.pearson$u, dims = 2, partial_pca = FALSE, num_threads = 7, verbose = TRUE)
tsne.deviance = Rtsne::Rtsne(model.deviance$u, dims = 2, partial_pca = FALSE, num_threads = 7, verbose = TRUE)
tsne.glmpca = Rtsne::Rtsne(model.glmpca$u, dims = 2, partial_pca = FALSE, num_threads = 7, verbose = TRUE)
tsne.nnlm = Rtsne::Rtsne(model.nnlm$u, dims = 2, partial_pca = FALSE, num_threads = 7, verbose = TRUE)
tsne.nmf = Rtsne::Rtsne(model.nmf$u, dims = 2, partial_pca = FALSE, num_threads = 7, verbose = TRUE)
tsne.cmf = Rtsne::Rtsne(model.cmf$u, dims = 2, partial_pca = FALSE, num_threads = 7, verbose = TRUE)
tsne.svdreg = Rtsne::Rtsne(model.svdreg$u, dims = 2, partial_pca = FALSE, num_threads = 7, verbose = TRUE)
tsne.airwls = Rtsne::Rtsne(model.airwls$u, dims = 2, partial_pca = FALSE, num_threads = 7, verbose = TRUE)
tsne.newton = Rtsne::Rtsne(model.newton$u, dims = 2, partial_pca = FALSE, num_threads = 7, verbose = TRUE)
tsne.msgd = Rtsne::Rtsne(model.msgd$u, dims = 2, partial_pca = FALSE, num_threads = 7, verbose = TRUE)
tsne.csgd = Rtsne::Rtsne(model.csgd$u, dims = 2, partial_pca = FALSE, num_threads = 7, verbose = TRUE)
tsne.bsgd = Rtsne::Rtsne(model.bsgd$u, dims = 2, partial_pca = FALSE, num_threads = 7, verbose = TRUE)

df = rbind(
  data.frame(x = tsne.pearson$Y[,1],  y = tsne.pearson$Y[,2],  group = cells$Group, batch = as.factor(cells$Batch), model = rep("Pearson", length = n)),
  data.frame(x = tsne.deviance$Y[,1], y = tsne.deviance$Y[,2], group = cells$Group, batch = as.factor(cells$Batch), model = rep("Deviance", length = n)),
  data.frame(x = tsne.nnlm$Y[,1],     y = tsne.nnlm$Y[,2],     group = cells$Group, batch = as.factor(cells$Batch), model = rep("NNLM", length = n)),
  data.frame(x = tsne.nmf$Y[,1],      y = tsne.nmf$Y[,2],      group = cells$Group, batch = as.factor(cells$Batch), model = rep("NMF", length = n)),
  data.frame(x = tsne.cmf$Y[,1],      y = tsne.cmf$Y[,2],      group = cells$Group, batch = as.factor(cells$Batch), model = rep("CMF", length = n)),
  data.frame(x = tsne.glmpca$Y[,1],   y = tsne.glmpca$Y[,2],   group = cells$Group, batch = as.factor(cells$Batch), model = rep("glmPCA", length = n)),
  # data.frame(x = tsne.svdreg$Y[,1],   y = tsne.svdreg$Y[,2],   group = cells$Group, batch = as.factor(cells$Batch), model = rep("SVDReg", length = n)),
  data.frame(x = tsne.airwls$Y[,1],   y = tsne.airwls$Y[,2],   group = cells$Group, batch = as.factor(cells$Batch), model = rep("AIRWLS", length = n)),
  data.frame(x = tsne.newton$Y[,1],   y = tsne.newton$Y[,2],   group = cells$Group, batch = as.factor(cells$Batch), model = rep("Newton", length = n)),
  # data.frame(x = tsne.msgd$Y[,1],     y = tsne.msgd$Y[,2],     group = cells$Group, batch = as.factor(cells$Batch), model = rep("M-SGD", length = n)),
  data.frame(x = tsne.csgd$Y[,1],     y = tsne.csgd$Y[,2],     group = cells$Group, batch = as.factor(cells$Batch), model = rep("C-SGD", length = n)),
  data.frame(x = tsne.bsgd$Y[,1],     y = tsne.bsgd$Y[,2],     group = cells$Group, batch = as.factor(cells$Batch), model = rep("B-SGD", length = n)))

df$model = factor(df$model, levels = c("Pearson", "Deviance", "NNLM", "NMF", "CMF", "glmPCA", "AIRWLS", "Newton", "C-SGD", "B-SGD"))
levels(df$model) = c("Pearson", "Deviance", "NNLM", "NMF", "CMF", "glmPCA", "AIRWLS", "Newton", "C-SGD", "B-SGD")
levels(df$group) = 1:5
levels(df$batch) = 1:2

# Cell-type clustering
colors = c("#ff7c04", "#387cbc", "#e81c1c", "#50ac4c", "#a04ca4")
plt = ggplot(data = df, mapping = aes(x = x, y = y, color = group)) +
  geom_point(alpha = 0.3) + facet_wrap(vars(model), nrow = 2, ncol = 5) +
  labs(color = "Cell-type") + theme_minimal() +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1.0))) +
  theme(text = element_text(size = 20),
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_blank())

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

# Batch effect
plt = ggplot(data = df, mapping = aes(x = x, y = y, color = batch)) +
  geom_point(alpha = 0.2) + facet_wrap(vars(model), nrow = 2, ncol = 5) +
  labs(color = "Batch    ") + theme_minimal() +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1.0))) +
  theme(text = element_text(size = 20),
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_blank())

print(plt)

if (SAVE) {
  filename = paste("tsne-cell-batch-", 100*na.freq, ".pdf", sep = "")
  path = "img"
  width = 5
  height = 2
  zoom = 8
  ggsave(filename = filename, plot = plt, path = "img",
         width = zoom * width, height = zoom * height, units = "cm")
}


## SILHOUETTE SCORES ----
sil.pearson = cluster::silhouette(as.numeric(cells$Group), dist(tsne.pearson$Y))
sil.deviance = cluster::silhouette(as.numeric(cells$Group), dist(tsne.deviance$Y))
sil.glmpca = cluster::silhouette(as.numeric(cells$Group), dist(tsne.glmpca$Y))
sil.nmf = cluster::silhouette(as.numeric(cells$Group), dist(tsne.nmf$Y))
sil.cmf = cluster::silhouette(as.numeric(cells$Group), dist(tsne.cmf$Y))
sil.nnlm = cluster::silhouette(as.numeric(cells$Group), dist(tsne.nnlm$Y))
sil.airwls = cluster::silhouette(as.numeric(cells$Group), dist(tsne.airwls$Y))
sil.newton = cluster::silhouette(as.numeric(cells$Group), dist(tsne.newton$Y))
sil.csgd = cluster::silhouette(as.numeric(cells$Group), dist(tsne.csgd$Y))
sil.bsgd = cluster::silhouette(as.numeric(cells$Group), dist(tsne.bsgd$Y))
sil.msgd = cluster::silhouette(as.numeric(cells$Group), dist(tsne.msgd$Y))

plt.pearson = factoextra::fviz_silhouette(sil.pearson) + ylim(-1,+1) + ggtitle(paste("Pearson (", round(mean(sil.pearson[,3]), 2), ")", sep = ""))
plt.deviance = factoextra::fviz_silhouette(sil.deviance) + ylim(-1,+1) + ggtitle(paste("Deviance (", round(mean(sil.deviance[,3]), 2), ")", sep = ""))
plt.glmpca = factoextra::fviz_silhouette(sil.glmpca) + ylim(-1,+1) + ggtitle(paste("glmPCA (", round(mean(sil.glmpca[,3]), 2), ")", sep = ""))
plt.nmf = factoextra::fviz_silhouette(sil.nmf) + ylim(-1,+1) + ggtitle(paste("NMF (", round(mean(sil.nmf[,3]), 2), ")", sep = ""))
plt.cmf = factoextra::fviz_silhouette(sil.cmf) + ylim(-1,+1) + ggtitle(paste("CMF (", round(mean(sil.cmf[,3]), 2), ")", sep = ""))
plt.nnlm = factoextra::fviz_silhouette(sil.nnlm) + ylim(-1,+1) + ggtitle(paste("NNLM (", round(mean(sil.nnlm[,3]), 2), ")", sep = ""))
plt.airwls = factoextra::fviz_silhouette(sil.airwls) + ylim(-1,+1) + ggtitle(paste("AIRWLS (", round(mean(sil.airwls[,3]), 2), ")", sep = ""))
plt.newton = factoextra::fviz_silhouette(sil.newton) + ylim(-1,+1) + ggtitle(paste("Newton (", round(mean(sil.newton[,3]), 2), ")", sep = ""))
plt.csgd = factoextra::fviz_silhouette(sil.csgd) + ylim(-1,+1) + ggtitle(paste("C-SGD (", round(mean(sil.csgd[,3]), 2), ")", sep = ""))
plt.bsgd = factoextra::fviz_silhouette(sil.bsgd) + ylim(-1,+1) + ggtitle(paste("B-SGD (", round(mean(sil.bsgd[,3]), 2), ")", sep = ""))
plt.msgd = factoextra::fviz_silhouette(sil.msgd) + ylim(-1,+1) + ggtitle(paste("M-SGD (", round(mean(sil.msgd[,3]), 2), ")", sep = ""))

ggpubr::ggarrange(
  plt.pearson, plt.deviance, plt.glmpca, plt.nmf, plt.cmf, plt.nnlm,
  plt.airwls, plt.newton, plt.csgd, plt.bsgd, plt.msgd,
  nrow = 3, ncol = 4, legend = "bottom", common.legend = TRUE
)


