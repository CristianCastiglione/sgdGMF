
## WORKSPACE SETUP ----

## Clean the workspace
rm(list = ls())
graphics.off()

## Load the package
devtools::load_all()

## Load the utility functions
source("sim/utilities.R")

## GLOBAL VARIABLES ----
SETTING = 1
SAVE = TRUE
SHOW = TRUE

colors = c("#F8766D", "#619CFF", "#00BA38", "#ff9933", "#d97ff2")

## PLOTTING FUNCTIONS ----
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

  # colors = c("#ff7c04", "#387cbc", "#e81c1c", "#50ac4c", "#a04ca4")
  plt = ggplot(data = df, mapping = aes(x = x, y = y, color = group, pch = batch)) +
    geom_point(alpha = 0.5, size = 2) + facet_wrap(vars(model)) +
    labs(color = "Cell-type", pch = "Batch") + theme_gray() +
    # scale_colour_manual(values = colors) +
    # scale_fill_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1.0))) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 13),
          strip.text = element_text(size = 13),
          axis.title = element_blank())

#    theme(text = element_text(size = 20),
#          legend.position = "bottom",
#          legend.text = element_text(size = 10),
#          legend.title = element_text(size = 15),
#          axis.title = element_blank(),
#          axis.text = element_blank())

  if (SAVE) {
    filename = paste("example_sbubble_n5000_m500_d5_i25.pdf", sep = "")
    path = "img/splatter"
    width = 1
    height = 1
    zoom = 8
    ggsave(filename = filename, plot = plt, path = path,
           width = zoom * width, height = zoom * height, units = "cm")
  }

  return (plt)
}

plot.sil.score = function (sil, main = "Silhouette plot", by = 1) {
  n = nrow(sil)
  sil = cluster::sortSilhouette(sil)
  df = data.frame(index = 1:5000, cluster = as.factor(sil[,1]), width = sil[,3])
  mn = mean(df$width)
  df = df[seq(1, n, by = by),]
  plt = ggplot(data = df, map = aes(x = index, y = width, fill = cluster, color = cluster)) +
    geom_area() + ylim(-1,+1) + theme_minimal() +
    geom_hline(yintercept = mn, col = 2, lty = 2) +
    # scale_colour_manual(values = colors) +
    # scale_fill_manual(values = colors) +
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
    # scale_colour_manual(values = colors) +
    # scale_fill_manual(values = colors) +
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

plot.stat.summary = function (s = "bubble", n = 5000, m = 500, d = 5, i = 25) {

  fileid = join.string("_s", s, "_n", n, "_m", m, "_d", d, "_i", i)
  filepath = join.path("data", "splatter")
  filename = join.string("summary", fileid, ".csv")

  models = c("CMF", "NMF", "NNLM", "glmPCA", "NBWaVE", "AIRWLS", "Newton", "C-SGD", "B-SGD")
  variables = c("Time", "Deviance", "Error", "Silhouette")

  df = read.table(file = join.path(filepath, filename), header = TRUE, sep = ";")
  df = df[df$Set == "Test", ]
  df = df[, c("Model", "Time", "Dev", "RSS", "Sil")]

  df = data.frame(
    Model = factor(rep(df$Model, times = 4), levels = models),
    Val = c(log10(df$Time), 100 * df$Dev, 100 * df$RSS, df$Sil),
    Var = factor(rep(variables, each = nrow(df)), levels = variables))

  levels(df$Model) = c("CMF", "NMF", "NNLM", "glmPCA", "NBWaVE", "AIRWLS", "Newton", "C-SQN", "B-SQN")

  plt = ggplot(data = df, map = aes(x = Model, y = Val, color = Model, fill = Model)) +
    geom_boxplot(alpha = 0.5) + facet_grid(rows = vars(Var), scales = "free_y") +
    ggtitle(join.string("n = ", 5000, ", m = ", 500, ", d = ", 5)) +
    theme_gray() +
    theme(axis.title = element_blank(),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 12.5, angle = 45, hjust = 1),
          strip.text = element_text(size = 13),
          plot.title = element_text(size = 13),
          legend.position = "none",
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 10))

  return (plt)
}

## LOAD DATA ----
filename = "example_sbubble_n5000_m500_d5_i25.RData"
filepath = join.path("data", "splatter", filename)
load(file = filepath)

## DATA EXTRACTION ----
logcounts = as.data.frame(logcounts(sim))
counts = as.data.frame(counts(sim))
cells = as.data.frame(colData(sim))
genes = as.data.frame(rowData(sim))
meta = metadata(sim)
groups = as.numeric(as.factor(cells$Group))
batches = as.numeric(as.factor(cells$Batch))

# PCA embedding
pca = RSpectra::svds(scale(t(as.matrix(logcounts))), k = 10)$u

# t-SNE embedding
tsne = Rtsne::Rtsne(pca, dims = 2, verbose = TRUE, num_threads = 8)$Y
# tsne = Rtsne::Rtsne(scale(t(as.matrix(logcounts))), dims = 2,
#                     verbose = TRUE, num_threads = 8)$Y

df = data.frame(
  x = c(minmax(pca[,1]), minmax(tsne[,1])),
  y = c(minmax(pca[,2]), minmax(tsne[,2])),
  group = as.factor(rep(groups, 2)),
  batch = as.factor(rep(batches, 2)),
  embedding = as.factor(rep(c("PCA", "tSNE"), each = n))
)


ggplot(data = df, map = aes(x = x, y = y, color = group, pch = batch)) +
  geom_point(alpha = 0.9) + facet_grid(cols = vars(embedding)) +
  # scale_colour_manual(values = colors) + theme_grey() +
  labs(x = "PC1", y = "PC2", color = "Cell-type", pch = "Batch")

# par(mfrow = c(1,2))
# plot(pca, col = groups, pch = batches, xlab = "PC1", ylab = "PC2", main = "PCA: cell-types")
# plot(tsne, col = groups, pch = batches, xlab = "PC1", ylab = "PC2", main = "tSNE: cell-types")
# par(mfrow = c(1,1))

## TRAIN-TEST SPLIT ----
X = model.matrix(~ Batch, data = cells)
Z = matrix(1, nrow = m, ncol = 1)
Y = matrix(NA, nrow = n, ncol = m)
Y[] = t(counts)

test = Y
train = Y
ctrain = Y

test[data$test] = NA
train[data$train] = NA
ctrain = naive.completion(train)

# nullres = null.residuals(Y, X, Z, family)
# pca = RSpectra::svds(scale(as.matrix(nullres)), k = 2)$u
# tsne = Rtsne::Rtsne(scale(as.matrix(nullres)), dims = 2,
#                     verbose = TRUE, num_threads = 8)$Y
# #
# par(mfrow = c(1,2))
# plot(pca, col = groups, pch = batches, xlab = "PC1", ylab = "PC2", main = "PCA: cell-types")
# plot(tsne, col = groups, pch = batches, xlab = "PC1", ylab = "PC2", main = "tSNE: cell-types")
# par(mfrow = c(1,1))


## MODEL FIT ----
model.pearson = fit.pearson(y = ctrain, x = NULL, z = NULL, ncomp = ncomp, family = family)
# model.deviance = fit.deviance(y = ctrain, x = X, z = Z, ncomp = ncomp, family = family)
# model.glmpca = fit.glmpca(y = ctrain, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, tol = 1e-05)
model.nbwave = fit.nbwave(y = ctrain, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, tol = 1e-04)
# model.cmf = fit.cmf(y = train, x = X, z = NULL, ncomp = ncomp, family = family, verbose = FALSE, maxiter = 500)
# model.nmf = fit.nmf(y = ctrain, x = NULL, z = NULL, ncomp = ncomp, family = family, verbose = TRUE)
# model.nnlm = fit.nnlm(y = train, x = NULL, z = NULL, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 2000)
# model.airwls = fit.C.airwls(y = train, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, stepsize = 0.9)
# model.newton = fit.C.newton(y = train, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, stepsize = 0.2)
model.csgd = fit.C.csgd(y = train, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 2000, stepsize = 0.01)
model.bsgd = fit.C.bsgd(y = train, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 2000, stepsize = 0.01)

## TSNE PROJECTION ----
plt.tsne = plot.tsne.grid(list(
  # list(model = "Log-Counts", tsne = tsne, group = groups, batch = batches),
  list(model = "Pearson", tsne = model.pearson$tsne, group = groups, batch = batches),
  # list(model = "Deviance", tsne = model.deviance$tsne, group = groups, batch = batches),
  # list(model = "NNLM", tsne = model.nnlm$tsne, group = groups, batch = batches),
  # list(model = "NMF", tsne = model.nmf$tsne, group = groups, batch = batches),
  # list(model = "CMF", tsne = model.cmf$tsne, group = groups, batch = batches),
  # list(model = "glmPCA", tsne = model.glmpca$tsne, group = groups, batch = batches),
  list(model = "NB-WaVE", tsne = model.nbwave$tsne, group = groups, batch = batches),
  # list(model = "SVDreg", tsne = model.svdreg$tsne, group = groups, batch = batches),
  # list(model = "AIRWLS", tsne = model.airwls$tsne, group = groups, batch = batches),
  # list(model = "Newton", tsne = model.newton$tsne, group = groups, batch = batches),
  # list(model = "M-SGD", tsne = model.msgd$tsne, group = groups, batch = batches),
  list(model = "C-SQN", tsne = model.csgd$tsne, group = groups, batch = batches),
  # list(model = "R-SGD", tsne = model.rsgd$tsne, group = groups, batch = batches),
  list(model = "B-SQN", tsne = model.bsgd$tsne, group = groups, batch = batches)
), by = 5)

print(plt.tsne)



## SILHOUETTE SCORES ----
{
  sil.pearson = cluster::silhouette(groups, dist(model.pearson$tsne))
  # sil.deviance = cluster::silhouette(groups, dist(model.deviance$tsne))
  # sil.glmpca = cluster::silhouette(groups, dist(model.glmpca$tsne))
  sil.nbwave = cluster::silhouette(groups, dist(model.nbwave$tsne))
  # sil.nmf = cluster::silhouette(groups, dist(model.nmf$tsne))
  # sil.cmf = cluster::silhouette(groups, dist(model.cmf$tsne))
  # sil.nnlm = cluster::silhouette(groups, dist(model.nnlm$tsne))
  # sil.svdreg = cluster::silhouette(groups, dist(model.svdreg$tsne))
  # sil.airwls = cluster::silhouette(groups, dist(model.airwls$tsne))
  # sil.newton = cluster::silhouette(groups, dist(model.newton$tsne))
  # sil.msgd = cluster::silhouette(groups, dist(model.msgd$tsne))
  sil.csgd = cluster::silhouette(groups, dist(model.csgd$tsne))
  # sil.rsgd = cluster::silhouette(groups, dist(model.rsgd$tsne))
  sil.bsgd = cluster::silhouette(groups, dist(model.bsgd$tsne))
}

plt.sil = plot.sil.grid(list(
  list(model = "Pearson", sil = sil.pearson),
  # list(model = "Deviance", sil = sil.deviance),
  # list(model = "glmPCA", sil = sil.glmpca),
  list(model = "NBWaVE", sil = sil.nbwave),
  # list(model = "NMF", sil = sil.nmf),
  # list(model = "CMF", sil = sil.cmf),
  # list(model = "NNLM", sil = sil.nnlm),
  # list(model = "SVDreg", sil = sil.svdreg),
  # list(model = "AIRWLS", sil = sil.airwls),
  # list(model = "Newton", sil = sil.newton),
  # list(model = "M-SGD", sil = sil.msgd),
  list(model = "C-SGD", sil = sil.csgd),
  # list(model = "R-SGD", sil = sil.rsgd),
  list(model = "B-SGD", sil = sil.bsgd)
), by = 5)

print(plt.sil)

## LOAD SUMMARY STATS ----

plt.stat = plot.stat.summary()

plt = ggpubr::ggarrange(
  plt.stat + ggtitle("Summary statistics"),
  plt.tsne + ggtitle("tSNE projection"),
  nrow = 1, ncol = 2, widths = c(1,2))

if (SHOW) {
  print(plt)
}

if (SAVE) {
  filename = paste("example_sbubble_n5000_m500_d5_i25.pdf", sep = "")
  path = "img/splatter"
  width = 2
  height = 1.25
  zoom = 14
  ggsave(filename = filename, plot = plt, path = path,
         width = zoom * width, height = zoom * height, units = "cm")
}

## END OF FILE ----



# model selection via cross-validation
model.bsgd = sgdGMF::sgdgmf.cv(
  Y = Y, X = X, Z = Z, ncomp = 10, family = family, method = "bsgd", cvopt = list(nfolds = 4),
  control = list(size = c(100,20), frequency = 500, stepsize = 0.01, rate1 = 0.1, rate2 = 0.01))

model.bsgd$cv.stat |>
  ggplot(map = aes(x = as.factor(ncomp), y = bic, fill = as.factor(ncomp))) +
  geom_boxplot() + theme(legend.position = "none") +
  labs(x = "Matrix rank", y = "Deviance", title = "Out-of-sample deviance path")

library(dplyr)
df = model.bsgd$cv.stat |> group_by(ncomp) |> summarize(
  aic.min = min(aic), bic.min = min(bic), dev.min = min(dev),
  aic.max = max(aic), bic.max = max(bic), dev.max = max(dev),
  aic.sd = sd(aic), bic.sd = sd(bic), dev.sd = sd(dev),
  aic = mean(aic), bic = mean(bic), dev = mean(dev)) |>
  as.data.frame()

ggplot(data = df, map = aes(x = ncomp, y = dev, ymin = dev.min, ymax = dev.max)) +
  geom_line() + geom_point() + theme(legend.position = "none") +
  geom_ribbon(alpha = 0.3) +
  labs(x = "Matrix rank", y = "Deviance", title = "Out-of-sample deviance path")





