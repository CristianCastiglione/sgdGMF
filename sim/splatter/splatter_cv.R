
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
  plt = ggplot(data = df, mapping = aes(x = x, y = y, color = group)) +
    geom_point(alpha = 0.5, size = 1) + facet_wrap(vars(model)) +
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
      model = rep(l$model, times = n),
      index = seq(n),
      cluster = as.numeric(s[,1]),
      width = as.numeric(s[,3]),
      mean = rep(m, times = n))
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

## MODEL SELECTION ----
method = "bsgd"
cvopt = list(nfolds = 5)
control = list(maxiter = 1000, size = c(100,20), frequency = 250, rate0 = 0.01)

gmf.fit = sgdGMF::sgdgmf.cv(Y, X, Z, ncomp = 10, family = family,
                            method = method, cvopt = cvopt, control = control)

cvstat = data.frame(
  vars = as.factor(rep(c("BIC", "Deviance"), each = nrow(gmf.fit$cv.stat))),
  vals = c(gmf.fit$cv.stat$bic, gmf.fit$cv.stat$dev),
  fold = c(gmf.fit$cv.stat$fold, gmf.fit$cv.stat$fold),
  ncomp = as.factor(c(gmf.fit$cv.stat$ncomp, gmf.fit$cv.stat$ncomp)))

plt.gof = ggplot(data = cvstat, map = aes(x = ncomp, y = vals, color = ncomp, fill = ncomp)) +
  geom_boxplot(alpha = 0.3) + facet_grid(rows = vars(vars), scale = "free_y") +
  # facet_wrap(vars(vars), nrow = 2, ncol = 1, scale = "free_y") +
  theme(legend.position = "none") + labs(x = "Matrix rank", y = "Goodness of fit") +
  ggtitle("Model selection criterion")

print(plt.gof)

## SEQUENTIAL MODEL FIT ----
gmf.fit.list = list()
for (ncomp in 1:10) {
  cat(" Rank =", ncomp, "\n")
  # Fit the model
  fit = sgdGMF::sgdgmf.fit(Y, X, Z, ncomp = ncomp, family = family,
                           method = method, control = control)

  # Store the estimated model
  gmf.fit.list[[ncomp]] = fit
}

for (ncomp in 1:10) {
  cat(" Rank =", ncomp, "\n")
  fit = gmf.fit.list[[ncomp]]

  # Compute a 2D tSNE embedding
  tsn = scale(Rtsne::Rtsne(as.matrix(fit$U), dims = 2,
                           parallel = TRUE, num_threads = 8)$Y)

  # Compute the silhouette obs-by-obs
  sil = cluster::silhouette(groups, dist(tsn))

  dfsil = as.data.frame(cbind(ncomp, sil[,-2]))
  colnames(dfsil) = c("ncomp", "group", "value")
  dfsil$ncomp = as.factor(dfsil$ncomp)
  dfsil$group = as.factor(dfsil$group)

  # Average the silhouette group-by-group
  avgsil = dfsil %>% group_by(ncomp, group) %>%
    reframe(value = mean(value)) %>% as.data.frame()

  # Store the results
  fit$tsne = tsn
  fit$sil = sil
  fit$avgsil = avgsil
  gmf.fit.list[[ncomp]] = fit
}

## EIGENVALUE ESTIMATION ----
# Vt = sgdGMF::sgdgmf.fit(Y, X, Z, ncomp = 10, family = family,
#                         method = method, control = control)$V
#
# Ct = crossprod(Vt) / nrow(Vt)
# eigval = eigen(Ct)$values
# barplot(eigval[1:10])
#
# Dt = diag(1/sqrt(diag(Ct)))
# Rt = Dt %*% Ct %*% Dt
# eigval = eigen(Rt)$values
# barplot(eigval[1:10])
#
# heatmap(Rt, symm = TRUE, scale = "none", distfun = function(x) as.dist(1 - x),
#         col = hcl.colors(50, palette = "RdBu"))

## SUMMARY PLOTS ----

# Silhouette summary
df.sil = do.call("rbind", lapply(gmf.fit.list, function(x) x$avgsil))
plt.sil = ggplot(data = df.sil, map = aes(x = ncomp, y = value, color = ncomp, fill = ncomp)) +
  geom_boxplot(alpha = 0.3) + labs(x = "Matrix rank", y = "Silhouette") +
  ggtitle("Silhouette path") + theme(legend.position = "none")

print(plt.sil)

# tSNE embeddings
plt.tsne = plot.tsne.grid(list(
  list(model = "d = 2", tsne = gmf.fit.list[[2]]$tsne, group = groups, batch = batches),
  list(model = "d = 3", tsne = gmf.fit.list[[3]]$tsne, group = groups, batch = batches),
  list(model = "d = 4", tsne = gmf.fit.list[[4]]$tsne, group = groups, batch = batches),
  list(model = "d = 5", tsne = gmf.fit.list[[5]]$tsne, group = groups, batch = batches),
  list(model = "d = 6", tsne = gmf.fit.list[[6]]$tsne, group = groups, batch = batches),
  list(model = "d = 10", tsne = gmf.fit.list[[10]]$tsne, group = groups, batch = batches)
), by = 5)

print(plt.tsne)

# Silhouette profile
plt.sil = plot.sil.grid(list(
  list(model = "rank = 2", sil = gmf.fit.list[[2]]$sil),
  list(model = "rank = 3", sil = gmf.fit.list[[3]]$sil),
  list(model = "rank = 4", sil = gmf.fit.list[[4]]$sil),
  list(model = "rank = 5", sil = gmf.fit.list[[5]]$sil),
  list(model = "rank = 6", sil = gmf.fit.list[[6]]$sil),
  list(model = "rank = 10", sil = gmf.fit.list[[10]]$sil)
), by = 5)

print(plt.sil)

# Summary

dfstat = data.frame(
  vars = as.factor(c(rep("BIC", 100), rep("Deviance", 100), rep("Silhouette", 50))),
  vals = c(gmf.fit$cv.stat$bic, gmf.fit$cv.stat$dev, df.sil$value),
  ncomp = as.factor(c(gmf.fit$cv.stat$ncomp, gmf.fit$cv.stat$ncomp, df.sil$ncomp))
)

plt.gof = ggplot(data = dfstat, map = aes(x = ncomp, y = vals, color = ncomp, fill = ncomp)) +
  geom_boxplot(alpha = 0.3) + facet_grid(rows = vars(vars), scale = "free_y") +
  labs(x = "Matrix rank", y = "Goodness of fit") +
  ggtitle("Model selection criterion") +
  theme(legend.position = "none",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 13),
        strip.text = element_text(size = 13),
        axis.title.y = element_blank())

plt.gof

plt = ggpubr::ggarrange(plt.gof + ggtitle("Model selection criteria"),
                        plt.tsne + ggtitle("tSNE projection"),
                        nrow = 1, ncol = 2, widths = c(2,3))
plt

if (SAVE) {
  filename = paste("model_selection.pdf", sep = "")
  path = "img/splatter"
  width = 2
  height = 1
  zoom = 14
  ggsave(filename = filename, plot = plt, path = path,
         width = zoom * width, height = zoom * height, units = "cm")
}
