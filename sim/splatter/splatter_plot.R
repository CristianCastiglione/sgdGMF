
## WORKSPACE SETUP ----

## Clean the workspace
rm(list = ls())
graphics.off()

## Load our functions
devtools::load_all()

## Load the packages
suppressPackageStartupMessages({
  library(gllvm)
  library(glmpca)

  library(scry)
  library(splatter)
  library(scater)

  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
})

## Global variables
SAVE = TRUE
SHOW = FALSE

## Some utility functions
join.path = function (...) paste(..., sep = "/")
join.string = function (...) paste(..., sep = "")

## SIMULATION RESULTS ----

load.data = function () {

  models = c("CMF", "NMF", "NNLM", "glmPCA", "NBWaVE", "AIRWLS", "Newton", "C-SGD", "B-SGD")
  setting = c("bubble", "linear", "branch")
  # nrows = c(1000, 2500, 5000, 7500)
  # ncols = c(100, 250, 500, 750)
  rcdims = c(10, 25, 50, 75)
  ncomps = c(5)
  niters = c(25)

  h = 1
  files = list()
  for (st in setting) {
    for (dm in rcdims) {
      for (nc in ncomps) {
        for (ni in niters) {
          strh = join.string("_s", st, "_n", 100*dm, "_m", 10*dm, "_d", nc, "_i", ni)
          idh = list(setting = st, rcdim = dm, ncomp = nc, niter = ni, str = strh)
          files[[h]] = idh
          h = h + 1
        }
      }
    }
  }

  H = h - 1
  error = list()
  for (h in 1:H) {
    fileh = files[[h]]
    filepath = join.path("sim", "splatter")
    filename = join.string("summary", fileh$str, ".csv")
    fileid = join.path(filepath, filename)

    errorh = list()
    errorh$setting = fileh$setting
    errorh$rcdim = fileh$rcdim
    errorh$ncomp = fileh$ncomp
    errorh$path = filepath
    errorh$file = filename
    errorh$str = fileh$str
    errorh$df = read.table(file = fileid, header = TRUE, sep = ";")

    errorh$df$Model = factor(errorh$df$Model, levels = models)
    errorh$df$Data = factor(errorh$df$Set, levels = c("Full", "Train", "Test"))
    errorh$df$Set = factor(errorh$setting, levels = setting)
    errorh$df$Iter = as.numeric(errorh$df$Iter)
    errorh$df$Time = as.numeric(errorh$df$Time)
    errorh$df$RSS = as.numeric(errorh$df$RSS)
    errorh$df$Cos = as.numeric(errorh$df$Cos)
    errorh$df$Dev = as.numeric(errorh$df$Dev)
    errorh$df$Sil = as.numeric(errorh$df$Sil)
    errorh$df$Dim = as.numeric(errorh$rcdim)
    errorh$df$NComp = as.numeric(errorh$ncomp)

    error[[h]] = errorh
  }

  return (error)
}

plot.summary = function (error) {
  sim.dim.str = function (n, m, d)
    paste0("n = ", n, ", m = ", m, ", d = ", d)

  H = length(error)
  for (h in 1:H) {
    errorh = error[[h]]

    dm = errorh$rcdim
    nc = errorh$ncomp
    df = errorh$df[, c("Model", "Data", "Set", "Time", "Dev", "RSS", "Sil")]

    df = data.frame(
      Model = rep(df$Model, times = 4),
      Data = rep(df$Data, times = 4),
      # Set = rep(df$Set, times = 4),
      Val = c(log10(df$Time), 100*df$Dev, 100*df$RSS, df$Sil),
      Var = rep(c("Exe-Time", "Res-Dev", "RMSE", "Silhouette"), each = nrow(df)))

    df$Var = factor(df$Var, levels = c("Exe-Time", "Res-Dev", "RMSE", "Silhouette"))

    plt = ggplot(data = df, map = aes(x = Model, y = Val, color = Data, fill = Data)) +
      geom_boxplot(alpha = 0.5) + facet_grid(rows = vars(Var), scales = "free_y") +
      ggtitle(sim.dim.str(100*dm, 10*dm, nc)) +
      theme_bw() + theme(axis.title = element_blank())

    if (SHOW) print(plt)

    if (SAVE) {
      filepath = join.path("img", "splatter")
      filename = join.string("summary", errorh$str, ".pdf")
      zoom = 8
      width = 1
      height = 0.8
      ggplot2::ggsave(
        file = filename, path = filepath, plot = plt,
        width = zoom * width, height = zoom * height)
    }
  }
}

error = load.data()
plot.summary(error)




models = c("CMF", "NMF", "NNLM", "glmPCA", "NBWaVE", "AIRWLS", "Newton", "C-SGD", "B-SGD")
setting = c("bubble", "linear", "branch")
nrows = c(1000, 2500, 5000, 7500)
ncols = c(100, 250, 500, 750)
ncomps = c(5)
niters = c(25)

h = 1
fileid = list()
for (st in setting) {
  for (nc in ncols) {
    for (nf in ncomps) {
      for (ni in niters) {
        strh = join.string("_s", st, "_n", 10*nc, "_m", nc, "_d", nf, "_i", ni)
        idh = list(setting = st, nrows = 10*nc, ncols = nc, ncomp = nf, niter = ni, str = strh)
        fileid[[h]] = idh
        h = h + 1
      }
    }
  }
}

H = h - 1
error = list()
for (h in 1:H) {
  filepath = join.path("sim", "splatter")
  filename = join.string("summary", fileid[[h]]$str, ".csv")
  errorh = list()
  errorh$setting = fileid[[h]]$setting
  errorh$nrows = fileid[[h]]$nrows
  errorh$mcols = fileid[[h]]$ncols
  errorh$ncomp = fileid[[h]]$ncomp
  errorh$path = filepath
  errorh$file = filename
  errorh$df = read.table(file = join.path(filepath, filename), header = TRUE, sep = ";")

  errorh$df$Model = factor(errorh$df$Model, levels = models)
  errorh$df$Set = factor(errorh$df$Set, levels = c("Full", "Train", "Test"))
  errorh$df$Time = as.numeric(errorh$df$Time)
  errorh$df$RSS = as.numeric(errorh$df$RSS)
  errorh$df$Cos = as.numeric(errorh$df$Cos)
  errorh$df$Dev = as.numeric(errorh$df$Dev)
  errorh$df$Sil = as.numeric(errorh$df$Sil)
  errorh$df$Iter = as.numeric(errorh$df$Iter)

  error[[h]] = errorh
}

for (h in 1:H) {
  nr = error[[h]]$nrows
  nc = error[[h]]$mcols
  nf = error[[h]]$ncomp
  df = error[[h]]$df
  sim.dim.str = function (n, m, d) paste0("(n = ", n, ", m = ", m, ", d = ", d, ")")

  ggplt = ggplot(data = df, mapping = aes(x = Model, color = Set, fill = Set)) +
    labs(fill = "", color = "") + theme_bw() + theme(axis.title = element_blank())
  plt = ggpubr::ggarrange(
    ggplt + geom_boxplot(mapping = aes(y = log10(Time)), alpha = 0.5) +
      ggtitle(paste("Log-execution time", sim.dim.str(nr, nc, nf))),
    ggplt + geom_boxplot(mapping = aes(y = Dev), alpha = 0.5) +
      ggtitle(paste("Residual deviance", sim.dim.str(nr, nc, nf))),
    ggplt + geom_boxplot(mapping = aes(y = RSS), alpha = 0.5) +
      ggtitle(paste("Residual sum of squares", sim.dim.str(nr, nc, nf))),
    ggplt + geom_boxplot(mapping = aes(y = Sil), alpha = 0.5) +
      ggtitle(paste("Average silhouette score", sim.dim.str(nr, nc, nf))),
    nrow = 4, ncol = 1, legend = "bottom", common.legend = TRUE)

  if (SHOW) print(plt)

  if (SAVE) {
    filepath = join.path("img", "splatter")
    filename = join.string("summary", fileid[[h]]$str, ".pdf")
    zoom = 7
    width = 1
    height = 1.2
    ggplot2::ggsave(
      file = filename, path = filepath, plot = plt,
      width = zoom * width, height = zoom * height)
  }
}


h = 1
{
  nr = error[[h]]$nrows
  nc = error[[h]]$mcols
  nf = error[[h]]$ncomp
  df = error[[h]]$df

  df$Time[!is.finite(df$Time)] = NA
  df$Dev[!is.finite(df$Dev)] = NA
  df$RSS[!is.finite(df$RSS)] = NA
  df$Cos[!is.finite(df$Cos)] = NA

  df = df[df$Set == "Test",]
  df = df |>
    dplyr::group_by(Model) |>
    dplyr::summarize(
      Time = mean(Time, na.rm = TRUE),
      Dev = mean(Dev, na.rm = TRUE),
      RSS = mean(RSS, na.rm = TRUE),
      Cos = mean(Cos, na.rm = TRUE))


  sim.dim.str = function (n, m, d) paste0("(n = ", n, ", m = ", m, ", d = ", d, ")")

  plt = ggplot(data = df, mapping = aes(x = log10(Time), y = Dev, color = Model, fill = Model)) +
    geom_point(size = 5) + theme_bw() +
    labs(x = "Log-execution time", y = "Residual deviance", fill = "Model", color = "Model")

  print(plt)

  h = h+1
}

dat = data.frame(
  Model = c(), Time = c(), RSS = c(),
  Cos = c(), Dev = c(), Sil = c(),
  Set = c(), Iter = c(), Data = c(),
  n = c(), m = c(), d = c())
for (h in 1:H) {
  nr = error[[h]]$nrows
  nc = error[[h]]$mcols
  nf = error[[h]]$ncomp
  df = error[[h]]$df

  df$Data = rep(error[[h]]$setting, times = nrow(df))
  df$n = rep(nr, times = nrow(df))
  df$m = rep(nc, times = nrow(df))
  df$d = rep(nf, times = nrow(df))

  df$Time[!is.finite(df$Time)] = NA
  df$Dev[!is.finite(df$Dev)] = NA
  df$RSS[!is.finite(df$RSS)] = NA
  df$Cos[!is.finite(df$Cos)] = NA

  dat = rbind(dat, df)
}

str(dat)

df = dat
df = df[df$Set == "Test",]
df = df |>
  dplyr::group_by(Model, Data, n, m, d) |>
  dplyr::summarize(
    Time = mean(Time, na.rm = TRUE),
    Dev = mean(Dev, na.rm = TRUE),
    RSS = mean(sqrt(RSS), na.rm = TRUE),
    Cos = mean(Cos, na.rm = TRUE),
    Sil = mean(Sil, na.rm = TRUE))

df$Data = as.factor(df$Data)
levels(df$Data) = c("Branch", "Bubble", "Linear")



library(ggrepel)
plt = ggplot(data = df, mapping = aes(x = log10(Time), y = Dev, color = Model, fill = Model)) +
  geom_line() + geom_point(size = 3) + facet_grid(cols = vars(Data)) +
  geom_text_repel(map = aes(label = as.character(n/100)), color = "black", size = 2.5) +
  # geom_label_repel(map = aes(label = as.character(n)), color = "white") +
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = "Log-execution time", y = "Residual deviance", fill = "Model", color = "Model")

plt

plt.time = ggplot(data = df, mapping = aes(x = n, y = log10(Time), color = Model, fill = Model)) +
  geom_line() + geom_point(size = 3) + facet_grid(cols = vars(Data)) +
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = "Matrix dimensions", y = "Log-execution time", fill = "Model", color = "Model")

plt.dev = ggplot(data = df, mapping = aes(x = n, y = Dev, color = Model, fill = Model)) +
  geom_line() + geom_point(size = 3) + facet_grid(cols = vars(Data)) +
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = "Matrix dimensions", y = "Residual deviance", fill = "Model", color = "Model")

plt.rss = ggplot(data = df, mapping = aes(x = n, y = RSS, color = Model, fill = Model)) +
  geom_line() + geom_point(size = 3) + facet_grid(cols = vars(Data)) +
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = "Matrix dimensions", y = "Residual sum of squared", fill = "Model", color = "Model")

plt.sil = ggplot(data = df, mapping = aes(x = n, y = Sil, color = Model, fill = Model)) +
  geom_line() + geom_point(size = 3) + facet_grid(cols = vars(Data)) +
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = "Matrix dimensions", y = "Average silhouette", fill = "Model", color = "Model")

plt = ggpubr::ggarrange(plt.time, plt.dev, plt.rss, plt.sil, nrow = 4, ncol = 1,
                        common.legend = TRUE, legend = "bottom")

print(plt)


df = df[,c("Model", "Data", "n", "Time", "Dev", "RSS", "Sil")]
colnames(df) = c("Model", "Setting", "Dimensions", "Time", "Deviance", "RMSE", "Silhouette")


bigdf = data.frame(Model = c(), Setting = c(), Dimensions = c())
bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "Dimensions")]))
bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "Dimensions")]))
bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "Dimensions")]))
bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "Dimensions")]))

bigdf$Values = c(log10(df$Time), 100*df$Deviance, 100*df$RMSE, df$Silhouette)

bigdf$Variables = rep(c("Time", "Deviance", "RMSE", "Silhouette"), each = nrow(df))
bigdf$Variables = as.factor(bigdf$Variables)

plt = ggplot(data = bigdf, map = aes(x = Dimensions, y = Values, color = Model, fill = Model)) +
  # geom_bar(stat="identity", color = "white", position=position_dodge()) +
  geom_line() + geom_point(size = 3) +
  facet_grid(cols = vars(Setting), rows = vars(Variables), scales = "free_y") +
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = "", y = "", fill = "Model", color = "Model")

print(plt)


