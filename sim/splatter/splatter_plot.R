
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
SAVE = FALSE
SHOW = TRUE

## Some utility functions
join.path = function (...) paste(..., sep = "/")
join.string = function (...) paste(..., sep = "")

## UTILITY FUNCTIONS ----

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
    filepath = join.path("data", "splatter")
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

merge.data = function (error) {
  H = length(error)
  df = data.frame(
    Model = c(), Time = c(), RSS = c(), Cos = c(), Dev = c(), Sil = c(),
    Set = c(), Iter = c(), Data = c(), Dim = c(), NComp = c())
  for (h in 1:H) {
    dfh = error[[h]]$df
    dfh$Time[!is.finite(dfh$Time)] = NA
    dfh$Dev[!is.finite(dfh$Dev)] = NA
    dfh$RSS[!is.finite(dfh$RSS)] = NA
    dfh$Cos[!is.finite(dfh$Cos)] = NA
    df = rbind(df, dfh)
  }

  return (df)
}

plot.error.list = function (error) {
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

plot.error.data = function (error) {

  df = error[error$Data == "Test", ]
  df = df |>
    dplyr::group_by(Model, Set, Dim, NComp) |>
    dplyr::summarize(
      Time = mean(Time, na.rm = TRUE),
      Dev = mean(Dev, na.rm = TRUE),
      RSS = mean(RSS, na.rm = TRUE),
      Cos = mean(Cos, na.rm = TRUE),
      Sil = mean(Sil, na.rm = TRUE)) |>
    as.data.frame()

  # Remap the labels of Set
  df$Set = as.factor(df$Set)
  levels(df$Set) = c("Branch", "Bubble", "Linear")

  # Remap the labels of Dim
  df$Dim[df$Dim == 10] = 1
  df$Dim[df$Dim == 25] = 2
  df$Dim[df$Dim == 50] = 3
  df$Dim[df$Dim == 75] = 4

  # Chenge the column names
  df = df[, c("Model", "Set", "Dim", "Time", "Dev", "RSS", "Sil")]
  colnames(df) = c("Model", "Setting", "Dim", "Time", "Deviance", "RMSE", "Silhouette")

  # Reorganize the data-frame into a flat representation
  bigdf = data.frame(Model = c(), Setting = c(), Dim = c())
  bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "Dim")]))
  bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "Dim")]))
  bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "Dim")]))
  bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "Dim")]))

  bigdf$Values = c(log10(df$Time), 100*df$Deviance, 100*df$RMSE, df$Silhouette)
  bigdf$Variables = rep(c("Time", "Deviance", "RMSE", "Silhouette"), each = nrow(df))
  bigdf$Variables = as.factor(bigdf$Variables)

  # Plot the grid of summaries
  plt = ggplot(data = bigdf, map = aes(x = Dim, y = Values, color = Model, fill = Model)) +
    # geom_bar(stat="identity", color = "white", position=position_dodge()) +
    geom_line() + geom_point(size = 3) +
    facet_grid(cols = vars(Setting), rows = vars(Variables), scales = "free_y") +
    theme_bw() + theme(legend.position = "bottom") +
    labs(x = "", y = "", fill = "Model", color = "Model") +
    scale_x_continuous(breaks = 1:4, labels = c("10", "25", "50", "75"))

  return (plt)
}

## SIMULATION RESULTS ----

error.list = load.data()
error.data = merge.data(error.list)

plot.error.list(error.list)
plot.error.data(error.data)

## END OF FILE ----













