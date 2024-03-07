
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
  library(latex2exp)
})

## Global variables
SAVE = TRUE
SHOW = TRUE

## Some utility functions
join.path = function (...) paste(..., sep = "/")
join.string = function (...) paste(..., sep = "")

## UTILITY FUNCTIONS ----

load.flat.data = function () {

  models = c("CMF", "NMF", "NNLM", "glmPCA", "NBWaVE", "AIRWLS", "Newton", "C-SGD", "B-SGD")
  setting = c("bubble", "linear", "branch")
  rcdims = c(10, 25, 50, 75, 100)
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

load.deep.data = function () {

  models = c("CMF", "NMF", "NNLM", "glmPCA", "NBWaVE", "AIRWLS", "Newton", "C-SGD", "B-SGD")
  setting = c("bubble", "linear", "branch")
  rcdims = c(50)
  ncomps = c(5, 10, 15, 20, 25)
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

merge.list = function (error) {
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

plot.list = function (error) {
  sim.dim.str = function (n, m, d)
    paste0("n = ", n, ", m = ", m, ", d = ", d)

  H = length(error)
  for (h in 1:H) {
    errorh = error[[h]]

    dm = errorh$rcdim
    nc = errorh$ncomp
    df = errorh$df[, c("Model", "Data", "Set", "Time", "Dev", "RSS", "Sil")]

    df = df[df$Data == "Test", ]

    labels = c("Time", "Deviance", "Error", "Silhouette")

    df = data.frame(
      Model = rep(df$Model, times = 4),
      Data = rep(df$Data, times = 4),
      # Set = rep(df$Set, times = 4),
      Val = c(log10(df$Time), 100*df$Dev, 100*df$RSS, df$Sil),
      Var = rep(labels, each = nrow(df)))

    df$Var = factor(df$Var, levels = labels)

    plt = ggplot(data = df, map = aes(x = Model, y = Val, color = Model, fill = Model)) +
      geom_boxplot(alpha = 0.5) + facet_grid(rows = vars(Var), scales = "free_y") +
      ggtitle(sim.dim.str(100*dm, 10*dm, nc)) +
      theme_gray() +
      theme(axis.title = element_blank(),
            axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
            strip.text = element_text(size = 15),
            legend.position = "none",
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 12))
      # theme(text = element_text(size = 20),
      #       legend.position = "bottom",
      #       legend.text = element_text(size = 10),
      #       legend.title = element_text(size = 15),
      #       axis.title = element_blank(),
      #       axis.text.x = element_text(angle = 45))

    if (SHOW) print(plt)

    if (SAVE) {
      filepath = join.path("img", "splatter")
      filename = join.string("summary", errorh$str, ".pdf")
      zoom = 8
      width = .6
      height = 1
      ggplot2::ggsave(
        file = filename, path = filepath, plot = plt,
        width = zoom * width, height = zoom * height)
    }
  }
}

plot.flat.data = function (error, setting = c("Bubble", "Linear", "Branch")) {

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
  df$Dim[df$Dim == 100] = 5

  # Chenge the column names
  df = df[, c("Model", "Set", "Dim", "Time", "Dev", "RSS", "Sil")]
  colnames(df) = c("Model", "Setting", "Dim", "Time", "Deviance", "RMSE", "Silhouette")

  df = df[df$Setting %in% c(setting), ]

  # Reorganize the data-frame into a flat representation
  bigdf = data.frame(Model = c(), Setting = c(), Dim = c())
  bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "Dim")]))
  bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "Dim")]))
  bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "Dim")]))
  bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "Dim")]))

  bigdf$Values = c(df$Time, 100*df$Deviance, 100*df$RMSE, df$Silhouette)
  bigdf$Variables = rep(c("Time", "Deviance", "RMSE", "Silhouette"), each = nrow(df))
  bigdf$Variables = as.factor(bigdf$Variables)

  # Plot the grid of summaries
  plt = ggplot(data = bigdf, map = aes(x = Dim, y = Values, color = Model, fill = Model)) +
    # geom_bar(stat="identity", color = "white", position=position_dodge()) +
    geom_line() + geom_point(shape = 21, size = 3, colour = "white") +
    facet_grid(cols = vars(Setting), rows = vars(Variables), scales = "free_y") +
    # theme_bw() +
    theme(legend.position = "bottom") +
    labs(x = "", y = "", fill = "Model", color = "Model") +
    scale_x_continuous(breaks = 1:5, labels = c("10", "25", "50", "75", "100"))

  if (SHOW) print(plt)

  if (SAVE) {
    filepath = join.path("img", "splatter")
    filename = join.string("summary_flat_sim.pdf")
    zoom = 3
    width = 1
    height = .8
    ggplot2::ggsave(
      file = filename, path = filepath, plot = plt,
      width = zoom * length(setting) * width, height = zoom * 4 * height)
  }

  return (plt)
}

plot.deep.data = function (error, setting = c("Bubble", "Linear", "Branch")) {

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

  df$Dev[which(df$Dev > 0.100)] = NA
  df$RSS[which(df$RSS > 0.042)] = NA

  # Chenge the column names
  df = df[, c("Model", "Set", "NComp", "Time", "Dev", "RSS", "Sil")]
  colnames(df) = c("Model", "Setting", "NComp", "Time", "Deviance", "RMSE", "Silhouette")

  df = df[df$Setting %in% c(setting), ]

  # Reorganize the data-frame into a flat representation
  bigdf = data.frame(Model = c(), Setting = c(), NComp = c())
  bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "NComp")]))
  bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "NComp")]))
  bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "NComp")]))
  bigdf = rbind(bigdf, as.data.frame(df[,c("Model", "Setting", "NComp")]))

  bigdf$Values = c(df$Time, 100*df$Deviance, 100*df$RMSE, df$Silhouette)
  bigdf$Variables = rep(c("Time", "Deviance", "RMSE", "Silhouette"), each = nrow(df))
  bigdf$Variables = as.factor(bigdf$Variables)

  # Plot the grid of summaries
  plt = ggplot(data = bigdf, map = aes(x = NComp, y = Values, color = Model, fill = Model)) +
    # geom_bar(stat="identity", color = "white", position=position_dodge()) +
    geom_line() + geom_point(shape = 21, size = 3, colour = "white") +
    facet_grid(cols = vars(Setting), rows = vars(Variables), scales = "free_y") +
    # theme_bw() +
    theme(legend.position = "bottom") +
    labs(x = "", y = "", fill = "Model", color = "Model")

  if (SHOW) print(plt)

  if (SAVE) {
    filepath = join.path("img", "splatter")
    filename = join.string("summary_deep_sim.pdf")
    zoom = 3
    width = 1
    height = .8
    ggplot2::ggsave(
      file = filename, path = filepath, plot = plt,
      width = zoom * length(setting) * width, height = zoom * 4 * height)
  }

  return (plt)
}

plot.contrast = function (flat, deep, setting = "Bubble") {

  flat = flat[flat$Data == "Test", ]
  flat = flat[flat$Set == setting, ]
  flat = flat |>
    dplyr::group_by(Model, Set, Dim, NComp) |>
    dplyr::summarize(
      Time = mean(Time, na.rm = TRUE),
      Dev = mean(Dev, na.rm = TRUE),
      RSS = mean(RSS, na.rm = TRUE),
      Cos = mean(Cos, na.rm = TRUE),
      Sil = mean(Sil, na.rm = TRUE)) |>
    as.data.frame()

  deep = deep[deep$Data == "Test", ]
  deep = deep[deep$Set == setting, ]
  deep = deep |>
    dplyr::group_by(Model, Set, Dim, NComp) |>
    dplyr::summarize(
      Time = mean(Time, na.rm = TRUE),
      Dev = mean(Dev, na.rm = TRUE),
      RSS = mean(RSS, na.rm = TRUE),
      Cos = mean(Cos, na.rm = TRUE),
      Sil = mean(Sil, na.rm = TRUE)) |>
    as.data.frame()

  flat = flat[, c("Model", "Dim", "Time", "Dev", "RSS", "Sil")]
  deep = deep[, c("Model", "NComp", "Time", "Dev", "RSS", "Sil")]

  colnames(flat) = c("Model", "Dim", "Time", "Deviance", "RMSE", "Silhouette")
  colnames(deep) = c("Model", "Dim", "Time", "Deviance", "RMSE", "Silhouette")

  # Change the column names
  df = rbind(flat, deep)

  df$Setting = factor(c(
    rep("$n \\uparrow, m \\uparrow, d = 5$", times = nrow(flat)),
    rep("$n = 5000, m = 500, d \\uparrow$", times = nrow(deep))),
    levels = c("$n \\uparrow, m \\uparrow, d = 5$", "$n = 5000, m = 500, d \\uparrow$"))

  # Reorganize the data-frame into a flat representation
  bigdf = data.frame(Model = c(), Setting = c(), Dim = c())
  bigdf = rbind(bigdf, as.data.frame(df[, c("Model", "Setting", "Dim")]))
  bigdf = rbind(bigdf, as.data.frame(df[, c("Model", "Setting", "Dim")]))
  bigdf = rbind(bigdf, as.data.frame(df[, c("Model", "Setting", "Dim")]))
  bigdf = rbind(bigdf, as.data.frame(df[, c("Model", "Setting", "Dim")]))

  bigdf$Values = c(log10(df$Time), 100*df$Deviance, 100*df$RMSE, df$Silhouette)
  bigdf$Variables = rep(c("Time", "Deviance", "Error", "Silhouette"), each = nrow(df))
  bigdf$Variables = factor(bigdf$Variables, levels = c("Time", "Deviance", "Error", "Silhouette"))

  parser = function(string) TeX(string)

  # Plot the grid of summaries
  plt = ggplot(data = bigdf, map = aes(x = Dim, y = Values, color = Model, fill = Model)) +
    # geom_bar(stat="identity", color = "white", position=position_dodge()) +
    geom_line() + geom_point(shape = 21, size = 4, colour = "white") +
    facet_grid(cols = vars(Setting), rows = vars(Variables), scales = "free",
               labeller = as_labeller(parser, default = label_parsed)) +
    theme_gray() + labs(x = "", y = "", fill = "Model", color = "Model") +
    theme(legend.position = "right",
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 13),
          strip.text = element_text(size = 13),
          axis.title = element_blank())
    # theme(legend.position = "right")

  if (SHOW) print(plt)

  if (SAVE) {
    filepath = join.path("img", "splatter")
    filename = join.string("summary_", setting,"_sim.pdf")
    zoom = 4
    nrows = 4
    ncols = 2
    width = 1 * ncols
    height = .5 * nrows
    ggplot2::ggsave(
      file = filename, path = filepath, plot = plt,
      width = zoom * width, height = zoom * height)
  }

  return (plt)
}

## SIMULATION RESULTS ----

flat.sim.list = load.flat.data()
deep.sim.list = load.deep.data()

flat.sim.df = merge.list(flat.sim.list)
deep.sim.df = merge.list(deep.sim.list)

plot.list(flat.sim.list)
plot.list(deep.sim.list)

plot.flat.data(flat.sim.df)
plot.deep.data(deep.sim.df)

plot.contrast(flat.sim.df, deep.sim.df, setting = "bubble")
plot.contrast(flat.sim.df, deep.sim.df, setting = "linear")
plot.contrast(flat.sim.df, deep.sim.df, setting = "branch")

## END OF FILE ----


