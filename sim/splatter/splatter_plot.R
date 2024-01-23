
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
})

## Some utility functions
join.path = function (...) paste(..., sep = "/")
join.string = function (...) paste(..., sep = "")

## SIMULATION RESULTS ----

fileid = c(
  "07-09-2023-11-42",
  "07-09-2023-12-59",
  "07-09-2023-13-26",
  "07-09-2023-13-56",
  "07-09-2023-14-26",
  "07-09-2023-14-56",
  "07-09-2023-15-28",
  "07-09-2023-16-50")

fileid = c("22-09-2023-10-21")

error = NULL
for (h in 1:length(fileid)) {
  filepath = join.path("tests", "splatter")
  filename = join.string("splatter_error_summary_", fileid[h], ".csv")
  error = rbind(error, read.table(file = join.path(filepath, filename), header = TRUE, sep = ";"))

}

models = c(
  "glmPCA-B",
  "glmPCA-M",
  "glmPCA-MB",
  "GMF-AIRWLS",
  "GMF-Newton",
  "GMF-M-SGD",
  "GMF-C-SGD",
  "GMF-B-SGD")

error$model = factor(error$model, levels = models)
error$model = as.factor(error$model)
error$time = as.numeric(error$time)
error$RSS = as.numeric(error$RSS)
error$Cos = as.numeric(error$Cos)
error$Dev = as.numeric(error$Dev)
error$iter = as.numeric(error$iter)


limits = log10(range(error$time[-c(1238, 1246)]))
plt = ggplot(data = error[-c(1238, 1246),], mapping = aes(x = model, color = set, fill = set)) +
  # lims(y = c(0, 1)) +
  labs(fill = "", color = "") + theme_bw()
plt = ggpubr::ggarrange(
  plt + geom_boxplot(mapping = aes(y = log10(time)), alpha = 0.5) +
    labs(x = "Models", y = "Execution time") + lims(y = limits),
  plt + geom_boxplot(mapping = aes(y = Dev), alpha = 0.5) +
    labs(x = "Models", y = "Residual deviance"),
  plt + geom_boxplot(mapping = aes(y = RSS), alpha = 0.5) +
    labs(x = "Models", y = "Residual sum of squares"),
  # plt + geom_boxplot(mapping = aes(y = Cos), alpha = 0.5) +
  #   labs(x = "Models", y = "Cosine distance"),
  nrow = 3, ncol = 1, legend = "right", common.legend = TRUE)

print(plt)

