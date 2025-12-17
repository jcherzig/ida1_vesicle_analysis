library(readxl)
library(writexl)
library(ggplot2)
library(pracma)
library(magrittr)
library(tidyverse)
library(gsignal)
library(patchwork)

source("Puncta_analysis_normalised_alldata_functions.R")

# set the fixed limits on the y axis
y_lim_1 <- -1000
y_lim_2 <- 20000

# neuron tail and middle will likely require different axis limits so run separately
filelist_tail <- list.files(pattern = paste(".*_tail.xlsx",sep="|"))
filelist_middle <- list.files(pattern = paste(".*_middle.xlsx",sep="|"))

# filelist_all <- list.files(pattern = paste(".*_tail.xlsx",".*_middle.xlsx",sep="|"))


lapply(filelist_middle,analyse_intensities)
