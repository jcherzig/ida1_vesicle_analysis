library(readxl)
library(writexl)
library(ggplot2)
library(pracma)
library(magrittr)
library(tidyverse)
library(gsignal)
library(patchwork)

# this script carries out the same analysis as 'Puncta_analysis_normalised_alldata_v2.R' but
# separates the neuron end and body

source("Puncta_analysis_normalised_neuron_end_functions.R")

n = 0.2 # set the fraction of total neuron length you want to analyse at the end

# set the fixed limits on the y axis
y_lim_1 <- -1000
y_lim_2 <- 25000

filelist_tail <- list.files(pattern = paste(".*_tail.xlsx",sep="|")) # we only wish to analyse the tail

lapply(filelist_tail,analyse_intensities)
