library(readxl)
library(writexl)
library(ggplot2)
library(pracma)
library(magrittr)
library(tidyverse)
library(gsignal)
library(patchwork)

### NOTE:
# for this script, the sheet name to be analysed has to be manually set in the lapply call at the bottom
# typically this is 'axon' or 'dendrite'

source("Puncta_analysis_normalised_alldata_nlp21_functions.R")

# set the fixed limits on the y axis
y_lim_1 <- -1000
y_lim_2 <- 6000

filelist_all <- list.files(pattern = paste(".*_nlp21.xlsx",sep="|"))

lapply(filelist_all,analyse_intensities,sheet="axon")
