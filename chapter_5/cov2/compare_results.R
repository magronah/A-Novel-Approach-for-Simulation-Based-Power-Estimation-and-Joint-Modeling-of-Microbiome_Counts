setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
source("func2.R")
fig_path=   paste0("fig/")
############################################################
##power analysis
pvalues  =  pvalue_fun(dd)

######Fisher's Method#
pvals <- c(0.01, 0.02, 0.03) # your p-values from simulations
dd =  data.frame(matrix(abs(rnorm(10)),2,5))
apply(dd, 1, combined_z)
######Stouffer's Z (Z-transform method):