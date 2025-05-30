library(glmmTMB)
library(ggplot2)
library(dplyr)
library(foreach)
library(patchwork)
library(here)
##############################################################
path <- paste0("real_data_analysis/atlass/results/")
##############################################################
fig_path = "fig/"
source("func.R")
##############################################################
# Define filenames to load
mod_files <- c("mod_rr.rds","mod_us.rds")
############################################################
# Load models 
models <- load_models(path, mod_files)
names(models)  =   c("rr", "us")
#####################################################
cc        =   commandArgs(trailingOnly  = TRUE)
i         =   as.integer(cc[1])
i = 2
#####################################################
mod   =  models[[i]]
condlik3 <- sum(dpois(model.response(model.frame(mod)),
                      lambda = fitted(mod),
                      log = TRUE))

lev_path  <- paste0(path,"leverage/",names(models)[i], "/")
files <-   list.files(lev_path, full.names = TRUE)

res = foreach(i = files, .combine = "c") %do% {readRDS(i)}
trace_hat <-   sum(res[res<1])
#####################################################
caic =  c(clik = condlik3, cdf = trace_hat, caic = 2*(-condlik3 + trace_hat))
saveRDS(caic, file = paste0(path,"caic_",names(models)[i],".rds"))
#####################################################
sub_levs  = res[res < 1]
sd_sub <- sd(res[res <1 ])
tt  = 2
dd  = model.frame(mod)

(nrow(dd)^2*(sd_sub^2))/tt^2
m_sub <- mean(sub_levs)
m_sub*nrow(dd)

