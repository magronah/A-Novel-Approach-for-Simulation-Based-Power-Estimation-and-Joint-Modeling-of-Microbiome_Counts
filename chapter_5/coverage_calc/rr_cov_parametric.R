setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
library(foreach)
source("load_glmmTMB.R")


path  <-   "~/scratch/dataset/new_sim/coverage/rr_res/"
files <-   list.files(path, full.names = TRUE, pattern = "^param")

res = foreach(i = files) %do% {
     readRDS(i)
}

names(res)    =  paste0("parametric",1:length(res))
saveRDS(res, file = paste0(getwd(), "/reproducible/new_sim/coverage_calc/data/rr_cov_parametric.rds"))

