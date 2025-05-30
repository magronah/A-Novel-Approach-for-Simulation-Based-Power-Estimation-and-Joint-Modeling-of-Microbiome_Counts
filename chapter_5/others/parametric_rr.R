setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
library(foreach)
source("load_glmmTMB.R")

path  <-   "~/scratch/dataset/new_sim/rr_res/"
files <-   list.files(path, full.names = TRUE)

res = foreach(i = files, .combine ="cbind") %do% {
     readRDS(i)
}

dd            =  as.data.frame(res)
rownames(dd)  =  paste0("taxon",1:nrow(dd))
colnames(dd)  =  paste0("sim",1:ncol(dd))

saveRDS(dd, file = paste0(getwd(), "/reproducible/new_sim/rr_parametric.rds"))


