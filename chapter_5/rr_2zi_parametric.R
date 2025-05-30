setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
library(foreach)
source("load_glmmTMB.R")

directory_path <-   "~/scratch/dataset/new_sim/rr_mod2zi/"
files          <-   list.files(directory_path, full.names = TRUE)

res = foreach(i = files, .combine ="cbind") %do% {
  mod  =  readRDS(i)
  (ranef(mod)$cond$taxon$grouptreat)
}

dd =  as.data.frame(res)
rownames(dd)  =  paste0("taxon",1:nrow(dd))
colnames(dd)  =  paste0("sim",1:ncol(dd))


saveRDS(dd, file = paste0(getwd(), "/reproducible/new_sim/data/rr2zi_parametric.rds"))


