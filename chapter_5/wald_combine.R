setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
library(foreach)
source("load_glmmTMB.R")


directory_path <-   "~/scratch/dataset/glmmTMB_wald/"
files          <-   list.files(directory_path, full.names = TRUE)

res = foreach(i = files, .combine ="rbind") %do% {
   readRDS(i)
}

dd =  as.data.frame(res)

saveRDS(dd, file = paste0(getwd(), "/reproducible/simp_single/data/wald_conf.rds"))


###


