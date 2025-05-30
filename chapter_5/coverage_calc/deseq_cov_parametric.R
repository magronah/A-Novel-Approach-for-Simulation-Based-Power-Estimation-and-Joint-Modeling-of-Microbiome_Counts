setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")

library(foreach)

directory_path <-   "~/scratch/dataset/new_sim/coverage/deseq_res/" 
files          <-   list.files(directory_path, full.names = TRUE)

res = foreach(i = files, .combine = "cbind") %do% {
     res        =    readRDS(i)
     res$log2FoldChange
}

res            =    as.data.frame(res)
colnames(res)  =    paste0("nsim", 1:ncol(res))
rownames(res)  =    paste0("taxon",1:nrow(res))

saveRDS(res, file = paste0(getwd(),"/reproducible/new_sim/coverage_calc/data/deseq_cov_parametric.rds"))

