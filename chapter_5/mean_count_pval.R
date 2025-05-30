setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
library(foreach)
source("reproducible/new_sim/func.R")
path    =    paste0(getwd(),"/reproducible/new_sim/data/")
##########################################################################  
if(FALSE){
  
directory_path <-   "~/scratch/dataset/new_sim/deseq_res2/" 
files          <-   list.files(directory_path, full.names = TRUE)

mean_count = foreach(i = files, .combine = "cbind") %do% {
  res        =    readRDS(i)
  res$baseMean
}

mean_count            =    as.data.frame(mean_count)
colnames(mean_count)  =    paste0("nsim", 1:ncol(mean_count))
rownames(mean_count)  =    paste0("taxon",1:nrow(mean_count))

saveRDS(mean_count, file = paste0(path,"mean_count_dd.rds"))

#######################################################
### deseq p value
pval_deseq = foreach(i = files, .combine = "cbind") %do% {
  res        =    readRDS(i)
  res$padj
}

pval_deseq            =    as.data.frame(pval_deseq)
colnames(pval_deseq)  =    paste0("nsim", 1:ncol(pval_deseq))
rownames(pval_deseq)  =    paste0("taxon",1:nrow(pval_deseq))

saveRDS(pval_deseq, file = paste0(path,"pval_deseq2.rds"))
#######################################################
### rr p value
rr_est      =  readRDS(paste0(path, "rr_parametric_new.rds"))
pval_rr     =  pval_cal_fun(rr_est)

saveRDS(pval_rr, file = paste0(path,"pval_rr.rds"))
###########################################################
### us p value
us_est      =  readRDS(paste0(path, "us_only_parametric.rds"))
pval_us     =  pval_cal_fun(us_est)

saveRDS(pval_us, file = paste0(path,"pval_us.rds"))
}
###########################################################
path1  <-   "~/scratch/dataset/new_sim/nbmm_pval/"
files1  <-   list.files(path1, full.names = TRUE)

res1 = foreach(i = files1, .combine ="cbind") %do% {
  readRDS(i)$padj
}

dd1            =  as.data.frame(res1)
rownames(dd1)  =  paste0("taxon",1:nrow(dd1))
colnames(dd1)  =  paste0("sim",1:ncol(dd1))
saveRDS(dd1, file = paste0(getwd(), "/reproducible/new_sim/data/pval_nbmm.rds"))
###########################################################
path2  <-   "~/scratch/dataset/new_sim/zinbmm_pval/"
files2 <-   list.files(path2, full.names = TRUE)

res2 = foreach(i = files2, .combine ="cbind") %do% {
  readRDS(i)$padj
}

dd2            =  as.data.frame(res2)
rownames(dd2)  =  paste0("taxon",1:nrow(dd2))
colnames(dd2)  =  paste0("sim",1:ncol(dd2))

saveRDS(dd2, file = paste0(getwd(), "/reproducible/new_sim/data/pval_zinbmm.rds"))
###########################################################
