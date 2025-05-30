setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
source("reproducible/simp_single/func.R")
source("reproducible/simp_single/param.R")
source("reproducible/simp_single/Load_Packages.R")
source("reproducible/simp_single/initial_param0.R")
source("load_glmmTMB.R")
path = paste0(getwd(),"/reproducible/simp_single/data/")
################################################################
data =   readRDS(paste0(getwd(),"/reproducible/simp_single/data/otu_meta_dd_pow.rds"))
################################################################
cc      =   commandArgs(trailingOnly  = TRUE)
i       =   as.integer(cc[1]) 

countdata   =   data[[i]]$otu_table
met_dd      =   data[[i]]$meta_dd
res         =   deseqfun(countdata,met_dd,alpha_level=0.1,ref_name="control")

saveRDS(res, file=paste0("~/scratch/dataset/power/deseq_res/mod", i, ".rds"))  


