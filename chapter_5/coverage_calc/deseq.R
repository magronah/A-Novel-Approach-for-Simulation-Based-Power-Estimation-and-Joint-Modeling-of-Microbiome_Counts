setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
source("reproducible/new_sim/func.R")
source("reproducible/new_sim/param.R")
source("reproducible/new_sim/Load_Packages.R")
source("reproducible/new_sim/initial_param0.R")
path = paste0(getwd(),"/reproducible/new_sim/coverage_calc/data/")
###########################################################
data      =   readRDS(paste0(path,"otu_meta_list_nozi2.rds"))
################################################################
cc      =   commandArgs(trailingOnly  = TRUE)
i       =   as.integer(cc[1]) 
dd      =   data[[i]]
################################################################
countdata  =   dd$countdata
met_dd     =   dd$met_data

res        =   deseqfun(countdata,met_dd,alpha_level=0.1,ref_name="control")
saveRDS(res, file=paste0("~/scratch/dataset/new_sim/coverage/deseq_res/mod",i, ".rds")) 

