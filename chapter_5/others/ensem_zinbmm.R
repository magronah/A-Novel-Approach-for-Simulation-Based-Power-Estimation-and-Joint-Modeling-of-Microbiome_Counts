setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
source("reproducible/new_sim/func.R")
source("reproducible/new_sim/param.R")
source("reproducible/new_sim/Load_Packages.R")
source("reproducible/new_sim/initial_param0.R")
path = paste0(getwd(),"/reproducible/new_sim/data/")
####################################################################
data      =   readRDS(paste0(path,"otu_meta_list_withzi.rds"))
################################################################
cc      =   commandArgs(trailingOnly  = TRUE)
i       =   as.integer(cc[1]) 
dd      =   data[[i]]
################################################################
countdata  =   dd$countdata
met_dd     =   dd$met_data
met_dd$dummy = factor(1)

mod        =   mms(y = countdata, fixed = ~group,
                   random = ~ 1|dummy, zi_fixed = ~1,
                      data = met_dd, method = "zinb")

saveRDS(mod, file=paste0("~/scratch/dataset/new_sim/zinbmm_mod/mod",i, ".rds")) 
