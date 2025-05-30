path        =   paste0(getwd(),"/reproducible/new_sim/real_data/")
library(DESeq2)
source(paste0(path,"fun.R"))
##############################################################
data        =   readRDS(paste0(path,"autism/data/data.rds"))
met_data    =   readRDS(paste0(path,"autism/data/metadata.rds"))
##############################################################
#use PRJNA168470, PRJNA589343 and PRJNA687773
names(data)
nam         =    "PRJNA687773"
count_dd    =   data[[nam]]
met_dd      =   met_data[[nam]] 
names(met_dd)  = c("subject", "group")
dd_filt      =   filter_fun(count_dd, met_dd,abund_thresh=10, sample_thresh=5)
##############################################################
# dd_filt_list = list()
# for(i in 1:7){
#  count_dd     =   data[[i]]
#  met_dd       =   met_data[[i]]
#  names(met_dd)  = c("subject", "group")
#  dd_filt_list[[i]]   =   filter_fun(count_dd, met_dd,abund_thresh=5, sample_thresh=3)
# }
# 
# lapply(dd_filt_list, dim)
##############################################################

