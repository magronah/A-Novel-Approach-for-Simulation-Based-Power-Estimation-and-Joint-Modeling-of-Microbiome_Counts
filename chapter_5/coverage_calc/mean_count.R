setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
path = paste0(getwd(),"/reproducible/new_sim/coverage_calc/data/")
###################################################################
dd        =    readRDS(paste0(path,"otu_meta_list_nozi.rds"))
count_dd  =    lapply(dd, function(x){as.data.frame(x["countdata"])})
colmean   =    as.data.frame(lapply(count_dd, function(x){colMeans(x)}))
saveRDS(colmean, file = paste0(path,"mean_counts.rds"))
#####################################################################
