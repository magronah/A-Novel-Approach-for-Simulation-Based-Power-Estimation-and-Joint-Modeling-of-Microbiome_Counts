#setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
library(tibble)
library(DESeq2)

path        =   paste0(getwd(),"/reproducible/new_sim/real_data/")
source(paste0(path,"fun.R"))
source(paste0(path,"autism/load_data.R"))

res        =   deseqfun(dd_filt,met_dd,ref_name="ASD")
saveRDS(res, file=paste0(path,"autism/data/deseq_mod.rds"))


dim(dd_filt)
