setwd("/home/agronahm/projects/def-bolker/agronahm/longitudinal_RR/")
library(RhpcBLASctl)
library(Matrix)
library(huge)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(glmmTMB)
library(foreach)
#library(here)
#########################################################
source("func.R")
#########################################################
#for(param_index  in  1:2){
param_index  =  2
source("simulation_analysis_new/initial_param0.R")
#########################################################
path  =  "simulation_analysis_new/sim_data/"
#########################################################
dd       =   metadata(ntaxa, nIndiv=nsubj, ntime) 
dd$time  =   dd$time/ntime
#########################################################
parms    =   lst(beta,betazi,theta,thetazi)
#########################################################
#simulate from the full model
form <- ~1 +  
   #ar1(factor(time) + 0| taxon:subject) +
   (group*time|taxon) +
   (1|nugget) +
   (1|subject:taxon) +
   (taxon + 0 | subject:time) 
# conclusion: remove the ar1 term and argue that there might 
# be confounding the trend term. 
# leave the rest the same
#########################################################
#### simulate fixed bs
pars0  <- simulate_new(form, nsim = 1, seed = seed,
                          newdata = dd,
                          ziformula  = ~1 + (1|taxon),
                          newparams = parms,
                          return_val = "pars",
                          family = poisson)

stopifnot(theta == pars0[names(pars0)=="theta"])
b_true_all      =   pars0[names(pars0)=="b"]
true_b_index    =   seq(n_gt,(n_gt*ntaxa),n_gt)
true_param      =   b_true_all[true_b_index]
true_param_dd   =   data.frame(param_name	  =  paste0("taxon",1:ntaxa),
                               true_param =  true_param)
#########################################################
#now simulate with bs fixed to bval
true_pars1  =  c(parms, lst(b = b_true_all))
pars1 <- simulate_new(form, nsim = 1, seed = NULL,
                      newdata    =  dd,
                      ziformula  = ~1 + (1|taxon),
                      newparams  =  true_pars1,
                      family     =  poisson, 
                      return_val = "pars"
)
stopifnot(b_true_all  == pars1[names(pars1) == "b"])
####################################################################  
sim_data <- simulate_new( form,
                           newdata   =  dd,
                           ziformula  = ~1 + (1|taxon),
                           newparams =  true_pars1,
                           return_val = "sim",
                           family    =  poisson,
                           nsim      =  nsim,
                           seed   =  seed)
###############################################################
res_dd = foreach(i = 1:length(sim_data))  %do% {
  dd$count = sim_data[[i]]
  dd
}
names(res_dd)  =   paste0("sim", 1:nsim)
###############################################################
split_dd = lapply(res_dd, function(lst){
  res = split(lst, lst$time)
  otu_meta_lst_fun(res,ntime)
  })
###############################################################
add_norm_list <- lapply(split_dd, function(sim) {
  df  = lapply(names(sim), function(time_name) {
    time <- sim[[time_name]]   
    otu_table  <- time$countdata
    otu_count  <- t(otu_table)
    meta_data   <- time$met_data
    
    dds <- DESeqDataSetFromMatrix(otu_count, meta_data, ~group)
    dds <- DESeq(dds, sfType = "poscounts", minReplicatesForReplace = Inf)
    normalizer <- sizeFactors(dds)
    
    time_num <- as.numeric(gsub("time", "", time_name))/ntime
    
    dd <- data.frame(
      subject    =  factor(names(normalizer), levels = 1:nsubj), 
      normalizer = normalizer,
      time       = time_num
    )
    metadata   =  left_join(meta_data,dd, by = "subject")  
    list(metadata   =  metadata,
        countdata  =  otu_table)
  })
  names(df) = names(sim)
  df
})
###########################################################
combine_metadata <- function(sim_data) {
  bind_rows(lapply(sim_data, function(time_point) time_point$metadata))
}
metadata_combined <- lapply(add_norm_list, combine_metadata)

combine_countdata <- function(sim_data) {
  bind_rows(lapply(sim_data, function(time_point) time_point$countdata))
}
countdata_combined <- lapply(add_norm_list, combine_countdata)
##################################################################
dd_list = list()
for(i in 1:length(metadata_combined)){
    dd1  =  metadata_combined[[i]]
    dd2  =  res_dd[[i]]
    dd_list[[i]]=left_join(dd1,dd2, by = c("subject","time","group"))
}
names(dd_list)  = paste0("sim",1:nsim)
##################################################################
file_path  =  paste0(path,nsubj,"_",ntaxa,"_",ntime,"/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}
##################################################################
saveRDS(dd_list, file=paste0(file_path,  "dd_long.rds"))
saveRDS(countdata_combined, file=paste0(file_path, "countdata.rds"))
saveRDS(metadata_combined, file=paste0(file_path, "metadata.rds"))
saveRDS(true_param_dd, file=paste0(file_path, "true_param.rds"))
#}
