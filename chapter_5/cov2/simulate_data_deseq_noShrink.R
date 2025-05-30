setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")
library(reformulas)
library(Matrix)
library(huge)
library(tidyverse)
library(dplyr)
library(glmmTMB)
library(foreach)
library(tibble)
library(DESeq2)
source("func2.R")
source("initial_param0.R")
####################################################################
path1 = paste0(nsubj,"_",ntaxa,"/")
path1
####################################################################
data	  =   readRDS(paste0(path1,"otu_meta_list_withzi_taxa.rds"))
################################################################
p  =   100 
## Note we are using the 100th simulation as the refernce model. 
dd	=   data[[p]]
################################################################
countdata  =   dd$countdata
met_dd     =   dd$met_data
res        =    deseqfun(countdata,met_dd,ref_name="control",do_shrinkage = "no")
################################################################
path  =  paste0("cov2/",nsubj,"_",ntaxa,"/")

if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE) 
  cat("Folder created at:", path, "\n")
} else {
  cat("Folder already exists at:", path, "\n")
}

####################################################################
dd     =   meta_data(ntaxa, nsubj)
sim_list = list()
for(i in 1:ntaxa) {
  
  beta      =  log(2)*as.numeric(res[i,][c("intercept","log2FoldChange")])
  betadisp  =   1/as.numeric(res[i,]["dispersion"])
  sim_pars  =   lst(beta, betadisp)
  
  ddd       =   dd %>% filter(taxon == i)
  
  ####################################################################
  ## check if rr is being mapped correctly
  seed      =   i
  form      =   formula(paste0("taxon", i, "~ group"))
  pars      =   simulate_new(RHSForm(form, as.form = TRUE), nsim = 1, seed = seed,
                             newdata    =  ddd,
                             newparams  =  sim_pars,
                             family     =  nbinom2,
                             return_val = "pars")
  
  #####################################################################
  #now simulate data using b_est
  sim_list[[i]]  =  simulate_new(RHSForm(form, as.form = TRUE), nsim = nsim, seed = seed,
                                 newdata   =  ddd,
                                 newparams =  sim_pars,
                                 family    =  nbinom2,
                                 return_val = "sim")
  
}
names(sim_list)  =  paste0("taxon",1:ntaxa)

result <- lapply(1:nsim, function(j) {
  data.frame(do.call(cbind, lapply(sim_list, function(lst) lst[[j]])))
})

sim_dd  <- lapply(1:nsim, function(k) { 
  list(countdata =  result[[k]], 
       met_data   =  ddd)})

names(sim_dd)  =  paste0("sim", 1:nsim)
saveRDS(sim_dd, file=paste0(path,"sim_data/deseq_noShrink_otu_meta_list_withzi_taxa.rds"))
##############################################################################
# for(r in 1:ntaxa){
#   pp=sim_list[[r]][[1]]
#   stopifnot(sim_dd[[1]][[1]][,r] == pp)
# }


