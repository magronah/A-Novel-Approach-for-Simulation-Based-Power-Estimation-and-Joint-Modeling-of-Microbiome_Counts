setwd("/home/agronahm/projects/def-bolker/agronahm/longitudinal_RR/")
library(glmmTMB)
library(ggplot2)
library(dplyr)
library(patchwork)
library(nlme)
library(MASS)
library(Matrix)
library(NBZIMM)
#library(here)
##############################################################
fig_path = "fig/"
source("func.R")
#########################################################
path_vec       =   c("real_data_analysis/atlass/results/",
                     "real_data_analysis/pregnancy/results/")

mod_names      =   c("us", "rr")
data_names     =   c("atlass","pregnancy")

group_lab_vec  =  c("bmi_groupobese:time", "pregnant1:GA_Days")

confint_list   =   confint_List  =  list()
for(i in 1:length(path)){
  path   =   path_vec[i]
  for(j in 1:length(mod_names)){
    mod_name           =   mod_names[j]
    mod                =   readRDS(paste0(path,"mod_",mod_name,".rds"))
    mean_count         =   readRDS(paste0(path,"mean_count.rds"))
    ntaxa              =   readRDS(paste0(path,"ntaxa.rds"))  
    
         confint       =   wald_confint2(mod, ntaxa =  ntaxa,
                                         mean_count =  mean_count,
                                         mod_name   =  mod_name, 
                                         path       =  path)
         
    saveRDS(confint,file = paste0(path,"confint_",mod_name,".rds"))
  }
}
######################################################################
if(FALSE){
group_lab_vec  =  c("bmi_groupobese:time", "pregnant:GA_Days")
for(i in 1:length(path_vec)){
  path          =    path_vec[i]
  group_label   =    group_lab_vec[i]
  
  mod_nb   <-  readRDS(paste0(path,"mod_nb.rds"))
  mod_znb  <-  readRDS(paste0(path,"mod_znb.rds"))
  mean_count         =   readRDS(paste0(path,"mean_count.rds"))
  ntaxa              =   readRDS(paste0(path,"ntaxa.rds")) 
  
  confint_nb <- nb_znb_confint(mod_nb, mean_count,
                               group_label = group_label,
                               conf_level = .95)
  
  saveRDS(confint_nb,file = paste0(path,"confint_nb.rds"))
  
  confint_znb <- nb_znb_confint(mod_znb, mean_count,
                                group_label = group_label,
                                conf_level = .95)
  
  saveRDS(confint_znb,file = paste0(path,"confint_znb.rds"))
}





library(RhpcBLASctl)
library(glmmTMB)
library(ggplot2)
library(dplyr)
library(patchwork)
library(AICcmodavg)
library(Matrix)
library(here)
##############################################################
source("func.R")
#########################################################
# Define filenames to load
filenames <- c(
  "mod_rr.rds",
  "mod_us.rds",
  "nbmm_aicc.rds",
  "zinbmm_aicc.rds"
)
######################################################################
# Load models for autism data
atlass_path <- paste0("atlass/results/")
pregnancy_path <- paste0("pregnancy/results/")

# Load models for dataset
atlass_models <- load_models(atlass_path, filenames)
pregnancy_models <- load_models(pregnancy_path, filenames)

# Assigning names 
mod_names <- c("RR","US","NB","ZNB")
names(atlass_models)  =  names(pregnancy_models)  =  mod_names
mod_list    =   lst(atlass_models,pregnancy_models)
###########################################################
## us and rr models only
us_rr_models <- lapply(mod_list, function(x){
  x[grep("^(RR|US)", names(x))]
})
###########################################################
atlass_us_rr      =   us_rr_models$atlass_models
pregnancy_us_rr   =   us_rr_models$pregnancy_models
###########################################################
atlass_us_rr      =   us_rr_models$atlass_models
pregnancy_us_rr   =   us_rr_models$pregnancy_models
############################################################
    ntaxa_file    <-  "ntaxa.rds"
atlass_ntaxa      <-  load_models(atlass_path, ntaxa_file)
pregnancy_ntaxa   <-  load_models(pregnancy_path,  ntaxa_file)
############################################################
datasets <- list(
  atlass  = list(model = atlass_us_rr, 
                 ntaxa = atlass_ntaxa, 
                 path = atlass_path),
  
  pregnancy = list(model = pregnancy_us_rr,  
                   ntaxa = pregnancy_ntaxa, 
                   path = pregnancy_path)
)
############################################################
#name = "atlass"
for (name in names(datasets)) {
  model <- datasets[[name]]$model
  ntax <- as.numeric(datasets[[name]]$ntaxa)
  path_dir <- datasets[[name]]$path
  
  mu_count  =   readRDS(paste0(path_dir,"mean_count.rds"))
  
  confint_us <- wald_confint2(mod  = model$US,
                             ntaxa = ntax,
                             mean_count = mu_count,
                             mod_name = "us",
                             path = path_dir)
  
  saveRDS(confint_us, file = paste0(path_dir, "CI_us.rds"))
  ###########################################################
  confint_rr <- wald_confint2(mod  = model$RR,
                             ntaxa = ntax,
                             mean_count = mu_count,
                             mod_name = "rr",
                             path = path_dir)
  
  saveRDS(confint_rr, file = paste0(path_dir, "CI_rr.rds"))
}
###########################################################
## nbmm and zinbmm models only
filenames2 <- c(
    "mod_nbmm.rds",
    "mod_zinbmm.rds"
  )
    
# Load models for dataset
atlass_models <- load_models(atlass_path, filenames2)
pregnancy_models <- load_models(pregnancy_path, filenames2)

# Assigning names 
mod_names <- c("NB","ZNB")
names(atlass_models)  =  names(pregnancy_models)  =  mod_names
mod_list    =   lst(atlass_models,pregnancy_models)

datasets_nbmm <- list(
  atlass  = list(model = atlass_models, path = atlass_path),
  pregnancy= list(model = pregnancy_models, path = pregnancy_path)
)
######################################################
# i = 1
for(i in 1:length(datasets_nbmm)){
  mod_nbmm   <-   datasets_nbmm[[i]]$model$NB
  mod_zinbmm <-   datasets_nbmm[[i]]$model$ZNB
  #################################################
  path       <-   datasets_nbmm[[i]]$path
  mu_count   <-   readRDS(paste0(path,"mean_count.rds"))
  #################################################
  name_vec   <-  mod_nbmm$variables$dist
  grp_label <-   grep("[:].*(time|days|day|week|weeks|times)", name_vec, 
                    ignore.case = TRUE, 
                    value = TRUE)
  #################################################
  CI_nbmm   <- zinbmm_confint(mod_nbmm, mean_count = mu_count,
                              group_label = grp_label,
                              conf_level = .95)

  CI_zinbmm  <- zinbmm_confint(mod_zinbmm, mean_count = mu_count,
                               group_label = grp_label,
                               conf_level = .95)
  #################################################  
  saveRDS(CI_nbmm, file   = paste0(path, "CI_nbmm.rds"))
  saveRDS(CI_zinbmm, file = paste0(path, "CI_zinbmm.rds"))
}
}
