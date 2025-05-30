setwd("/home/agronahm/projects/def-bolker/agronahm/longitudinal_RR/")
library(RhpcBLASctl)
library(Matrix)
library(huge)
library(tidyverse)
library(dplyr)
library(NBZIMM)
library(glmmTMB)
library(foreach)
library(doParallel)
###########################################################
source("func.R")
###########################################################
param_index   =   1
source("simulation_analysis/initial_param0.R")
path = paste0("simulation_analysis/sim_data/",nsubj,"_",ntaxa,"/")
path
###########################################################
nb             =    readRDS(paste0(path,"nb.rds"))
zinb           =    readRDS(paste0(path,"zinb.rds"))
us             =    readRDS(paste0(path,"us.rds"))
rr             =    readRDS(paste0(path,"rr.rds"))
true_param     =    readRDS(paste0(path,"true_param.rds"))
###########################################################
compute_ci <- function(est_col, sd_err_col, true_parm, 
                       conf_level,
                       taxa_names) {
  z_score <- qnorm(conf_level + (1 - conf_level) / 2)
  
  dd <- data.frame(
    est_param = est_col,
    lwr = est_col - z_score * sd_err_col,
    upr = est_col + z_score * sd_err_col,
    width = 2 * z_score * sd_err_col
  )

  dd$taxon  <-   taxa_names
  df       <-  inner_join(dd, true_param, by  =  "taxon")
  df$cov   <-  if_else(df$true_param >= df$lwr & df$true_param <= df$upr, 
                    1,0)
  df
}

sd_err  =   nb$sdr
est     =   nb$est
############################################################
conf_level <- 0.95
result_list <- lapply(seq_len(ncol(est)), function(i) {
  compute_ci(est[, i], sd_err[, i],
             true_param, 
             conf_level,
             taxa_names = rownames(est))
  
})
############################################################
ll     <-  length(result_list)
cov_dd <- (result_list
           %>% setNames(paste0("cov", 1:ll))
           %>% purrr::map_dfr(pull, cov) 
           %>% data.frame()
)
coverage   <-  mean(rowMeans(cov_dd))

############################################################
mod_res  =   lst(nb, zinb, rr, us)
power    =   as.data.frame(lapply(mod_res, function(x){
                                  mean(colMeans(x$padj))}))
############################################################
nn     =   10
ggplot(power, aes(model, power, group =1)) +
  geom_point()  +
   custom_theme(nn)
############################################################
View(result_list[[1]])
# Combine into a single data frame with column index identifier
result_df <- do.call(rbind, result_list)
result_df$column_id <- rep(seq_len(ncol(est)), each = nrow(est))
######################################################################

confint_zinb   =   function(mod, mean_count,
                            group_label = "grouptreat",
                            conf_level = .95){
  for(i in nrow(sd_err)){
    
  }
  dd         =   data.frame(est_param = est)
  z_score    =   qnorm(conf_level + (1 - conf_level)/2)
  
  dd$lwr     =   est  -  z_score*sd_err
  dd$upr     =   est  +  z_score*sd_err
  dd$width   =   dd$upr  - dd$lwr
  
  dd$pvalue      =   pvalue  
  dd$param_name  =   names(mod$fit)
  dd
}
