setwd("/home/agronahm/projects/def-bolker/agronahm/Michael-n-Ben-Repo/")
library(DESeq2)
library(tidyverse)
library(dplyr)
library(rlist)
library(RhpcBLASctl)
library(glmmTMB)
library(gtools)
library(foreach)
library(latex2exp)
library(reformulas)
###############################################################
source("reproducible/power/utils.R")
source("reproducible/power/fitting_fun.R")
path1  =  "reproducible/power/datasets/"
countdata_sim_list =   readRDS(paste0(path1,"countdata_sim_compare.rds"))

path = "reproducible/power/datasets/"
###########################Read otu dataset###########################
otu_dataset_list     =    readRDS(paste0(path,"otu_dataset_list.rds"))
filtered_otu_list    =    otu_dataset_list[["filtered_otu"]]
sub =  c("PRJNA168470", "PRJNA355023", "PRJNA589343","PRJNA687773")
rr_sim_lst        =   foreach(i = sub) %do%{
  readRDS(paste0(path1,i,"_rr.rds"))
}
####################################################################
HMP_lst          =    read_data(countdata_sim_list, "HMP")
metaSPARSim_lst  =    read_data(countdata_sim_list, "metaSPARSim")
scaled_lst       =    read_data(countdata_sim_list, "scaled")
rr_sim_list       =    lapply(rr_sim_lst, function(x){x$sim$countdata_list})
##########################################################
filtered_otu_lst   =   filtered_otu_list[sub]
HMP_list           =   HMP_lst[sub]
metaSPARSim_list   =   metaSPARSim_lst[sub]
scaled_list        =   scaled_lst[sub]

plt =  foreach(i = 1:length(HMP_list)) %do% {
  nam       =    names(filtered_otu_lst)[i]
  countdata_obs     =    filtered_otu_lst[[i]]
  
  HMP       =    HMP_list[[i]]
  metaSPARSim =    metaSPARSim_list[[i]]
  scaled      =    scaled_list[[i]] 
  rrsim       =    rr_sim_list[[i]][[1]]
  
  countdata_sim       =    list(HMP = HMP, metaSPARSim = metaSPARSim,
                                MixGaussSim =  scaled,
                                RRSim   =  rrsim)
  ###################################################################
  real_mean  =   rowMeans(countdata_obs)
  real_var   =   rowVars(as.matrix(countdata_obs))
  ###################################################################
  mm    =  lapply(countdata_sim, function(x){rowMeans(x)})
  varr  =  lapply(countdata_sim, function(x){rowVars(as.matrix(x))})
  ###################################################################
  ks_mean  <-  lapply(mm, function(x){ks.test(real_mean, x)})
  ks_var   <-  lapply(varr, function(x){ks.test(real_var, x)})
  ###################################################################
  ks_mean_est <-  lapply(ks_mean, function(x){c(stat = x$statistic, pval = x$p.value)})
  ks_var_est  <-  lapply(ks_var, function(x){c(stat = x$statistic, pval = x$p.value)})
  ###################################################################
  best_mean_pval <- which.max(unlist(lapply(ks_mean_est, function(x){x["pval"]})))
  best_mean_stat <- which.min(unlist(lapply(ks_mean_est, function(x){x["stat.D"]})))
  
  best_var_pval <- which.max(unlist(lapply(ks_var_est, function(x){x["pval"]})))
  best_var_stat <- which.min(unlist(lapply(ks_var_est, function(x){x["stat.D"]})))
  ###################################################################
  okabe_ito_colors = c("#556B2F", "#E23D28", "#0000FF","#E69F00","#000000")
  p1   =   compare_dataset(countdata_sim,countdata_obs,method = "var") 
  p1   =   p1  + ggtitle(nam) + scale_color_manual(values = okabe_ito_colors)   
  
  p2   =   compare_dataset(countdata_sim,countdata_obs,method = "mean")
  p2   =   p2  + ggtitle(nam) + scale_color_manual(values = okabe_ito_colors)   
  
  p3   =   compare_dataset(countdata_sim,countdata_obs,method = "compare_zeros")
  p3   =   p3 + ggtitle(nam) + scale_color_manual(values = okabe_ito_colors)  
  
  pp   =   list(var_plt  =  p1,  
                mean_plt = p2, 
                compare_zeros  = p3,
                ks_mean_est =  ks_mean_est,
                ks_var_est = ks_var_est,
                best_mean =  c(pval =  best_mean_pval,
                               stat =  best_mean_stat),
                best_var  = c(pval  =  best_var_pval,
                              stat  =  best_var_stat)
  )
  pp
}


ks_mean = read_data(plt, "ks_mean_est")
ks_var  = read_data(plt, "ks_var_est")

names(ks_mean) =  sub

px = lapply(ks_mean,  function(pr){
  t(data.frame(pr)) %>% 
    data.frame %>% 
    rownames_to_column("Method")}
  )

ppr = do.call(rbind, px)
ppr$dataset  =  rep(sub, each = 4)

p1 = ggplot(ppr, aes(x = dataset, y = (pval), color = Method )) +
  geom_line(aes(x = dataset, y = pval, group = Method ), linewidth = 1) +
  geom_point(size = 3) +
  labs(x = " ", y = "P_Value (log scale)") + 
  scale_y_log10() +
  custom_theme(10)

p11 = ggplot(ppr, aes(x = dataset, y = stat.D, color = Method )) +
  geom_line(aes(x = dataset, y = stat.D, group = Method ), linewidth = 1) +
  geom_point(size = 3) +
  labs(x = " ") +
   scale_y_log10() +
  custom_theme(10) 
######################################################################
# names(ks_mean) <- c("PRJNA168470", "PRJNA355023", "PRJNA589343", "PRJNA687773")

# Combine into one data.frame
df <- do.call(rbind, lapply(seq_along(ks_mean), function(i) {
  dataset_name <- names(ks_mean)[i]
  ks <- ks_mean[[i]]
  
  # For each method in the dataset
  do.call(rbind, lapply(names(ks), function(method) {
    data.frame(
      dataset = dataset_name,
      Method = method,
      KS_Statistic = ks[[method]]["stat.D"],
      p_value = ks[[method]]["pval"],
      row.names = NULL
    )
  }))
}))


p1 = ggplot(df, aes(x = dataset, y = p_value, color = Method )) +
  geom_line(aes(x = dataset, y = p_value, group = Method ), linewidth = 1) +
  geom_point(size = 3) +
  labs(x = " ") + 
  scale_y_log10() +
  custom_theme(10)

p11 = ggplot(df, aes(x = dataset, y = KS_Statistic, color = Method )) +
  geom_line(aes(x = dataset, y = KS_Statistic, group = Method ), linewidth = 1) +
  geom_point(size = 3) +
  labs(x = " ") +
  scale_y_log10() +
  custom_theme(10) 

library(patchwork)
p1|p11
