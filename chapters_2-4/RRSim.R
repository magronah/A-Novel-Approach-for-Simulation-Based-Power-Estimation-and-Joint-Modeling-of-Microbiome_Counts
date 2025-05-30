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
##########################################################
path  =  "reproducible/power/datasets2/"
countdata_list_obs  =  readRDS(file = paste0(path,"data.rds"))
metadata_list_obs   =  readRDS(file = paste0(path,"metadata.rds"))
###################################################
cc =   commandArgs(trailingOnly  = TRUE)
i         =   as.integer(cc[1])

countdata  =   countdata_list_obs[[i]]
metadata    =   metadata_list_obs[[i]]
names(metadata) <- c("subject", "group")
ntaxa =  nrow(countdata)
##########################################################
#filter low abundance taxa
filter_low_otu    = filter_fun(countdata,metadata,abund_thresh=5,
                                 sample_thresh=3,
                                 subject_label = "subject",
                                 group_label = "group",
                                 ntaxa)
  
countdata_filt    =  filter_low_otu$countdata_filt
##########################################################
dlong    =    dd_long_fun(countdata_filt, metadata, subject_label = "subject")
ddlong   =    left_join(dlong, metadata, by ="subject")
#############################################################  
mm     =   dd_wide_fun(ddlong)
mm_dd  =  (mm$countdata)
#############################################################    
for(i in 1:nrow(countdata_filt)){
    test=ddlong %>% 
      filter(taxon == (ddlong$taxon)[i])%>% 
      data.frame() %>% 
      select(count)
     stopifnot(unlist(countdata_filt[i,]) == mm_dd[i,])
}

stopifnot(as.numeric(unlist(test)) == as.numeric(unlist(countdata_filt[i,])))
##########################################################
dflong  =   otu_meta_fun(ddlong)
mod     <-  RR_fit(dflong)
sim_rr  =  RRSim(mod, dflong,nsim = 1, seed =100)
rr_dd = list(mod = mod, sim = sim_rr)
  
nam  =   names(countdata_list_obs)[i]
saveRDS(rr_dd, file = paste0(path,nam,"sim_mod_rrzi.rds"))
####################################################################
sub = c(1,2,6,7)
if(FALSE){
  path1  =  "reproducible/power/datasets/"
countdata_sim_list =   readRDS(paste0(path1,"countdata_sim_compare.rds"))
rr_sim_lst        =   readRDS(paste0(path1,"sim_mod_rr.rds"))
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
plt =  foreach(i = 1:length(HMP_list), .errorhandling = "pass") %do% {
  nam       =    names(filtered_otu_lst)[i]
  countdata_obs     =    filtered_otu_lst[[i]]
  
    HMP       =    HMP_list[[i]]
  metaSPARSim =    metaSPARSim_list[[i]]
  scaled      =    scaled_list[[i]] 
  #rrsim       =    rr_sim_list[[i]][[1]]

  countdata_sim       =    list(HMP = HMP, metaSPARSim = metaSPARSim,
                                MixGaussSim =  scaled#,
                                #RRSim   =  rrsim
                                )
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

mean_plt   = read_data(plt, "mean_plt")
var_plt    = read_data(plt, "var_plt")
zeros_plt  = read_data(plt, "compare_zeros")


(mean_plt[[1]]|var_plt[[1]]) +   plot_layout(guides = "collect")  
(mean_plt[[2]]|var_plt[[2]]) +   plot_layout(guides = "collect")  
(mean_plt[[3]]|var_plt[[3]]) +   plot_layout(guides = "collect")  &
  theme(legend.position='bottom')



best_mean = read_data(plt, "best_mean")
best_var  = read_data(plt, "best_var")

ks_mean = read_data(plt, "ks_mean_est")
ks_var  = read_data(plt, "ks_var_est")
####################################################################
(zeros_plt[[1]]|zeros_plt[[2]]|zeros_plt[[3]]) +
  plot_layout(guides = "collect")  &
  theme(legend.position='bottom')
#It

mu_plt  = (mean_plt[[1]]|mean_plt[[2]]|mean_plt[[3]]) + 
        plot_layout(guides = "collect")   &
  theme(legend.position='bottom')

var_plt  = (var_plt[[1]]|var_plt[[2]]|var_plt[[3]]) +  
  plot_layout(guides = "collect")  &
  theme(legend.position='bottom')

width =  15; height =  5; dpi = 300
ggsave("reproducible/power/figures/mean_dist.png", plot = mu_plt, 
       width = width, 
       height = height, 
       dpi = dpi)

ggsave("reproducible/power/figures/var_dist.png", plot = var_plt, 
       width = width, 
       height = height, 
       dpi = dpi)
####################################################################
nt =  lapply(countdata_sim, function(x){
  nt_subjects <- metadata_s %>% 
    dplyr::filter(group == "ASD") %>% 
    dplyr::pull(subject)
  x[, colnames(x) %in% nt_subjects]
  
}
  )

# Extract subjects labeled as "NT" from metadata_s
nt_subjects <- metadata_s %>% 
  dplyr::filter(group == "ASD") %>% 
  dplyr::pull(subject)

# Subset countdata_sim$sim1 to include only the columns corresponding to NT subjects
countdata_nt <- countdata_filt[, colnames(countdata_filt) %in% nt_subjects]
p1   =   compare_dataset(nt,countdata_nt,method = "var") 
p2   =   compare_dataset(nt,countdata_nt,method = "mean") 

p1|p2


#####################################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
source("reproducible/power/utils.R")

ppr = plt[[1]]$ks_mean_est
names(plt$pp)

nn = 20; wd =  1.5
ks_data <- data.frame(
  Dataset = rep(c("PRJNA168470", "PRJNA589343", "PRJNA687773"), each = 4),
  Method = rep(c("HMP", "metaSPARSim", "MixGaussSim", "RRSim"), times = 3),
  KS_Statistic = c(0.106, 0.060, 0.069, 0.065,
                   0.189, 0.042, 0.058, 0.210,
                   0.242, 0.027, 0.048, 0.222),
  p_value = c(0.00098, 0.1687, 0.0777, 0.1160,
              1.10e-09, 0.6739, 0.2582, 6.17e-12,
              3.18e-23, 0.9052, 0.2540, 1.34e-19)
)

# Ensure Dataset and Method are factors
ks_data$Dataset <- factor(ks_data$Dataset,
                          levels = c("PRJNA168470", "PRJNA589343", 
                                     "PRJNA687773"))
ks_data$Method <- factor(ks_data$Method, 
                         levels = c("HMP", "metaSPARSim", 
                                    "MixGaussSim", "RRSim"))

# Plot KS Statistic and P-value using lines
p1 = ggplot(ks_data, aes(x = Dataset, y = p_value, color = Method )) +
  geom_line(aes(x = Dataset, y = p_value, group = Method ), linewidth = wd) +
  geom_point(size = 3) +
  labs(x = " ") +
  custom_theme(nn)

p11 = ggplot(ks_data, aes(x = dataset, y = stats, color = Method )) +
  geom_line(aes(x = Dataset, y = KS_Statistic, group = Method ), linewidth = wd) +
  geom_point(size = 3) +
  labs(x = " ") +
  custom_theme(nn) 


p111 = (p1|p11) +  plot_layout(guides = "collect") #& theme(legend.position='bottom')
width =  18; height =  7; dpi = 300

ggsave("reproducible/power/figures/ks_logmean.png", plot = p111, 
       width = width, 
       height = height, 
       dpi = dpi)


ks_data2 <- data.frame(
  Dataset = rep(c("PRJNA168470", "PRJNA589343", "PRJNA687773"), each = 4),
  Method = rep(c("HMP", "metaSPARSim", "MixGaussSim", "RRSim"), times = 3),
  KS_Statistic = c(0.207, 0.451, 0.129, 0.138,
                   0.297, 0.290, 0.097, 0.250,
                   0.308, 0.284, 0.097, 0.236),
  P_Value = c(4.01e-13, 1.28e-60, 2.27e-5, 4.55e-6,
              2.13e-23, 2.24e-22, 0.0073, 9.72e-17,
              2.63e-37, 6.58e-32, 0.00043, 3.47e-22)
)

# Ensure Dataset and Method are factors
ks_data2$Dataset <- factor(ks_data2$Dataset,
                          levels = c("PRJNA168470", "PRJNA589343", 
                                     "PRJNA687773"))
ks_data2$Method <- factor(ks_data2$Method, 
                         levels = c("HMP", "metaSPARSim", 
                                    "MixGaussSim", "RRSim"))

# Plot KS Statistic and P-value using lines
p2 = ggplot(ks_data2, aes(x = Dataset, y = P_Value, color = Method )) +
  geom_line(aes(x = Dataset, y = P_Value, group = Method ), linewidth = wd) +
  geom_point(size = 3) +
  #scale_y_log10() +  # Use a log10 scale for better visualization
  labs(x = " ",y = "P-Value (log scale)") +
  custom_theme(nn)

p22 = ggplot(ks_data2, aes(x = Dataset, y = KS_Statistic, color = Method )) +
  geom_line(aes(x = Dataset, y = KS_Statistic, group = Method ), linewidth = wd) +
  geom_point( size = 3) +
  labs(x = " ") +
  custom_theme(nn)

p222 = (p2|p22) +  plot_layout(guides = "collect") 
#width =  15; height =  5; dpi = 300
ggsave("reproducible/power/figures/ks_logvar.png", plot = p222, 
       width = width, 
       height = height, 
       dpi = dpi)

}
