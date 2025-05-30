source("reproducible/power/Load_Packages.R")
source("reproducible/power/simulate_fun.R")
source("reproducible/power/utils.R")
source("reproducible/power/fitting_fun.R")
source("reproducible/power/read_dataset.R")
path = "reproducible/power/datasets2/sample_size_cal/"

#countdata_list_obs  =  readRDS(file = paste0(path,"data.rds"))
#metadata_list_obs   =  readRDS(file = paste0(path,"metadata.rds"))
##############################################################
nsim = 100; notu=1000; disp_scale = 0.3
otu_data =  list()
nsamp_vec = seq(30,110,20)

for(i in 1:7){
  nam       =   names(dispersion_param_list)[i] 
  dispersion_param      =   dispersion_param_list[[i]]
  logmean_param         =   logmean_param_list[[i]]
  logfoldchange_param   =   logfoldchange_param_list[[i]]
  
  dd = foreach(k = 1:length(nsamp_vec)) %do% {
    source("reproducible/power/simulate_fun.R")
    nsamp     =   nsamp_vec[[k]]
    countdata_sim_fun(logmean_param, logfoldchange_param, 
                      dispersion_param,
                      nsamp_per_group=nsamp,
                      disp_scale = disp_scale,
                      ncont = NULL,ntreat = NULL,
                      notu, nsim, 
                      maxlfc_iter = 10000,
                      seed = 100)
    
  }
  names(dd) = paste0("sample", nsamp_vec)  
  otu_data[[i]] = dd
}
names(otu_data) = names(dispersion_param_list)  

#saveRDS(otu_data, file = paste0(path,"sim_data_list.rds"))

path1       =   "reproducible/power/datasets2/"
path_deseq  =    paste0(path1,"deseq_results/")

deseq_sample_size  = foreach(i = 1:7) %:%
  foreach(j = 1:length(nsamp_vec),.combine = "rbind") %do% {
  nam        =   names(dispersion_param_list)[i]
  nsam       =   nsamp_vec[j]
  ds         =   readRDS(paste0(path_deseq,nsam,"_samples/",nam,".rds"))
  deseq_list =   do.call(rbind, lapply(ds, function(x) x$deseq_estimate))
  deseq_list$sample_size = rep(nsam, nrow(deseq_list))
  deseq_list 
}

path_comb =    paste0(path1,"contour_data/")

true_data  = foreach(i = 1:7) %:%
  foreach(j = 1:length(nsamp_vec),.combine = "rbind") %do% {
    nam        =   names(dispersion_param_list)[i]
    nsam       =   nsamp_vec[j]
    ds         =   readRDS(paste0(path_comb,nsam,"_samples/combined_data_",nam,"_.rds"))
    ds$sample_size = rep(nsam, nrow(ds))
    ds 
  }

#############################################################
for(j in 1:7){
  nam       =   names(dispersion_param_list)[j] 
  #dd_otu    =   otu_data[[j]]
  deseq_est_list  =   deseq_sample_size[[j]]
  
  #sim_logmean_list      =   read_data(dd, "logmean_list") 
  #sim_logfoldchange_list       =   read_data(dd, "logfoldchange_list") 
  true_logmean = (true_data[[j]])$lmean_abund
  true_logfoldchange = true_data[[j]]$abs_lfc
  #length(true_logfoldchange)
  ################################################################
  power_est       =    power_fun_ss(deseq_est_list, true_logfoldchange,
                                    true_logmean, nsamp_vec, alpha=0.1)
  gam_fit_ss      =    power_est$gam_mod 
  combined_data   =    power_est$combined_data
  ################################################################
  path_gam <- paste0(path,"gam_fit/")
  if (!file.exists(path_gam)) {
    dir.create(path_gam, recursive = TRUE)
    print(paste("Folder", path_gam, "created successfully."))
  } else {
    print(paste("Folder", path_gam, "already exists."))
  }
  
  saveRDS(gam_fit_ss, file = paste0(path_gam,nam,".rds"))
  
  path_comb_data <- paste0(path,"combined_data/")
  if (!file.exists(path_comb_data)) {
    dir.create(path_comb_data, recursive = TRUE)
    print(paste("Folder", path_comb_data, "created successfully."))
  } else {
    print(paste("Folder", path_comb_data, "already exists."))
  }
  
  saveRDS(combined_data, file = paste0(path_comb_data,nam,".rds"))
}

#############################################################
## load the fitted models
path_gam      =   paste0(path,"gam_fit/")
path_simdata  =   paste0(path,"combined_data/")
gam_mod_list      =   comb_list  =   list()

for(i in 1:7){
  nam         =   names(filtered_otu_list)[i]
  comb_list[[i]]    =   readRDS(paste0(path_simdata,nam,".rds"))
  gam_mod_list[[i]] =   readRDS(paste0(path_gam,nam,".rds"))}

names(gam_mod_list)   =   names(filtered_otu_list)
names(comb_list)  =   names(filtered_otu_list)


## Simulate data
# nsam = 100
# path1  =  "reproducible/power/datasets/contour_data/"
# nam       =   names(dispersion_param_list)[1]
# combined_data = readRDS(paste0(path1,nsam,"_samples/","combined_data_",nam,"_.rds"))
# power_estimate = readRDS(paste0(path1,nsam,"_samples/","power_estimate_",nam,"_.rds"))
# 
# combined_data$pvalue_reject <- factor(combined_data$pval_reject)
# 
# combined_data$pvalue_reject <- factor(ifelse(combined_data$pvalue_reject == 0, 
#                                              "No", "Yes"))
# 
# cont_breaks =  seq(0,1,0.1)
# gg_2dimc <- (ggplot(combined_data)
#              + aes(lmean_abund, abs_lfc)
#              + rasterise(geom_point(aes(color = pvalue_reject), alpha = 0.5))
#              + xlab(TeX("$\\log_2$(mean counts)")) 
#              + ylab(TeX("|$\\log_2$(fold change)|")) 
#              + scale_colour_manual(values = c("black", "red")) 
#              + labs(color="reject null hypothesis") 
#              + geom_contour(data = power_estimate,
#                             aes(z=power),lwd=1,
#                             breaks = cont_breaks)
#              + geom_label_contour(data = power_estimate,
#                                   aes(z= power,label = sprintf("%.3f", after_stat(level))),
#                                   breaks = cont_breaks
#              )
#              + custom_theme(16)
# )
# gg_2dimc
# 

# combined_data_list  =list()
# for(j in 1:7){
#          gam_mod        =   gam_mod_list[[j]]
#          combined_data  =   comb_list[[j]]
#   combined_data$power   =  predict(gam_mod,
#                             type = "response",
#                             newdata = data.frame(sample_size = combined_data$sample_size,
#                                        lmean_abund = combined_data$lmean_abund,
#                                        abs_lfc     = combined_data$abs_lfc))
# combined_data_list[[j]]  =  combined_data
# }
#names(combined_data_list) =   names(filtered_otu_list)
#saveRDS(combined_data_list, file = paste0(path,"data_with_power.rds"))

#############################################################
#log fold change, up to 2 is realistic, 
# 5 log mean abundance is ok
#############################################################
#names(dd_list) = names(filtered_otu_list)

ddd = do.call(rbind, dd_list)

ggplot(dd,aes(x= sample_size, y= power, color = as.factor(lfc))) +
  geom_point(size = 2) +
  geom_line(lwd=1)  +
  labs(color = "logfoldchange")+ 
  #facet_wrap(~ dataset, scales = "free_x") +
  custom_theme(16)
#################################################################
dd_list2 = list() #list to store the results
samp  = 50; lfc = 5
for(k in 1:7){
  nam      =   names(filtered_otu_list)[k]
  gam_fit_ss  =   gam_mod_list[[k]]
  
  dd = foreach(i = 1:length(lfc) ,.combine = "rbind") %do%{
    
   abs_lfc = rep(lfc,len)
    pp  = uniroot_lmb(target_power=targ,sample_size=samp,abs_lfc=abs_lfc,
                   model=gam_fit_ss,xmin=-5,xmax=15)
    
    d    =    data.frame(sample_size = pp, lmb=lmb, 
                         lfc = abs_lfc, power=target)
    
    d$dataset = rep(nam, nrow(d))
    d
  }
  dd_list2[[k]] = dd
}

for(r in 1:length(target)){
  targ  =   target[r]
  
  pp=uniroot_lmb(target_power=targ,sample_size=50,abs_lfc=5,
                 model=gam_fit_ss,xmin=-5,xmax=15)
  print(pp)
}


for(r in 1:length(target)){
  targ  =   target[r]
  pp=uniroot_abs_lfc(target_power=targ,sample_size=100,lmean_abund=6,
                 model=gam_fit_ss,xmin=0,xmax=15)
  print(pp)
}




plot(gam_fit_ss2)


# 
# 
# res =  foreach(i = 1:2)%do%{
#     
#     source("reproducible/power/Load_Packages.R")
#     
#     metadata_list   =  metadata_list_sampvec[[i]]
#     countdata_list  =  countdata_list_sampvec[[i]]
#     
#     dds1 =  deseq_fun_est(metadata_list[1:50],countdata_list[1:50],
#                           num_cores=10, ref_name= "control")
#     gc() 
#     
#     dds2 =  deseq_fun_est(metadata_list[51:100],countdata_list[51:100],
#                           num_cores=10, ref_name= "control")
#     gc() 
#     
#     dds  =  c(dds1,dds2)
#     dds
#   }
#   gc() 
#   
#   res2 =  foreach(i = 3:4)%do%{
#     
#     source("reproducible/power/Load_Packages.R")
#     
#     metadata_list   =  metadata_list_sampvec[[i]]
#     countdata_list  =  countdata_list_sampvec[[i]]
#     dds1 =  deseq_fun_est(metadata_list[1:50],countdata_list[1:50],
#                           num_cores=10, ref_name= "control")
#     gc()
#     dds2 =  deseq_fun_est(metadata_list[51:100],countdata_list[51:100],
#                           num_cores=10, ref_name= "control")
#     gc()
#     dds  =  c(dds1,dds2)
#     dds
#   }
#   res3 = c(res,res2)
#   
#   names(res3)  =  paste0("sample", nsamp_vec) 
#   nam         =  names(dispersion_param_list)[i]
#   
# 
# 
#   
# 
