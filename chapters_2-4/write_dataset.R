source("reproducible/power/Load_Packages.R")
source("reproducible/power/fitting_fun.R")
source("reproducible/power/simulate_fun.R")
source("reproducible/power/utils.R")
source("reproducible/power/read_dataset.R")
######################################################
path  =  "reproducible/power/datasets2/"

countdata_list_obs  =  readRDS(file = paste0(path,"data.rds"))
metadata_list_obs   =  readRDS(file = paste0(path,"metadata.rds"))
names(countdata_list_obs)
lapply(countdata_list_obs, dim)
######################################################
res = foreach(i = 1:7) %do% {
  source("reproducible/power/Load_Packages.R")
  source("reproducible/power/fitting_fun.R")
  source("reproducible/power/utils.R")
  
  nam  =   names(metadata_list_obs)[i]
  countdata  =   countdata_list_obs[[i]]
  metadata   =   metadata_list_obs[[i]]
  ######################################################
  #filter low abundance taxa
  filter_low_otu    = filter_fun(countdata,metadata,abund_thresh=5,
                                 sample_thresh=3)
  
  countdata_filt    =  filter_low_otu$countdata_filt
  control_count     =  filter_low_otu$control_count
  treat_count       =  filter_low_otu$treat_count

  #run deseq for log fold change estimates 
  deseq          =   deseqfun(countdata_filt,metadata,ref_name="NT")
  countdata_norm =   deseq$normalised_count  
  dispersion     =   deseq$dispersion
  logfoldchange  =   deseq$logfoldchange
  logmean        =   deseq$logmean
  
  ## fit mean log abundance 
  mean_fit         =    logmean_fit(logmean)
  logmean_param    =    mean_fit$logmean_param
  optim_logmean    =    mean_fit$components
  logmean_density  =    mean_fit$density_plot
  
  ## fit dispersion  
  dispers_fit        =    dispersion_fit(dispersion,logmean)
  dispersion_param   =    dispers_fit$param
  dispersion_confint =    dispers_fit$confint
  
  logfoldchange_param  =    logfoldchange_fit(logmean,logfoldchange)
  
  
  list(countdata_filted =   countdata_filt,
                 countdata_norm   =   countdata_norm,
                 control_count    =   control_count,
                 treat_count      =   treat_count,
                 dispersion       =   dispersion,
                 logfoldchange    =   logfoldchange,
                 logmean          =   logmean,
                 logmean_param    =   logmean_param,
                 optim_logmean    =   optim_logmean,
                 logmean_density  =   logmean_density,
                 dispersion_param =   dispersion_param,
                 dispersion_confint   = dispersion_confint,
                 logfoldchange_param  =   logfoldchange_param)
}
##################################################################
names(res)        =   names(countdata_list_obs)

otu_dataset_list  = list(filtered_otu    =  read_data(res, "countdata_filted"),
                         normalised_otu  =  read_data(res, "countdata_norm"),
                         control_otu     =  read_data(res, "control_count"),
                         treatment_otu   =  read_data(res, "treat_count"))

estimates_obs_data_list =  list(logfoldchange  =  read_data(res, "logfoldchange"),
                                dispersion     =  read_data(res, "dispersion"),
                                logmean        =  read_data(res, "logmean"))

fit_param_list  =  list(logmean_param      =  read_data(res, "logmean_param"),
                        optim_logmean        =  read_data(res, "optim_logmean"),
                        dispersion_param     =  read_data(res, "dispersion_param"),
                        dispersion_confint   =  read_data(res, "dispersion_confint"),
                        logfoldchange_param  =  read_data(res, "logfoldchange_param"))

saveRDS(otu_dataset_list, file = paste0(path,"otu_dataset_list.rds"))
saveRDS(estimates_obs_data_list, file = paste0(path,"estimates_obs_data_list.rds"))
saveRDS(fit_param_list, file = paste0(path,"fit_param_list.rds"))
########################################################################
samp_vec =  c(seq(130,200,20),200)
nsim = 100;  notu=1000
disp_scale = 0.3

for(j in 1:length(samp_vec)){
  nsamp  =  samp_vec[j]
  for(i in 1:7){
    nam  =   names(dispersion_param_list)[i] 
    dispersion_param      =   dispersion_param_list[[i]]
    logmean_param         =   logmean_param_list[[i]]
    logfoldchange_param   =   logfoldchange_param_list[[i]]
    
    ###### Data simulation 
    otu_data   =  countdata_sim_fun(logmean_param, logfoldchange_param, 
                                    dispersion_param,
                                    nsamp_per_group=nsamp,
                                    disp_scale = disp_scale,
                                    ncont = NULL,ntreat = NULL,
                                    notu, nsim, 
                                    maxlfc_iter = 10000,
                                    seed = 100)
    
    countdata_li           =   otu_data$countdata_list
    metadata_list          =   otu_data$metadata_list
    true_lmean_li          =   otu_data$logmean_list
    true_lfoldchange_li    =   otu_data$logfoldchange_list
    treat_countdata_list   =   otu_data$treat_countdata_list
    control_countdata_list =   otu_data$control_countdata_list
    
    countdata_list         =   foreach(k = 1:length(countdata_li)) %do% {
      countdata =  countdata_li[[k]]
      metadata  =  metadata_list[[k]]
      filter_low_otu    = filter_fun(countdata,metadata,abund_thresh=5,
                                     sample_thresh=3)
      countdata_filt    =  filter_low_otu$countdata_filt
      
      countdata_filt
    }
    names(countdata_list)  =   names(countdata_li)
 
    ################################################################################
    #Remove the filtered low abundance 
    true_lfoldchange_list   =   true_lmean_list  = list()
    
    for(k1 in 1:length(true_lfoldchange_li)){
      
      names_to_extract     =   rownames(countdata_list[[k1]])
      lfc      =   true_lfoldchange_li[[k1]]
      lmb      =   true_lmean_li[[k1]]
      
      lfc      =  lfc[names(lfc) %in% names_to_extract]
      lmb      =  lmb[names(lmb) %in% names_to_extract]
      
      stopifnot(names_to_extract ==  names(lfc))
      stopifnot(names_to_extract ==  names(lmb))
      true_lmean_list[[k1]]     =   lmb
      true_lfoldchange_list[[k1]] =   lfc
    }
    
    names(true_lfoldchange_list)  =   names(true_lfoldchange_li)
    names(true_lmean_list)        =   names(true_lmean_li)      
    
    
   p=unlist(lapply(true_lfoldchange_list,function(x) max(abs(x))))

    ###################################################
    dds1 =  deseq_fun_est(metadata_list[1:50],countdata_list[1:50],
                          num_cores=10, ref_name= "control")
    
    dds2 =  deseq_fun_est(metadata_list[51:100],countdata_list[51:100],
                          num_cores=10, ref_name= "control")
    
    dds  =  c(dds1,dds2)
    
    folder_path <- paste0(path,"deseq_results/",nsamp,"_samples/")
    
    #dds = readRDS(paste0(folder_path,nam,".rds"))
    # write a note explaining what you did
    # Check if the folder exists
    if (!file.exists(folder_path)) {
      # Create the folder if it doesn't exist
      dir.create(folder_path, recursive = TRUE)
      print(paste("Folder", folder_path, "created successfully."))
    } else {
      print(paste("Folder", folder_path, "already exists."))
    }
    
    saveRDS(dds, file = paste0(folder_path,nam,".rds"))
    
    deseq_est_list  =  read_data(dds, "deseq_estimate")
    
    my_gam  =  gam_fit(deseq_est_list,
                       true_lfoldchange_list,
                       true_lmean_list,
                       grid_len = 50,
                       alpha_level=0.1)
    
    combined_data  =  my_gam[["combined_data"]]
    power_estimate =  my_gam[["power_estimate"]]
    gam_model      =  my_gam[["fit_2d"]]
    
    folder_path1 <- paste0(path,"gam_model/",nsamp,"_samples/")
    
    # Check if the folder exists
    if (!file.exists(folder_path1)) {
      # Create the folder if it doesn't exist
      dir.create(folder_path1, recursive = TRUE)
      print(paste("Folder", folder_path1, "created successfully."))
    } else {
      print(paste("Folder", folder_path1, "already exists."))
    }
    
    saveRDS(gam_model, file = paste0(folder_path1,nam,".rds"))
    
    
    folder_path2 <- paste0(path,"contour_data/",nsamp,"_samples/")
    
    # # Check if the folder exists
    if (!file.exists(folder_path2)) {
      # Create the folder if it doesn't exist
      dir.create(folder_path2, recursive = TRUE)
      print(paste("Folder", folder_path2, "created successfully."))
    } else {
      print(paste("Folder", folder_path2, "already exists."))
    }
  
    saveRDS(combined_data, file = paste0(folder_path2,"combined_data_",nam,"_",".rds"))
    saveRDS(power_estimate, file = paste0(folder_path2,"power_estimate_",nam,"_",".rds"))
  }
  
}
################################################################################
#' check if my gam fit makes sense 
ddf = data.frame(lmean = unlist(true_lmean_list),abs= unlist(true_lfoldchange_list))
ggplot(ddf, aes(x = lmean, y= abs)) +
  geom_point()



for(j in 2:length(nsamp_vec)){
  
  j = 2
  ncores  =  10
  cl      =   makeCluster(ncores)
  registerDoParallel(cl)
  
  nsamp  =  nsamp_vec[[j]]# nsamp  =  50; 
    nsim = 100; disp_scale = 0.3
  notu   =  1000# mod  = list()

  #res = list()
for(i in 1:7){
  nam  =   names(metadata_list_obs)[i]

  dispersion_param      =   dispersion_param_list[[i]]
  logmean_param         =   logmean_param_list[[i]]
  logfoldchange_param   =   logfoldchange_param_list[[i]]
  ###### Data simulation 
  otu_data   =  countdata_sim_fun(logmean_param, logfoldchange_param, 
                                  dispersion_param,
                                  nsamp_per_group=nsamp,
                                  disp_scale = disp_scale,
                                  ncont = NULL,ntreat = NULL,
                                  notu, nsim, 
                                  seed = 100)
  
  
  countdata_list      =   otu_data$countdata_list
  metadata_list         =   otu_data$metadata_list
  true_lmean_list       =   otu_data$logmean_list
  true_lfoldchange_list    =   otu_data$logfoldchange_list
  treat_countdata_list    =   otu_data$treat_countdata_list
  control_countdata_list   =   otu_data$control_countdata_list
  
  ### compare variance of simulation and observed
  #compare_dataset(countdata_list[1:5], countdata_filt)
  
  #compare_dataset(countdata_sim_list = treat_countdata_list[1:5], 
  #                countdata_obs     =  treat_count_obs)
  
  #compare_dataset(control_countdata_list[1:5], 
   #               control_count_obs)
  ###################################################
  dds1 =  deseq_fun_est(metadata_list[1:50],countdata_list[1:50],
                        num_cores=10, ref_name= "control")

  dds2 =  deseq_fun_est(metadata_list[51:100],countdata_list[51:100],
                        num_cores=10, ref_name= "control")

  dds  =  c(dds1,dds2)

  folder_path <- paste0(path,"deseq_results/",nsamp,"_samples/")
  #dds = readRDS(paste0(folder_path,nam,".rds"))
  # write a note explaining what you did
  # Check if the folder exists
  if (!file.exists(folder_path)) {
    # Create the folder if it doesn't exist
    dir.create(folder_path, recursive = TRUE)
    print(paste("Folder", folder_path, "created successfully."))
  } else {
    print(paste("Folder", folder_path, "already exists."))
  }

  saveRDS(dds, file = paste0(folder_path,nam,".rds"))
  
  deseq_est_list  =  read_data(dds, "deseq_estimate")
  my_gam  =  gam_fit(deseq_est_list,
                     true_lfoldchange_list,
                     true_lmean_list,
                     grid_len = 50,
                     alpha_level=0.1)
  
  combined_data  =  my_gam[["combined_data"]]
  power_estimate =  my_gam[["power_estimate"]]
  
  
  folder_path1 <- paste0(path,"gam_model/",nsamp,"_samples/")
  dds = readRDS(paste0(folder_path,nam,".rds"))
  # write a note explaining what you did
  # Check if the folder exists
  if (!file.exists(folder_path1)) {
    # Create the folder if it doesn't exist
    dir.create(folder_path1, recursive = TRUE)
    print(paste("Folder", folder_path1, "created successfully."))
  } else {
    print(paste("Folder", folder_path1, "already exists."))
  }
  
  gam_model    =  my_gam[["fit_2d"]]
  saveRDS(gam_model, file = paste0(folder_path1,nam,".rds"))
   
  
  folder_path2 <- paste0(path,"contour_data/",notu,"_notu/")
  # 
  # # Check if the folder exists
  if (!file.exists(folder_path2)) {
    # Create the folder if it doesn't exist
    dir.create(folder_path2, recursive = TRUE)
    print(paste("Folder", folder_path2, "created successfully."))
  } else {
    print(paste("Folder", folder_path2, "already exists."))
  }

  saveRDS(combined_data, file = paste0(folder_path2,"combined_data_",nam,"_",".rds"))
  saveRDS(power_estimate, file = paste0(folder_path2,"power_estimate_",nam,"_",".rds"))
}

names(res) = names(filtered_otu_list)
saveRDS(res, file = paste0(path,"gam_model/gam_fit_",nsamp,"_samples",".rds"))

stopCluster(cl)
unregister_dopar() 
gc()
}
#nam = names(filtered_otu_list)[1]

# do it for some of the datasets. At least you have some of the 
#' 10 -12
#' 12 - 2 teach 
#' 2- 3 order for lava pizza pickup while working and then go and set up
#' 3:30- 4:30 send out the reminders again both by email and whatsapp
#' 

# res2 <- foreach(i = 1:7) %do% {
#   
#   source("reproducible/power/Load_Packages.R")
#   source("reproducible/power/fitting_fun.R")
#   source("reproducible/power/simulate_fun.R")
#   source("reproducible/power/utils.R")
#   
#   logmean   =  logmean_list[[i]]
#   logfoldchange  = logfoldchange_list[[i]]
#   
#   logfoldchange_fit(logmean,logfoldchange)
# }
# 
# names(res2)        =   names(countdata_list_obs)
# saveRDS(res2, file = paste0(path,"lfc_param_list.rds"))


# for(i in 1:7){
#   
# nam  =   names(metadata_list_obs)[i]
# countdata  =   countdata_list_obs[[i]]
# metadata   =   metadata_list_obs[[i]]
# ######################################################
# #filter low abundance taxa
# filter_low_otu    = filter_fun(countdata,metadata,abund_thresh=5,
#                              sample_thresh=3)
# 
# countdata_filt    =  filter_low_otu$countdata_filt
# control_count_obs =  filter_low_otu$control_count
# treat_count_obs   =  filter_low_otu$treat_count
# 
# #store filter_low_otu as one dataset as well as countdata_norm
# # and call it otu datasets list
# 
# ######################################################
# #run deseq for log fold change estimates 
# deseq          =   deseqfun(countdata_filt,metadata,ref_name="NT")
# countdata_norm =   deseq$normalised_count  
# dispersion     =   deseq$dispersion
# logfoldchange  =   deseq$logfoldchange
# logmean        =   deseq$logmean
# ######################################################
# ## fit mean log abundance 
# mean_fit         =    logmean_fit(logmean)
# logmean_param    =    mean_fit$logmean_param
# ######################################################
# ## fit log fold change  
# #saveRDS(logfoldchange_param, 
# #        file = paste0(path,"logfoldchange_param_",nam,".rds"))
# 
# logfoldchange_param = readRDS(paste0(path,"logfoldchange_param_",nam,".rds"))
# 
# #logfoldchange_param  =    logfoldchange_fit(logmean,logfoldchange)
#           par    =    logfoldchange_param$par
#           np     =    logfoldchange_param$np
#        sd_ord    =    logfoldchange_param$sd_ord
# ######################################################
# ## fit dispersion  
# dispers_fit       =    dispersion_fit(dispersion,logmean)
# dispersion_param  =    dispers_fit$param
# 
# ###### Data simulation 
# otu_data   =  countdata_sim_fun(logmean_param, logfoldchange_param, 
#                              dispersion_param,
#                              nsamp_per_group=nsamp,
#                              disp_scale = disp_scale,
#                              ncont = NULL,ntreat = NULL,
#                              notu, nsim, 
#                              seed = NULL)
# 
# 
#      countdata_list      =   otu_data$countdata_list
#    metadata_list         =   otu_data$metadata_list
#    true_lmean_list       =   otu_data$logmean_list
# true_lfoldchange_list    =   otu_data$logfoldchange_list
#  treat_countdata_list    =   otu_data$treat_countdata_list
# control_countdata_list   =   otu_data$control_countdata_list
# 
# ### compare variance of simulation and observed
# compare_dataset(countdata_list[1:5], countdata_filt)
# 
# compare_dataset(countdata_sim_list = treat_countdata_list[1:5], 
#                 countdata_obs     =  treat_count_obs)
# 
# compare_dataset(control_countdata_list[1:5], 
#                 control_count_obs)
# ###################################################
# dds1 =  deseq_fun_est(metadata_list[1:50],countdata_list[1:50], 
#                      num_cores=10, ref_name= "control")
# 
# dds2 =  deseq_fun_est(metadata_list[51:100],countdata_list[51:100], 
#                      num_cores=10, ref_name= "control")
# 
# ### You are forgetting to run dds2
# dds  =  c(dds1,dds2)
# 
# folder_path <- paste0(path,nsamp,"_samples/")
# 
# # Check if the folder exists
# if (!file.exists(folder_path)) {
#   # Create the folder if it doesn't exist
#   dir.create(folder_path, recursive = TRUE)
#   print(paste("Folder", folder_path, "created successfully."))
# } else {
#   print(paste("Folder", folder_path, "already exists."))
# }
# 
# saveRDS(dds, file = paste0(folder_path,"/",nam,".rds"))
# 
# dds = readRDS(paste0(folder_path,nam,".rds"))
# 
# deseq_est_list  =  read_data(dds, "deseq_estimate")
# deseq_estimate_list = dds
# 
# my_gam  =  gam_fit(deseq_est_list,
#                    true_lfoldchange_list,
#                    true_lmean_list,
#                    grid_len = 50,
#                    alpha_level=0.1)
# 
# combined_data  =  my_gam[["combined_data"]]
# power_estimate =  my_gam[["power_estimate"]]
# 
# saveRDS(combined_data, file = paste0(path,"combined_data_",nam,"_", nsamp,".rds"))
# saveRDS(power_estimate, file = paste0(path,"power_estimate_",nam,"_", nsamp,".rds"))
# }

cont_breaks = seq(0, 1, 0.1)
contour_plot_fun(combined_data, power_estimate, cont_breaks)


# library(ggrastr)
# cont_breaks    =  seq(0,1,0.1)
# plt = contour_plot_fun(combined_data,power_estimate,cont_breaks)



