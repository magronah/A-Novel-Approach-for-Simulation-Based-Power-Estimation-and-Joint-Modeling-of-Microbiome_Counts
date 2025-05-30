path = "reproducible/power/datasets/"
###########################Read otu dataset###########################
   otu_dataset_list     =    readRDS(paste0(path,"otu_dataset_list.rds"))
   filtered_otu_list    =    otu_dataset_list[["filtered_otu"]]
  normalised_otu_list   =    otu_dataset_list[["normalised_otu"]]
     control_otu_list   =    otu_dataset_list[["control_otu"]]
   treatment_otu_list   =    otu_dataset_list[["treatment_otu"]]

   nam = "PRJNA453621"
  dim(control_otu_list[[nam]])
  dim(treatment_otu_list[[nam]])
  
  dim(control_otu_list[[nam]]) +   dim(treatment_otu_list[[nam]])
# rowSums(filtered_otu_list[[1]])
##########Read logfoldchange, logmean and dispersion estimates#########
estimates_obs_data_list  =   readRDS( paste0(path,"estimates_obs_data_list.rds"))
logfoldchange_list       =   estimates_obs_data_list[["logfoldchange"]]
logmean_list             =   estimates_obs_data_list[["logmean"]]
dispersion_list          =   estimates_obs_data_list[["dispersion"]]
########################Read parameters################################
fit_param_list           =   readRDS(paste0(path,"fit_param_list.rds"))
logmean_param_list       =   fit_param_list[["logmean_param"]]
dispersion_param_list    =   fit_param_list[["dispersion_param"]]
logfoldchange_param_list =   fit_param_list[["logfoldchange_param"]]

###############Deseq's estimates for simulated data####################
#nsamp_vec     =    seq(50,300,50)

# nsamp_vec     =   seq(10,100,10)
# contour_data  =   power_estimate =  list()
# df1   =   df2 =   list()
# 
# folder_path  =   paste0(path,"contour_data/")
# 
# for(i in 1:7){
#   nam     =   names(filtered_otu_list)[i]
# 
#   for(j in 1:length(nsamp_vec)){
#    nsam    =   nsamp_vec[j]
#    dd1      =   readRDS(paste0(folder_path,nsam,"_samples/combined_data_",
#                        nam, "_.rds"))
#    dd2      =   readRDS(paste0(folder_path,nsam,"_samples/power_estimate_",
#                                nam, "_.rds"))
# 
#    dd1$sample_size = rep(nsam, nrow(dd1))
#    dd2$sample_size = rep(nsam, nrow(dd2))
# 
#    df1[[j]]   =  dd1
#    df2[[j]]   =  dd2
#   }
# 
#   names(df1)  =  names(df2)  =  paste0(nsamp_vec, "samples")
#   contour_data[[i]]    =  df1
#   power_estimate[[i]]  =  df2
# }
# 
# names(contour_data)  =  names(power_estimate)  =  names(filtered_otu_list)

# 
# 
# #####################Gam Fit Models##############################
# nsamp_vec    =  c(20,seq(50,200, 50))
# gam_model_samp    =  gm  =  list()
# 
# folder_path_mod  =   paste0(path,"gam_model/")
# 
# for(i in 1:7){
#   nam     =   names(filtered_otu_list)[i]    
#   
#   for(j in 1:length(nsamp_vec)){
#     nsam     =   nsamp_vec[j]
#     gm[[j]] =   readRDS(paste0(folder_path_mod,nsam,"_samples/", 
#                                 nam, ".rds"))
# 
#   }
#   
#   names(gm)  = paste0(nsamp_vec, "samples") 
#   gam_model_samp[[i]]    =  gm
# }     
# names(gam_model_samp)  =  names(filtered_otu_list)
# (gam_model_samp[[1]][[1]]$model)
# 
# #####################Gam Fit Models##############################
# notu_vec    = seq(500,1500, 500)
# gam_model_otu    =  gm  =  list()
# 
# folder_path_mod  =   paste0(path,"gam_model/")
# 
# for(i in 1:7){
#   nam     =   names(filtered_otu_list)[i]    
#   
#   for(j in 1:length(notu_vec)){
#     nsam     =   notu_vec[j]
#     gm[[j]] =   readRDS(paste0(folder_path_mod,nsam,"_notu/", 
#                                nam, ".rds"))
#     
#   }
#   
#   names(gm)  = paste0(notu_vec, "otu") 
#   gam_model_otu[[i]]    =  gm
# }     
# names(gam_model_otu)  =  names(filtered_otu_list)
# 

# path_deseq  =   paste0(path,"deseq_results/")
# deseq_sample_size  =  list()
# nsamp_vec  =  seq(10,100,20)
# 
# for(i in 1:7){
#   nam     =   names(filtered_otu_list)[i]
#   dds = foreach(j = 1:length(nsamp_vec),.combine = "rbind") %do% {
#     nsam       =   nsamp_vec[j]
#     dd         =   readRDS(paste0(path_deseq,nsam,"_samples/",nam,".rds"))
#     deseq_list =   do.call(rbind, lapply(dd, function(x) x$deseq_estimate))
#     deseq_list$sample_size = rep(nsam, nrow(deseq_list))
#     deseq_list
#   }
#   deseq_sample_size[[i]]    =  dds
# }
# names(deseq_sample_size)  =  names(filtered_otu_list)

#lapply(countdata_list_obs, dim)

# library(tibble)
# df1 <- enframe(logfoldchange_list) %>% 
#          unnest(value)  %>% 
#         setNames(c("name","logfoldchange")) %>%
#         as.data.frame()
# 
# df2 <- enframe(logmean_list) %>% 
#        unnest(value) %>% 
#        setNames(c("name","logmean")) %>%
#        as.data.frame()
# stopifnot(df1$name == df2$name)
# 
# df  = data.frame(df2,variance=  sqrt(df1$logfoldchange^2))
# 
# oka_col = c(
#   "#E23D28", 
#   "#0000FF",
#   "#E69F00" # Orange
# )
# ggplot(df, aes(logmean, variance)) +
#   geom_point(alpha=0.1)  + 
#   geom_smooth(aes(color="smooth line")) +
#   geom_smooth(method = "lm", formula = y ~poly(x,2), aes(color="quadratic"),lwd=1) +
#   geom_smooth(method = "lm", formula = y ~x, aes(color="linear"),lwd=1) +
#   xlab(TeX("$\\log_2$(mean count)")) + ylab(TeX("|$\\log_2$(fold change)|")) +
#   scale_color_manual(name = "Type",
#                      values = oka_col) +
#   facet_wrap(~name,scale="free")+
#   theme_bw(base_size = 10)+
#   labs(color = " ")#, text = element_text(size = 20))



# scale_location <- function(logmean,logfoldchange){
#   
#   dat <- data.frame(logmean=logmean, 
#                     variance=sqrt(logfoldchange^2))
#   
#   ggplot(dat, aes(logmean, variance)) +
#     geom_point(alpha=0.1)  + 
#     geom_smooth(aes(color="smooth line")) +
#     geom_smooth(method = "lm", formula = y ~poly(x,2), aes(color="quadratic"),lwd=1) +
#     geom_smooth(method = "lm", formula = y ~x, aes(color="linear"),lwd=1) +
#     xlab("logmean") + ylab(expression(sqrt(logfoldchange^2))) +
#     theme_bw(base_size = 10)#, text = element_text(size = 20))
# }






