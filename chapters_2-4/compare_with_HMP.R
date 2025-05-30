library(HMP)
library(DirichletMultinomial)
library(metaSPARSim)
library(ggokabeito)
######################################################
source("reproducible/power/read_dataset.R")
source("reproducible/power/Load_Packages.R")
source("reproducible/power/simulate_fun.R")
source("reproducible/power/utils.R")
path = "reproducible/power/datasets/"
####################################################################
countdata_list_obs  =  readRDS(file = paste0(path,"data.rds"))
metadata_list_obs   =  readRDS(file = paste0(path,"metadata.rds"))
pp=(countdata_list_obs[[1]])

rownames(pp) = paste0("taxon",nrow(pp):1)
colnames(pp) = paste0("subject",1:ncol(pp))
View(t(pp))
######################################################
ncores  =  10
cl      =   makeCluster(ncores)
registerDoParallel(cl)
hmp = metasparsim = RR =  list()
res <-  foreach(n = 1:7, .errorhandling = "pass") %dopar% {  
  library(HMP)
  library(DirichletMultinomial)
  library(metaSPARSim)
  
  source("reproducible/power/read_dataset.R")
  source("reproducible/power/Load_Packages.R")
  source("reproducible/power/simulate_fun.R")
  source("reproducible/power/utils.R")
  
  nam    =   names(metadata_list_obs)[n]
  countdata   =   filtered_otu_list[[n]]
  metadata    =   metadata_list_obs[[n]]
  
  logmean_param = logmean_param_list[[n]]
  logfoldchange_param = logfoldchange_param_list[[n]]
  dispersion_param  =  dispersion_param_list[[n]]
  ######################################################
  #HMP needs the data sets to have sample by taxa format.
  countdata_hmp =    t(countdata)
  shape         =    dirmult(countdata_hmp)$gamma
  size          =    rowSums(countdata_hmp)
  sim_hmp       =    Dirichlet.multinomial(size, shape)
  sim_hmp       =    as.data.frame(t(sim_hmp))
  hmp[[n]]      =    sim_hmp
  # ######################################################
  # #simulate from metaSPARSim
  # #metaSPARSim needs the data sets to have taxa by sample format.
  normalised_data =    normalised_otu_list[[n]]
  control_otu     =    control_otu_list[[n]]
  treatment_otu   =    treatment_otu_list[[n]]
  nt_index        =    which(colnames(normalised_data) %in%
                               colnames(control_otu))
  asd_index       =    which(colnames(normalised_data) %in%
                               colnames(treatment_otu))
  stopifnot(colnames(control_otu) ==
              colnames(normalised_data)[nt_index])
  stopifnot(colnames(treatment_otu) ==
              colnames(normalised_data)[asd_index])

  grp_index      =    list(NT = nt_index, ASD = asd_index)
  params         =    estimate_parameter_from_data(countdata, normalised_data, grp_index)
  names(params)  =    names(grp_index)
  sim_data       =    metaSPARSim(params)

  countdata_metaSPARSim = sim_data$counts
  metasparsim[[n]] =   countdata_metaSPARSim
  #################################################################
  #########my simulating from our method 
  nsim = 1;  notu = nrow(countdata)
  ncont = length(nt_index); ntreat = length(asd_index)

  ########simulate from scaled dispersion
  disp_scale = 0.3
  otu_data1   =  countdata_sim_fun(logmean_param, logfoldchange_param, 
                                   dispersion_param,
                                   nsamp_per_group=NULL,
                                   disp_scale = disp_scale,
                                   ncont = ncont,ntreat = ntreat,
                                   notu, nsim, 
                                   seed = 100)
  countdata_scaled =  otu_data1$countdata_list[[1]]
  ##############################################
  #simulate from unscaled dispersion
  # disp_scale  =  1
  # otu_data2   =  countdata_sim_fun(logmean_param, logfoldchange_param, 
  #                                  dispersion_param,
  #                                  nsamp_per_group=NULL,
  #                                  disp_scale = disp_scale,
  #                                  ncont = ncont,ntreat = ntreat,
  #                                  notu, nsim, 
  #                                  seed = 100)
  # countdata_unscaled  =  otu_data2$countdata_list[[1]]

  countdata_sim_list  =  list(metaSPARSim   =   countdata_metaSPARSim,
                              scaled        =   countdata_scaled,
                              #RR            =   RR,
                              #unscaled      =   countdata_unscaled,
                              HMP           =   sim_hmp)
  
  countdata_sim_list
}

stopCluster(cl)
unregister_dopar() 
names(res)  =  names(metadata_list_obs) 

path  =  "reproducible/power/datasets/sample_size_cal/datasets/"
saveRDS(res, file = paste0(path,"countdata_sim_compare.rds"))





# pp2=(var_plt1[[1]]|var_plt1[[2]]|var_plt1[[3]])/(var_plt1[[4]]|var_plt1[[6]]|var_plt1[[7]])  +  plot_layout(guides = "collect")  
# pp1/pp2

# we generated the contour plots with ...
# aroun 9:00 
(mean_plt[[1]]|mean_plt[[2]]|mean_plt[[3]])/(mean_plt[[4]]|mean_plt[[6]]|mean_plt[[7]])  +  plot_layout(guides = "collect")  
(var_plt[[1]]|var_plt[[2]]|var_plt[[3]])/(var_plt[[4]]|var_plt[[6]]|var_plt[[7]])  +  plot_layout(guides = "collect")  


(mean_plt[[1]]|mean_plt[[2]]|mean_plt[[3]])  +  plot_layout(guides = "collect")  
(var_plt[[1]]|var_plt[[2]]|var_plt[[3]]) +  plot_layout(guides = "collect")  
pp1/pp2


pp11 = (mean_plt[[4]]|mean_plt[[6]]|mean_plt[[7]])  +  plot_layout(guides = "collect")  
pp22 = (var_plt[[4]]|var_plt[[6]]|var_plt[[7]]) +  plot_layout(guides = "collect")  
pp11/pp22


lmean_component_list= readRDS(file = paste0(path,"lmean_component_list.rds"))
lmean_density_plt_list = readRDS(file = paste0(path,"lmean_density_plt_list.rds"))
logmean_param_list=  readRDS(file = paste0(path,"logmean_param_list.rds"))
dispersion_param_list =  readRDS(file = paste0(path,"dispersion_param_list.rds"))
dispersion_confint_list = readRDS(file = paste0(path,"dispersion_confint_list.rds"))
lfc_lm_plot_list  = readRDS(file = paste0(path,"lfc_lm_plot_list.rds"))
optimal_model_list = readRDS(file = paste0(path,"optimal_model_list.rds"))


lfc_plt_list  =  readRDS(file = paste0(path,"lfc_plt_list.rds"))
lfc_aic_list  = readRDS(file = paste0(path,"lfc_aic_list.rds"))
logfoldchange_param_list = readRDS(file = paste0(path,"logfoldchange_param_list.rds"))
countdata_filt_list = readRDS(file = paste0(path,"countdata_filt_list.rds"))
######################################################
source("R_files/Load_Packages.R")
source("R_files/debug_single_example_fun.R")
source("R_files/debug_variance_functions.R")
source("R_files/debug_initial_conditions.R")
######################################################
path     <- "Datasets/"
plt.path <- "Figures/"
######################################################
countdata_list = readRDS(file = "Datasets/data.rds")
metadata_list = readRDS(file = "Datasets/metadata.rds")
logmean_param_list=  readRDS(file = paste0(path,"logmean_param_list.rds"))
dispersion_param_list =  readRDS(file = paste0(path,"dispersion_param_list.rds"))
logfoldchange_param_list = readRDS(file = paste0(path,"logfoldchange_param_list.rds"))
######################################################
n=1; nam    =   names(metadata_list)[n]
countdata =   countdata_list[[n]]
metadata  =   metadata_list[[n]]
logmean_param = logmean_param_list[[n]]
logfoldchange_param = logfoldchange_param_list[[n]]
dispersion_param  =  dispersion_param_list[[n]]
######################################################
v = filter(countdata, metadata,filt_thresh=5)
names(v) 
countdata_obs = v$countdata
#countdat = data.frame(t(v$countdata))
class(countdata_obs)


nsim = 5; nsamp=ncol(countdata_obs)/2; notu=nrow(countdata_obs)
##########################################################
d1 = simulate_fun(logmean_param,logfoldchange_param,
                 dispersion_param,
                 nsamples_per_group=nsamp,
                 notu=notu,nsim=nsim,
                 scale_dispersion=0.3,
                 num_cores = 11,
                 seed = 100)

d1 = d1$countdata_list[1:2]

d2 = simulate_fun(logmean_param,logfoldchange_param,
                  dispersion_param,
                  nsamples_per_group=nsamp,
                  notu=notu,nsim=nsim,
                  scale_dispersion=1,
                  num_cores = 11,
                  seed = 100)

d2 = d2$countdata_list[1:2]
#########################################################



compare_dataset(d1,countdata_obs, data_name= "XYZ")

class(d1)
length(d1)
compare_dataset(countdata_obs,d2)
#########################################################
dlist <- list.append(contdat_list,obs_filt=countdata_filt)

dlist = list.append(contdat_list1,contdat_list,obs_filt=countdata_filt)
dlist <- list.append(contdat_list,obs_filt=countdata_filt)
vars <-(dlist |> purrr::map(~apply(., 1, stats::var)) |>
          purrr::map_dfr(~ tibble(var = .),.id = "type")) 

vars_obs <- vars[vars$type == "obs_filt",]
vars_others <- vars[vars$type != "obs_filt",]

vars_scaled <- vars[vars$type == "scaled",]
vars_unscaled <- vars[vars$type == "unscaled",]
#vars_HMP <- vars[vars$type == "HMP",]

p = ggplot(vars_others, aes(x = log(var), colour = type)) + 
  geom_density(lwd=1.5) +
  geom_density(data=vars_obs, aes(x = log(var)), linetype = "solid",lwd=3) +
  xlab(TeX("$\\log_{10}$(variance of taxa)"))  
p


dlist = list(obs_filt=countdata_obs, scaled=d1,unscaled=d2,HMP=d3)
mean_vec <-(dlist |> purrr::map(~apply(., 1, mean)) |>
          purrr::map_dfr(~ tibble(mean = .),.id = "type")) 

mean_obs <- mean_vec[mean_vec$type == "obs_filt",]
mean_others <- mean_vec[mean_vec$type != "obs_filt",]

p2 = ggplot(mean_others, aes(x = log(mean), colour = type)) + 
  geom_density(lwd=1.5) +
  geom_density(data=mean_obs, aes(x = log(mean)), linetype = "solid",lwd=3) +
  xlab(TeX("$\\log_{10}$(variance of taxa)"))  
p2

  
  
  # otu_table <- data.frame(t(countdata)) %>%
  #   rownames_to_column(var = "Samples")
  # 
  # combined_otu_table=left_join(metadata,otu_table, by = "Samples")  
  # 
  # otu_table_control <- combined_otu_table %>%
  #   dplyr::filter(Groups == "NT") %>%
  #   dplyr::select(-c(Groups)) %>%
  #   column_to_rownames("Samples")
  # 
  # otu_table_treatment <- combined_otu_table %>%
  #   dplyr::filter(Groups == "ASD") %>%
  #   dplyr::select(-c(Groups)) %>%
  #   column_to_rownames("Samples")
  # 
  # size_cont = rowSums(otu_table_control)
  # size_treat = rowSums(otu_table_treatment)
  