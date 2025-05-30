source("reproducible/power/Load_Packages.R")
source("reproducible/power/fitting_fun.R")
source("reproducible/power/simulate_fun.R")
source("reproducible/power/utils.R")
source("reproducible/power/read_dataset.R")
######################################################
path  =  "reproducible/power/datasets2/"
nsim = 3
notu_vec=seq(1000,10000,1000)
disp_scale = 0.3; nsamp = 100

i=1
dispersion_param      =   dispersion_param_list[[i]]
logmean_param         =   logmean_param_list[[i]]
logfoldchange_param   =   logfoldchange_param_list[[i]]
#j=1
true_lfoldchange_li  = list()
for(j in 1:length(notu_vec)){
  notu = notu_vec[j]
  ###### Data simulation 
  otu_data   =  countdata_sim_fun(logmean_param, logfoldchange_param, 
                                  dispersion_param,
                                  nsamp_per_group=nsamp,
                                  disp_scale = disp_scale,
                                  ncont = NULL,ntreat = NULL,
                                  notu, nsim, 
                                  maxlfc_iter = 10000,
                                  seed = 100)
  true_lfoldchange_li[[j]]    =   mean(unlist(otu_data$logfoldchange_list))
}
