setwd("/home/agronahm/Michael-n-Ben-Repo/")
source("reproducible/power/Load_Packages.R")
source("reproducible/power/simulate_fun.R")
source("reproducible/power/utils.R")
source("reproducible/power/fitting_fun.R")
source("reproducible/power/read_dataset.R")
path = "reproducible/power/datasets/"
###########################################################
## parameters for data simulation
nsim = 100; notu=1000; disp_scale = 0.3
nsamp_vec = seq(5,50,10)

i = 3
#index for a dateset
#########################################################
################Simulate data
              nam	=   names(dispersion_param_list)[i]
  dispersion_param	=   dispersion_param_list[[i]]
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


sim_logmean_list      =   read_data(dd, "logmean_list")
sim_logfoldchange_list       =   read_data(dd, "logfoldchange_list")
################################################################
##read the files containing the deseq results from
## running simulation in deseq
path_deseq  =   paste0(path,"deseq_results/")
deseq_est_list  = foreach(j = 1:length(nsamp_vec),.combine = "rbind") %do% {
      nsam	 =   nsamp_vec[j]
      ds         =   readRDS(paste0(path_deseq,nsam,"_samples/",nam,".rds"))
      deseq_list =   do.call(rbind, lapply(ds, function(x) x$deseq_estimate))
      deseq_list$sample_size = rep(nsam, nrow(deseq_list))
      deseq_list }
################################################################
## Run scam to fit gam
power_est       =   power_fun_ss(deseq_est_list, sim_logfoldchange_list,
                                   sim_logmean_list, nsamp_vec, alpha=0.1)

gam_fit_ss      =    power_est$gam_mod

summary(gam_fit_ss)
