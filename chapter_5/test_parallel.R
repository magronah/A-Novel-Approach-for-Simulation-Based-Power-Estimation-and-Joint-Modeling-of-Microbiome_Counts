#setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
source("reproducible/new_sim/Load_Packages.R")
source("reproducible/new_sim/func.R")
source("reproducible/new_sim/param.R")
source("reproducible/new_sim/initial_param0.R")
#source("load_glmmTMB.R")
library(glmmTMB)
library(RhpcBLASctl)
data_path = paste0(getwd(),"/reproducible/new_sim/")
###################################################################
ntaxa =  100; nsubj=50; beta; betadisp
theta_true =  c(
  ##log-sd of intercept, group,
  get_theta_logSD(n = 2, seed = seed),
  get_theta_corr(n  = 2,  seed=seed),
  get_theta_logSD(n = ntaxa, seed = seed),
  get_theta_corr(n  = ntaxa,  seed=seed)
)
form <- count ~ 1 + us(1 + group|taxon) + us(0 + taxon | subject)
#####################################################################
dd          =   met_data(ntaxa, nsubj)
true_pars   =   lst(beta, theta = theta_true, betadisp)
####################################################################  
dd$count    =  simulate_new(RHSForm(form, as.form = TRUE),
                           newdata   =  dd,
                           newparams =  true_pars,
                           family    =  nbinom2,
                           nsim      =  1,
                           seed   =  seed, 
                           se  = FALSE
                           )[[1]]
###############################################################################
n_list =  6; fit =  elapsed_time = list()
for(i in 1:length(n_list)){
  

par_ctrl <- glmmTMBControl(
  optCtrl = list(trace=1, eval.max = 1000,iter.max = 1000),
  parallel = list(n = 6,autopar = TRUE)
)


gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")

form2  =   update(form,.~.-us(0 + taxon | subject) +
                    rr(0 + taxon | subject, d=2))

options(glmmTMB_openmp_debug = TRUE)
blas_set_num_threads(1)

elapsed_time5 = system.time(
  fit[[i]]  <-  glmmTMB(form2, data = dd,
                   family  = nbinom2, 
                   prior   = gprior,
                   REML    = TRUE,
                   control = par_ctrl,
                   se      = FALSE
                   
  )
)

}

View(fit)
pp=(fit[[1]])
ranef(fit[[1]])
# elapsed_time
# elapsed_time2

#elapsed_time3
#saveRDS(elapsed_time, paste0(getwd(),"/reproducible/new_sim/elapsed_time1.rds"))
#print("with")


