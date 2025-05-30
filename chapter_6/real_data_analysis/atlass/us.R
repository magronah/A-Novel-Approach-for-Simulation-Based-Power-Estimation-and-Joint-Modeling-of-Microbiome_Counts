library(RhpcBLASctl)
library(glmmTMB)
library(here)
path  = "real_data_analysis/atlass/results/"
ddf   =  readRDS(paste0(path,"dd_long.rds"))
############################################################
form      =   count ~1  +  offset(normalizer) +
                        (bmi_group*time|taxon) +
                        (1|nugget) + 
                        (1|subject:taxon) 
############################################################
par_ctrl <- glmmTMBControl(
  parallel = list(n = 10, autopar = TRUE)
)

gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")

options(glmmTMB_openmp_debug = TRUE)
blas_set_num_threads(1)


tt = system.time(
  fit  <-  glmmTMB(form, data  = ddf,
                   family  = poisson, 
                   ziformula  = ~1,# + (1|taxon),
                   prior   = gprior,
                   REML    = FALSE,
                   control = par_ctrl
  )
)

saveRDS(fit, file=paste0(path,"mod_us.rds"))
saveRDS(tt, file=paste0(path,"runtime_us.rds"))
