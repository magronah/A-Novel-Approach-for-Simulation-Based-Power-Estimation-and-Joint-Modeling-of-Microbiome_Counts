library(RhpcBLASctl)
library(glmmTMB)
library(here)
path  = "real_data_analysis/atlass/results/"
ddf   =  readRDS(paste0(path,"dd_long.rds"))
############################################################
form      =   count ~1  +  offset(normalizer) +
                       (bmi_group*time|taxon) +
                       (1|nugget) + 
                       (1|subject:taxon) +
                       rr(taxon + 0 | subject:time,2) 
############################################################
par_ctrl <- glmmTMBControl(
  parallel = list(n = 10, autopar = TRUE),
  optCtrl  = list(eval.max=500, iter.max = 100)
)

gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")

options(glmmTMB_openmp_debug = TRUE)
blas_set_num_threads(1)

tt <- system.time( 
  fit  <-  glmmTMB(form, data  = ddf,
                   family  = poisson, 
                   ziformula  = ~1,
                   prior   = gprior,
                   REML    = FALSE,
                   control = par_ctrl
  )
)

saveRDS(fit, file=paste0(path,"mod_rr.rds"))
saveRDS(tt, file=paste0(path,"runtime_rr.rds"))


