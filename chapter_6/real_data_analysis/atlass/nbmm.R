library(here)
library(NBZIMM)
library(glmmTMB)
source("real_data_analysis/atlass/prep_data.R")
#########################################################
ttnb = system.time({
  modnb    =   mms(y = ddd, fixed = ~bmi_group*time  + offset(normalizer),
                   random = ~ 1|subject,
                   correlation = corAR1(form = ~ as.numeric(
                     factor(time, levels = unique(time)))| subject),
                   data = meta_dd,
                   method = "nb",
                   niter = 100)
})

ttznb = system.time({
  modznb    =   mms(y = ddd, fixed = ~bmi_group*time  + offset(normalizer),
                    random = ~ 1|subject,
                    correlation = corAR1(form = ~ as.numeric(
                      factor(time, levels = unique(time)))| subject),
                    zi_fixed = ~1,
                    data = meta_dd,
                    method = "zinb",
                    niter = 100)
})

common_elements <- Reduce(intersect, list(names(modnb$fit), 
                                          names(modznb$fit)))

#saveRDS(common_elements,file = "real_data_analysis/atlass/include_taxa.rds")
file_path  =  "real_data_analysis/atlass/results/"
saveRDS(modnb, file = paste0(file_path,"mod_nb.rds"))
saveRDS(modznb, file = paste0(file_path,"mod_znb.rds"))
saveRDS(ttnb, file = paste0(file_path,"runtime_nb.rds"))
saveRDS(ttznb, file = paste0(file_path,"runtime_znb.rds"))
#########################################################
par_list = lapply(mod$fit, function(x){
  list(beta = fixef(x), 
       sd   = log(as.numeric(VarCorr(x)["(Intercept)", "StdDev"])),
       theta = log(x$theta))
})
###########################################################
res   =  list()

for(i in 1:length(mod$fit)){
  y  =  ddd[,i]
  
  
  modA   =   glmmTMB(y ~ bmi_group*time + (1 | subject) + offset(normalizer), 
                     data = meta_dd, 
                     family = nbinom2())
  
  modB1  =   glmmTMB(y ~ bmi_group*time + (1 | subject) + offset(normalizer), 
                     data   = meta_dd, 
                     family = nbinom2(),
                     doFit  = FALSE)
  
  modB2 <- fitTMB(modB1, doOptim = FALSE)
  
  pars <- modB2$env$par
  f   =  par_list[[i]]
  pars[names(pars) == "beta"] <- f$beta
  pars[names(pars) == "theta"] <- f$sd
  pars[names(pars) == "betadisp"] <- f$theta
  
  ## drop RE parameters
  pars <- pars[names(pars) != "b"]
  
  par2 <- with(modA$obj$env, last.par.best[-random])
  
  num_params  =  attr(logLik(modA), "df")  
  aic      =   2*(modB2$fn(pars)) +  2*(num_params)
  correction  =   2*num_params*(num_params+1)/(length(y) - num_params -  1)
  
  res[[i]]  = list(zinbmm_LL    =    modB2$fn(pars),
                   glmmTMB_LL   =    modB2$fn(par2), 
                   params    =    cbind(NBZIMM = pars, glmmTMB = par2),
                   AIC   =    aic,
                   AICc    =    aic +  correction,
                   mod    =    modA)
}

nbmm_aicc   <-   sum(sapply(res, `[[`, "AICc"))
saveRDS(nbmm_aicc, file=paste0(file_path,"nbmm_aicc.rds"))
saveRDS(res, file=paste0(file_path,"nbmm_res.rds"))
