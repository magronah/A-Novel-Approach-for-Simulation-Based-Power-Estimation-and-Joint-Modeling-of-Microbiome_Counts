library(here)
library(glmmTMB)
library(NBZIMM)
path   =  "real_data_analysis/pregnancy/results/"
#########################################################
meta_df     =   readRDS(paste0(path, "metadata.rds"))
ddd         =   readRDS(paste0(path, "countdata.rds"))

ttnb = system.time({
mod_nb    =   mms(y = ddd, fixed = ~pregnant*GA_Days  + offset(normalizer),
               random = ~ 1|Subect_ID,
               correlation = corAR1(form = ~ as.numeric(
                 factor(GA_Days, levels = unique(GA_Days)))| Subect_ID),
               data = meta_df,
               method = "nb",
               niter = 100)
})

ttznb = system.time({
  mod_znb    =   mms(y = ddd, fixed = ~pregnant*GA_Days  + offset(normalizer),
                     random = ~ 1|Subect_ID,
                     correlation = corAR1(form = ~ as.numeric(
                       factor(GA_Days, levels = unique(GA_Days)))| Subect_ID),
                     zi_fixed = ~1,
                     data = meta_df,
                     method = "zinb",
                     niter = 100)
})

# common_elements <- Reduce(intersect, list(names(mod_nb$fit), 
#                                           names(mod_znb$fit)))
# 
# saveRDS(common_elements,file = paste0(path,"include_taxa.rds"))

saveRDS(mod_nb, file = paste0(path,"mod_nb.rds"))
saveRDS(mod_znb, file = paste0(path,"mod_znb.rds"))
saveRDS(ttnb, file = paste0(path,"runtime_nb.rds"))
saveRDS(ttznb, file = paste0(path,"runtime_znb.rds"))
###########################################################
par_list = lapply(mod$fit, function(x){
  list(beta = fixef(x), 
       sd   = log(as.numeric(VarCorr(x)["(Intercept)", "StdDev"])),
       theta = log(x$theta))
})
# the uni
###########################################################
res   =  list()
for(i in 1:length(mod$fit)){
  y  =  ddd[,i]
  
  modA   =   glmmTMB(y ~ pregnant*GA_Days + (1 | Subect_ID) + offset(normalizer), 
                     data = meta_dd, 
                     family = nbinom2())
  
  modB1  =   glmmTMB(y ~ pregnant*GA_Days + (1 | Subect_ID) + offset(normalizer), 
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
                   AIC      =    aic,
                   AICc    =    aic +  correction,
                   mod    =    modA)
}

nbmm_aicc   <-   sum(sapply(res, `[[`, "AICc"))
saveRDS(nbmm_aicc, file=paste0(file_path,"nbmm_aicc.rds"))
saveRDS(res, file=paste0(file_path,"nbmm_res.rds"))




