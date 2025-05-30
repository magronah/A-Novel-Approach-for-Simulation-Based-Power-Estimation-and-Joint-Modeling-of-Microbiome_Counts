library(here)
library(NBZIMM)
library(glmmTMB)
source("pregnancy/prep_data.R")
#########################################################
mod    =   mms(y = dd, fixed = ~pregnant*GA_Days  + offset(normalizer),
               random = ~ 1|Subect_ID,
               zi_fixed = ~1,
               data = meta_df,
               method = "zinb",
               niter = 100)

file_path = "pregnancy/results/"
saveRDS(mod, file = paste0(file_path,"mod_zinbmm.rds"))
###########################################################
par_list = lapply(mod$fit, function(x){
  list(beta = fixef(x), 
       sd   = log(as.numeric(VarCorr(x)["(Intercept)", "StdDev"])),
       theta = log(x$theta),
       betazi  = x$zi.fit[[1]])
})
###########################################################
res   =  list()
for(i in 1:length(mod$fit)){
  y  =  ddd[,i]
  
  modA   =   glmmTMB(y ~ pregnant*GA_Days + (1 | Subect_ID) + offset(normalizer), 
                     data = meta_dd, 
                     ziformula = ~1,
                     family = nbinom2())
  
  modB1  =   glmmTMB(y ~ pregnant*GA_Days + (1 | Subect_ID) + offset(normalizer), 
                     data   = meta_dd, 
                     ziformula = ~1,
                     family = nbinom2(),
                     doFit  = FALSE)
  
  modB2 <- fitTMB(modB1, doOptim = FALSE)
  
  pars <- modB2$env$par
  f   =  par_list[[i]]
  pars[names(pars) == "beta"] <- f$beta
  pars[names(pars) == "theta"] <- f$sd
  pars[names(pars) == "betadisp"] <- f$theta
  pars[names(pars) == "betazi"] <- f$betazi
  
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

zinbmm_aicc   <-   sum(sapply(res, `[[`, "AICc"))
saveRDS(zinbmm_aicc, file=paste0(file_path,"zinbmm_aicc.rds"))
saveRDS(res, file=paste0(file_path, "zinbmm_res.rds"))
