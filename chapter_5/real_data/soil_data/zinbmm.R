library(glmmTMB)
library(NBZIMM)
##############################################################
path   =   paste0(getwd(),"/real_data/soil_data/")
source(paste0(path,"prep_data.R"))
##############################################################
ddd   =   t(countdata)
###########################################################
tt = system.time({
  mod    =   mms(y = ddd, fixed = ~group + offset(normalizer),
                 random = ~ 1 | site,
                 zi_fixed = ~1,
                 data = meta_dd, 
                 method = "zinb",
                 niter = 100)
})
###########################################################
file_path  =  paste0(path,"results/")
if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}
###########################################################
saveRDS(tt, file=paste0(file_path,"zinbmm_runtime.rds"))
saveRDS(mod, file=paste0(file_path,"zinbmm_mod.rds"))

###########################################################
par_list = lapply(mod$fit, function(x){
  list(beta = fixef(x), 
       sd   = log(as.numeric(VarCorr(x)["(Intercept)", "StdDev"])),
       theta = log(x$theta),
       betazi  = x$zi.fit[[1]])
})
###########################################################
res1   =  list()
taxa_name  =  names(mod$fit)

for(i in 1:length(taxa_name)){
  y  =  dd[,taxa_name[i]]
  
  modA   =   glmmTMB(y ~ group + (1 | site) + offset(normalizer), 
                     data = meta_dd, 
                     ziformula = ~1,
                     family = nbinom2())
  
  
  modB1  =   glmmTMB(y ~ group + (1 | site) + offset(normalizer), 
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
  
  res1[[i]]  = list(zinbmm_LL    =    modB2$fn(pars),
                   glmmTMB_LL   =    modB2$fn(par2), 
                      params    =    cbind(NBZIMM = pars, glmmTMB = par2),
                          AIC   =    aic,
                        AICc    =    aic +  correction,
                         mod    =    modA)
}

zinbmm_aicc   <-   sum(sapply(res1, `[[`, "AICc"))
saveRDS(zinbmm_aicc, file=paste0(file_path,"zinbmm_aicc.rds"))
saveRDS(res1, file=paste0(file_path,"zinbmm_res.rds"))