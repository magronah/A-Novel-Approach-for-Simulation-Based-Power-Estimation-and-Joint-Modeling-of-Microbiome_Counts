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
                 data = meta_dd, 
                 method = "nb",
                 niter = 100)
})
################################################################
file_path  =  paste0(path,"results/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}
###########################################################
saveRDS(tt, file=paste0(file_path,"nbmm_runtime.rds"))
saveRDS(mod, file=paste0(file_path,"nbmm_mod.rds"))

###########################################################
par_list = lapply(mod$fit, function(x){
  list(beta = fixef(x), 
       sd   = log(as.numeric(VarCorr(x)["(Intercept)", "StdDev"])),
       theta = log(x$theta))
})
###########################################################
     res1   =  list()
taxa_name  =  names(mod$fit)

for(i in 1:length(taxa_name)){
  y  =  dd[,i]
  
  modA   =   glmmTMB(y ~ group + (1 | site) + offset(normalizer), 
                     data = meta_dd, 
                     family = nbinom2())
  
  
  modB1  =   glmmTMB(y ~ group + (1 | site) + offset(normalizer), 
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

  res1[[i]]  = list(nbmm_LL     =    modB2$fn(pars),
                   glmmTMB_LL  =    modB2$fn(par2), 
                       params  =    cbind(NBZIMM = pars, glmmTMB = par2),
                         AIC   =    aic,
                       AICc    =    aic +  correction,
                        mod    =    modA)
}
#########################################################
nbmm_aicc   <-   sum(sapply(res1, `[[`, "AICc"))
saveRDS(nbmm_aicc, file=paste0(file_path,"nbmm_aicc.rds"))
saveRDS(res1, file=paste0(file_path,"nbmm_res.rds"))

##' run those that NBMM could not fit 
##' (I suppose due to small count values-check this!) 
##' in glmmTMB


diff = setdiff(colnames(dd), taxa_name)
if(is.null(diff)){
  
  nbmm_aicc   <-   sum(sapply(res1, `[[`, "AICc"))
  saveRDS(nbmm_aicc, file=paste0(file_path,"nbmm_aicc.rds"))
  saveRDS(res1, file=paste0(file_path,"nbmm_res.rds"))
  
}else{
  res2 = list()
  
  for(i in 1:length(diff)){
    y  =  dd[,diff[i]]
    
    mod   =   glmmTMB(y ~ group + (1 | site) + offset(normalizer), 
                      data = meta_dd, 
                      family = nbinom2())
    
    num_params  =   attr(logLik(mod), "df")  
    LL    =   as.numeric(logLik(mod))
    aic      =   -2*LL +  2*num_params
    correction  =    2*num_params*(num_params+1)/(length(y) - num_params -  1)
    
    
    res2[[i]] = list(mod    =   mod,
                     AIC    =   aic,
                     AICc   =   aic + correction,
                     LL     =   LL)
  }
  ####################################################################
  nbmm_aicc <-  sum(sapply(res, `[[`, "AICc")) + sum(sapply(res2, `[[`, "AICc"))
  saveRDS(nbmm_aicc, file=paste0(file_path,"nbmm_aicc.rds"))
  res  = list(res1, res2)
  saveRDS(res, file=paste0(file_path,"nbmm_res.rds"))
}