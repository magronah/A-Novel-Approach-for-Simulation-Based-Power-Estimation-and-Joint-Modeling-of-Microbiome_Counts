library(glmmTMB)
library(NBZIMM)
##############################################################
path   =   paste0(getwd(),"/real_data/autism_data")
source(paste0(path,"/prep_data.R"))
##############################################################
ddd   =   t(countdata)
meta_dd$dummy  =  factor(1)
##############################################################
tt = system.time({
  mod    =   mms(y = ddd, fixed = ~group + offset(normalizer),
                 random = ~ 1 | dummy,
                 zi_fixed = ~1,
                 data = meta_dd, 
                 method = "zinb",
                 niter = 100)
})

###########################################################
file_path  =  paste0(path,"autism_data/results/")
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
  
  modA   =   glmmTMB(y ~ group + age + (1 | dummy) + offset(normalizer), 
                     data = meta_dd, 
                     ziformula = ~1,
                     family = nbinom2())
  
  
  modB1  =   glmmTMB(y ~ group + age + (1 | dummy) + offset(normalizer), 
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
#########################################################
##' run those that NBMM could not fit 
##' (I suppose due to small count values-check this!) 
##' in glmmTMB
diff = setdiff(colnames(ddd), taxa_name)
if(is.null(diff)){
  
  zinbmm_aicc   <-   sum(sapply(res1, `[[`, "AICc"))
  saveRDS(zinbmm_aicc, file=paste0(file_path,"zinbmm_aicc.rds"))
  saveRDS(res1, file=paste0(file_path,"zinbmm_res.rds"))
  
}else{
  
  res2 = list()
  for(i in 1:length(diff)){
    y  =  ddd[,diff[i]]
    
    mod   =   glmmTMB(y ~ group + age + (1 | dummy) + offset(normalizer), 
                      data = meta_dd, 
                      ziformula = ~1,
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
  zinbmm_aicc <-  sum(sapply(res1, `[[`, "AICc")) + sum(sapply(res2, `[[`, "AICc"))
  res    =   list(res1, res2)
  
  saveRDS(zinbmm_aicc, file=paste0(file_path,"zinbmm_aicc.rds"))
  saveRDS(res, file=paste0(file_path,"zinbmm_res.rds"))
}

length(names(zinbmm_mod$fit))
saveRDS(zinbmm_mod, file=paste0(path,"autism/data/zinbmm_mod",
                            nam,".rds"))
#saveRDS(mod_list, file=paste0(path,"autism/data/zinbmm_runtime",nam,".rds"))
###########################################################
par_zinbmm = lapply(zinbmm_mod$fit, function(x){
  list(beta = fixef(x), 
       sd   = log(as.numeric(VarCorr(x)["(Intercept)", "StdDev"])),
       theta = log(x$theta),
       betazi  = x$zi.fit[[1]] )
})
###########################################################
res   =  list()
n     =  names(zinbmm_mod$fit)

for(i in 1:length(n)){
  y  =  dd[,n[i]]
  
  modA   =   glmmTMB(y ~ group + (1 | dummy) + offset(normalizer), 
                     data = meta_dd, 
                     ziformula = ~1,
                     #prior = gprior,
                     family = nbinom2())
  
  
  modB1  =   glmmTMB(y ~ group + (1 | dummy) + offset(normalizer), 
                     data = meta_dd, 
                     ziformula = ~1,
                     #prior = gprior,
                     family = nbinom2(),
                     doFit  = FALSE)
  
  modB2 <- fitTMB(modB1, doOptim = FALSE)
  
  pars <- modB2$env$par
  f   =  par_zinbmm[[i]]
  pars[names(pars) == "beta"] <- f$beta
  pars[names(pars) == "theta"] <- f$sd
  pars[names(pars) == "betadisp"] <- f$theta
  pars[names(pars) == "betazi"] <- f$betazi
  
  
  ## drop RE parameters
  pars <- pars[names(pars) != "b"]
  
  par2 <- with(modA$obj$env, last.par.best[-random])

    
  num_params =  attr(logLik(modA), "df")  
  res[[i]]  = list(zinbmm_loglik    = modB2$fn(pars),
                   glmmTMB_loglik = modB2$fn(par2), 
                   params =  cbind(NBZIMM = pars, glmmTMB = par2),
                   AIC    = 2*(modB2$fn(pars)) +  2*(num_params),
                   mod    =   modA)
}
#####################################################################
gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")

diff = setdiff(colnames(ddd), (names(zinbmm_mod$fit)))

saveRDS(diff, file = paste0(path,"/taxa_exclude.rds"))
res2 = list()
for(i in 1:length(diff)){
  y  =  dd[,diff[i]]
  mod   =   glmmTMB(y ~ group + (1 | dummy) + offset(normalizer), 
                     data = meta_dd, 
                     ziformula = ~1,
                     prior = gprior,
                     family = nbinom2())
  
  num_params =  attr(logLik(mod), "df")  
  loglike = as.numeric(logLik(mod))
  
  res2[[i]] = list(mod = mod,
                   AIC    =  -2*(loglike) + 2*(num_params),
                   loglike = loglike)
}

####################################################################
zinbmm_aic   = c(AIC  =   sum(unlist((lapply(res, function(x){x[["AIC"]]})))) +
                 sum(unlist((lapply(res2, function(x){x[["AIC"]]})))),
               loglik  =   -sum(unlist((lapply(res, function(x){x[["nbmm_loglik"]]})))) +
                 sum(unlist((lapply(res2, function(x){x[["loglike"]]}))))
)
saveRDS(zinbmm_aic, file=paste0(path,"autism/data/zinbmm_aic",nam,".rds"))
