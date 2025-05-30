library(RhpcBLASctl)
library(glmmTMB)
###########################################################
path   =   paste0(getwd(),"/real_data/soil_data/")
source(paste0(path,"prep_data.R"))
##############################################################
ddd     =    t(countdata)
param   =    coef(pp$object)
result  =    pp$result
res     =    list()
#View(param) i =1

for(i in 1:ncol(ddd)){
  y      =   ddd[,i]
  modA   =   glmmTMB(y ~ group + site + offset(normalizer), 
                     data = meta_dd, 
                     family = nbinom2())
  
  
  modB1  =   glmmTMB(y ~ group + site + offset(normalizer), 
                     data   = meta_dd, 
                     family = nbinom2(),
                     doFit  = FALSE)
  
  modB2 <- fitTMB(modB1, doOptim = FALSE)
  
  pars <- modB2$env$par
  
  pars[names(pars) == "beta"]    <-  log(2)*as.numeric(param[i,])
  pars[names(pars) == "betadisp"] <- log(1/result$dispersion[[i]])
  
  par2     <-   modA$fit$par
  
  num_params  =   attr(logLik(modA), "df")  
  aic         =   2*(modB2$fn(pars)) +  2*(num_params)
  correction  =   2*num_params*(num_params+1)/(length(y) - num_params -  1)
  
  res[[i]]  = list(deseq_LL     =    modB2$fn(pars),
                   glmmTMB_LL  =    modB2$fn(par2), 
                   params  =    cbind(deseq = pars, glmmTMB = par2),
                   AIC     =    aic,
                   AICc    =    aic +  correction,
                   mod     =    modA)
}

####################################################################
file_path  =  paste0(path,"results/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}

deseq_aicc   <-   sum(sapply(res, `[[`, "AICc"))
saveRDS(deseq_aicc, file=paste0(file_path,"deseq_aicc.rds"))
saveRDS(res, file=paste0(file_path,"deseq_res.rds"))


# fit1  <-  glmmTMB(form, data = df,
#                  family  = nbinom2, 
#                  prior   = gprior,
#                  REML    = FALSE,
#                  doFit  = FALSE)
# modB2 <- fitTMB(fit1, doOptim = FALSE)
# par2 <- with(fit$obj$env, last.par.best[-random])
# modB2$fn(pars2)

