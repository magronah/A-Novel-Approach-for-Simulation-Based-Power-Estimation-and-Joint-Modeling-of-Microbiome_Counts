##############################################################
path   =   paste0(getwd(),"/real_data/soil_data")
source(paste0(path,"/prep_data.R"))
##############################################################
pp$result
res   =  list()
for(i in 1:ncol(dd)){
  y  =  dd[,i]
  
  modA   =   glmmTMB(y ~ group + offset(normalizer), 
                     data = meta_dd, 
                     family = nbinom2())
  
  
  modB1  =   glmmTMB(y ~ group + offset(normalizer), 
                     data = meta_dd, 
                     family = nbinom2(),
                     doFit  = FALSE)
  
  modB2 <- fitTMB(modB1, doOptim = FALSE)
  
  pars <- modB2$env$par
  pars[names(pars) == "beta"] <- c((pp$result$intercept[[i]])*log(2),
                                   log(2)*(pp$result$log2FoldChange[[i]])
  )
  pars[names(pars) == "betadisp"] <- log(1/pp$result$dispersion[[i]])
  
  
  ## drop RE parameters
  pars <- pars[names(pars) != "b"]
  
  par2 <-  modA$fit$par
  
  num_params =  attr(logLik(modA), "df")  
  res[[i]]  = list(deseq_loglik    = modB2$fn(pars),
                   glmmTMB_loglik = modB2$fn(par2), 
                   params =  cbind(deseq = pars, glmmTMB = par2),
                   AIC    = 2*(modB2$fn(pars)) +  2*(num_params),
                   mod    =   modA)
}


#####################################################################
gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")

diff = setdiff(colnames(dd), n)
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