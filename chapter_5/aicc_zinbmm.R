library(tibble)
library(DESeq2)
library(Matrix)
library(glmmTMB)
library(huge)
library(here)
####################################################################
source("func2.R")
source("initial_param0.R")
####################################################################
path = paste0(nsubj,"_",ntaxa,"/")
path
####################################################################
data      =   readRDS(paste0(path,"otu_meta_list_withzi_taxa.rds"))
######################################################################
file_path  =   paste0("~/scratch/data/",nsubj,"_",ntaxa,"/zinbmm/")
file       =   list.files(file_path, full.names = TRUE)
######################################################################
aicc_list   =  Res   =  list()

for(j in 1:length(file)){
  countdata     =    data[[j]]$countdata
  meta_dd       =    data[[j]]$met_data
        mod     =    readRDS(file[j])
  
  meta_dd$dummy   =    factor(1)
  ######################################################################
  normalizer = normalizer_fun(countdata, meta_dd, 
                              ref_name="control",
                              ntaxa = ntaxa)
  meta_dd$normalizer  =  normalizer
  ######################################################################
  res     =    list()
  ###########################################################
  par_list = lapply(mod$fit, function(x){
    list(beta = fixef(x), 
         sd   = log(as.numeric(VarCorr(x)["(Intercept)", "StdDev"])),
         theta = log(x$theta),
         betazi  = x$zi.fit[[1]])
  })
  
  ###########################################################
  for(i in 1:ncol(countdata)){
    y  =  countdata[, i]
    
    modA   =   glmmTMB(y ~ group + (1 | dummy) + offset(normalizer), 
                       data = meta_dd, 
                       family = nbinom2())
    
    
    modB1  =   glmmTMB(y ~ group + (1 | dummy) + offset(normalizer), 
                       data   = meta_dd, 
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
    
    res[[i]]  = list(nbmm_LL     =    modB2$fn(pars),
                     glmmTMB_LL  =    modB2$fn(par2), 
                     params  =    cbind(NBZIMM = pars, glmmTMB = par2),
                     AIC   =    aic,
                     AICc    =    aic +  correction,
                     mod    =    modA)
  }
  #########################################################
  Res[[j]]        <-    res
  aicc_list[[j]]   <-   sum(sapply(res, `[[`, "AICc"))
}

names(Res)  =  names(aicc_list)  =  paste0("sim",1:length(file))
saveRDS(aicc_list, file=paste0(path,"zinbmm_aicc.rds"))
saveRDS(Res, file=paste0(path,"zinbmm_res.rds"))
