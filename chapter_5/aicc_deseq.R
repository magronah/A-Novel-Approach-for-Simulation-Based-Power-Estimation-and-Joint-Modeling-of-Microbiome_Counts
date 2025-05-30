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
file_path  =   paste0("~/scratch/data/",nsubj,"_",ntaxa,"/deseq/")
file       =   list.files(file_path, full.names = TRUE)
######################################################################
aicc_list   =  Res   =  list()

for(j in 1:length(file)){
  countdata  =    data[[j]]$countdata
  met_data   =    data[[j]]$met_data
  result     =    readRDS(file[j])
  ######################################################################
  normalizer = normalizer_fun(countdata, met_data, 
                              ref_name="control",
                              ntaxa = ntaxa)
  met_data$normalizer  =  normalizer
  ######################################################################
  res     =    list()
  for(i in 1:ncol(countdata)){
    y      =   countdata[,i]
    modA   =   glmmTMB(y ~ group + offset(normalizer), 
                       data = met_data, 
                       family = nbinom2())
    
    
    modB1  =   glmmTMB(y ~ group + offset(normalizer), 
                       data   = met_data, 
                       family = nbinom2(),
                       doFit  = FALSE)
    
    modB2 <- fitTMB(modB1, doOptim = FALSE)
    
    pars <- modB2$env$par
    
    
    pars[names(pars) == "beta"]    <-  log(2)*c(result$intercept[i],
                                                result$log2FoldChange[i])
    
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
  
  Res[[j]]        <-    res
  aicc_list[[j]]   <-   sum(sapply(res, `[[`, "AICc"))
  
}

names(Res)  =  names(aicc_list)  =  paste0("sim",1:length(file))
###############################################################
saveRDS(aicc_list, file=paste0(path,"deseq_aicc.rds"))
saveRDS(Res, file=paste0(path,"deseq_res.rds"))
