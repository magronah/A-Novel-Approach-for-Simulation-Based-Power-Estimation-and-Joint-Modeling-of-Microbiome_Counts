library(glmmTMB)
library(ggplot2)
library(dplyr)
library(AICcmodavg) 
##############################################################
path   =   paste0(getwd(),"/real_data/atlass_data/results/")
##############################################################
     rr       =   readRDS(paste0(path,"rr_mod.rds"))
     rrzi     =   readRDS(paste0(path,"rrzi_mod.rds"))
     us       =   readRDS(paste0(path,"us_mod.rds"))
     uszi     =   readRDS(paste0(path,"uszi_mod.rds"))
     nbmm     =   readRDS(paste0(path,"nbmm_aicc.rds"))
   zinbmm     =   readRDS(paste0(path,"zinbmm_aicc.rds"))
    #deseq     =   readRDS(paste0(path,"zinbmm_aicc.rds"))
##############################################################
    mod     =   lst(rr,rrzi,us,uszi)
      dd    =   data.frame(aictab(mod, modnames = names(mod)))
##############################################################
df  <- dd %>%
        select(Modnames,AICc, Delta_AICc)
      
a   =  min(df$AICc) -  zinbmm
if(z < 0){
  df 
}
  

      
      
            