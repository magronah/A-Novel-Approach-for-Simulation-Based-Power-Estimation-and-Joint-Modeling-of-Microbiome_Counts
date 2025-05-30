library(glmmTMB)
library(ggplot2)
library(dplyr)
library(AICcmodavg) 
library(here) 
##############################################################
path1   =   paste0("real_data/")
source(paste0(path1,"/fun.R"))

path   =   paste0("real_data/CrohnD_data/results/")
##############################################################
     rr       =   readRDS(paste0(path,"rr_mod.rds"))
     rrzi     =   readRDS(paste0(path,"rrzi_mod.rds"))
     us       =   readRDS(paste0(path,"us_mod.rds"))
     uszi     =   readRDS(paste0(path,"uszi_mod.rds"))
     nbmm     =   readRDS(paste0(path,"nbmm_aicc.rds"))
   zinbmm     =   readRDS(paste0(path,"zinbmm_aicc.rds"))
    deseq     =   readRDS(paste0(path,"deseq_aicc.rds"))
##############################################################
    mod     =   lst(rr,rrzi,us,uszi)
      dd    =   data.frame(aictab(mod, modnames = names(mod)))
##############################################################
df1  <- dd %>% select(Modnames,AICc)
df2  <- data.frame(Modnames =  c("nbmm", "zinbmm", "deseq"),
                   AICc     =  c(nbmm,zinbmm, deseq))

df   =   rbind(df1, df2)
df$Delta_AICc  <-  df$AICc - min(df$AICc)
df$Modnames   = c("RRzi","RR","USzi","US","Nbmm","Zinbmm","DE" )
##############################################################
ggplot(df, aes(Modnames,Delta_AICc)) + 
  geom_point() +
  custom_theme(12) +
  labs(x ="", y = "Change in AICc")




      
      
            