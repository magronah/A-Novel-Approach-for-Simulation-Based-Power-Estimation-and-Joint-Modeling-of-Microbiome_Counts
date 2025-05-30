library(tidyr)
library(tidyverse)
library(tibble)
library(here)
library(DESeq2)
library(RhpcBLASctl)
library(glmmTMB)
library(ggplot2)
library(AICcmodavg)
##############################################################
path   =   paste0(getwd(),"/real_data/")
source(paste0(path,"/fun.R"))
data        =   readRDS(paste0(path,"aut_data.rds"))
metadata    =   readRDS(paste0(path,"aut_metadata.rds"))

nam         =   "PRJNA644763"
##############################################################
rr1     =  readRDS(paste0(path,"autism/data/rr_mod",nam,".rds"))
rrzi   =  readRDS(paste0(path,"autism/data/rrzi_mod",nam,".rds"))
us     =  readRDS(paste0(path,"autism/data/us_mod",nam,".rds"))
uszi   =  readRDS(paste0(path,"autism/data/uszi_mod",nam,".rds"))
#uszi_each   =  readRDS(paste0(path,"autism/data/uszi_each_mod",nam,".rds"))
#rr3     =  readRDS(paste0(path,"autism/data/rr_mod_rr=3",nam,".rds"))

nbmm   =  readRDS(paste0(path,"autism/data/nbmm_mod",nam,".rds"))
zinbmm   =  readRDS(paste0(path,"autism/data/zinbmm_mod",nam,".rds"))

##############################################################
models  =  lst(rr,rrzi,us,uszi)
aic_dd  =  aictab(models, modnames = names(models))
aic_dd
##############################################################
bbmle::AICtab(rr,rrzi,us,uszi, logLik =TRUE, base =TRUE)

aics1  =   aictab(nbmm)
aics2  =   aictab(zinbmm)

aics_nbmm <- sum(aics1[!is.na(aics1$AICc), ]$AICc)
aics_zinbmm <- sum(aics2[!is.na(aics2$AICc), ]$AICc)


aic_dd$nbmm  =  sum(aics_clean$AICc)
aic_dd$nbmm  =  sum(aics_zinbmm$AICc)

aic0   =  data.frame(AIC=c(AIC(rrzi),AIC(uszi),AIC(rr),AIC(us)))
rownames(aic0)= c("rrzi","uszi","rr","us")

library(MuMIn)
aic   =  data.frame(AICc=c(AICc(rrzi),AICc(uszi),AICc(rr),AICc(us)))
rownames(aic)= c("rrzi","uszi","rr","us")


library(ggplot2)
library(dplyr)

# Data
model_results <- data.frame(
  Model = c("mod1", "mod3", "mod2"),
  dAIC = c(0.0, 0.1, 0.5),
  df = c(4, 5, 5)
)

# Create the plot
ggplot(model_results, aes(x = reorder(Model, dAIC), y = dAIC, fill = Model)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = round(dAIC, 2), y = dAIC + 0.05), size = 5, vjust = 0) +
  labs(
    title = "Model Comparison based on dAIC",
    x = "Model",
    y = "Delta AIC (dAIC)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


#the zero inflation term is important 
#' us is better than rr
#' The reduced rank term is goo

rr    =  readRDS(paste0(path,"aut_data.rds"))

