library(tibble)
library(dplyr)
library(Matrix)
library(glmmTMB)
library(AICcmodavg)
library(huge)
library(here)
####################################################################
source("func2.R")
source("initial_param0.R")
####################################################################
path = paste0(nsubj,"_",ntaxa,"/")
path
######################################################################
file_path  =   paste0("~/scratch/data/",nsubj,"_",ntaxa,"/rr/")
file       =   list.files(file_path, full.names = TRUE)
######################################################################
mod_list   =   list()

for(i in 1:length(file)){
  mod_list[[i]]     =    readRDS(file[i])
}  
names(mod_list)     =    paste0("sim",1:length(file))
dd    =   data.frame(aictab(mod_list, modnames = names(mod_list)))
ddd   =   dd$AICc 
names(ddd)   =   dd$Modnames

saveRDS(ddd, file=paste0(path,"rr_aicc.rds"))

######################################################################
file_path  =   paste0("~/scratch/data/",nsubj,"_",ntaxa,"/rrzi/")
file       =   list.files(file_path, full.names = TRUE)
######################################################################
mod_list   =   list()

for(i in 1:length(file)){
  mod_list[[i]]     =    readRDS(file[i])
}  
names(mod_list)     =    paste0("sim",1:length(file))
dd    =   data.frame(aictab(mod_list, modnames = names(mod_list)))
ddd   =   dd$AICc 
names(ddd)   =   dd$Modnames

saveRDS(ddd, file=paste0(path,"rrzi_aicc.rds"))
######################################################################
file_path  =   paste0("~/scratch/data/",nsubj,"_",ntaxa,"/uszi/")
file       =   list.files(file_path, full.names = TRUE)
######################################################################
mod_list   =   list()

for(i in 1:length(file)){
  mod_list[[i]]     =    readRDS(file[i])
}  
names(mod_list)     =    paste0("sim",1:length(file))
dd    =   data.frame(aictab(mod_list, modnames = names(mod_list)))
ddd   =   dd$AICc 
names(ddd)   =   dd$Modnames

saveRDS(ddd, file=paste0(path,"uszi_aicc.rds"))
######################################################################

file_path  =   paste0("~/scratch/data/",nsubj,"_",ntaxa,"/us/")
file       =   list.files(file_path, full.names = TRUE)
######################################################################
mod_list   =   list()

for(i in 1:length(file)){
  mod_list[[i]]     =    readRDS(file[i])
}  
names(mod_list)     =    paste0("sim",1:length(file))
dd    =   data.frame(aictab(mod_list, modnames = names(mod_list)))
ddd   =   dd$AICc 
names(ddd)   =   dd$Modnames

saveRDS(ddd, file=paste0(path,"us_aicc.rds"))
######################################################################