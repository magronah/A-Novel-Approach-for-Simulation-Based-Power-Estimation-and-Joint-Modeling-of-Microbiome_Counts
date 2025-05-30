setwd("/home/agronahm/projects/def-bolker/agronahm/longitudinal_RR/")
library(foreach)
library(huge)
library(glmmTMB)
library(Matrix)
source("func.R")
############################################################
param_index  =  2
#########################################################
source("simulation_analysis/initial_param0.R")
path1  =  "simulation_analysis/sim_data/"

path = paste0("~/scratch/long/",nsubj,"_",ntaxa,"_",ntime,"/rr")
path

files <-   list.files(path, full.names = TRUE)

error_files <- c()

res = foreach(i = files, .combine = "cbind") %do% {
  tryCatch({
    mod  = readRDS(i)
    dd   = ranef(mod, condVar = FALSE)$cond$taxon$`grouptreatment:time`
    dd
  }, error = function(e) {
    message("Error in file: ", i)
    error_files <<- c(error_files, i)
    NULL  # Return NULL so foreach can keep going
  })
}

#res = foreach(i = files, .combine ="cbind", .errorhandeling = "remove") %do% {
#    mod  =   readRDS(i)
#    dd   =   ranef(mod,condVar = FALSE)$cond$taxon$`grouptreatment:time` 
#    dd 
#}


dd            =  as.data.frame(res)
rownames(dd)  =  paste0("taxon",1:nrow(dd))
files_sub     =  setdiff(files,error_files)

sim_numbers <- gsub("\\D+", "", basename(files_sub))
colnames(dd)  =  paste0("sim", sim_numbers)
saveRDS(dd, file = paste0(path1,nsubj,"_",ntaxa,"_",ntime,"/results/rr.rds"))



