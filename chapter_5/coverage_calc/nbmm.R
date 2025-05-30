setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
source("reproducible/new_sim/func.R")
source("reproducible/new_sim/param.R")
source("reproducible/new_sim/Load_Packages.R")
source("reproducible/new_sim/initial_param0.R")
path = paste0(getwd(),"/reproducible/new_sim/coverage_calc/data/")
####################################################################
data      =   readRDS(paste0(path,"otu_meta_list_nozi2.rds"))
################################################################
cc      =   commandArgs(trailingOnly  = TRUE)
n       =   as.integer(cc[1]) 
################################################################
step_size =   20
  start   =   (n - 1) * step_size + 1
  end     =   (n * step_size)
  dd	  =   data[start:end]

nam       =   start:end

for(i in 1:length(dd)){
 countdata  =   dd[[i]]$countdata
 met_dd     =   dd[[i]]$met_data
 met_dd$dummy = factor(1)
 mod        =   mms(y = countdata, fixed = ~group,
                   random = ~ 1|dummy,
                   data = met_dd, method = "nb")

saveRDS(mod, file=paste0("~/scratch/dataset/new_sim/coverage/nbmm_mod/mod",nam[i], ".rds"))   


nsim  =   length(mod$fit)
print(length(mod$fit))
print(names(mod$fit))

res   = foreach(j = 1:nsim, .combine="c", .packages = "NBZIMM") %do% {
    tryCatch({
         coef(mod$fit[[j]])[["grouptreat"]]
    }, error = function(e) {
      message(paste("Error in file:", i, "at simulation:", j, "-", e))
      return(NA)
    })
  }

names(res)  =  names(mod$fit)
saveRDS(res, file = paste0("~/scratch/dataset/new_sim/coverage/nbmm_res/res",nam[i],".rds" ))

}
