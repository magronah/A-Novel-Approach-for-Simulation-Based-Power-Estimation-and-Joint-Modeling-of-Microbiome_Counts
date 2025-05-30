#setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
path        =   paste0(getwd(),"/reproducible/new_sim/real_data/")
##############################################################
library(NBZIMM)
source(paste0(path,"fun.R"))
source(paste0(path,"autism/load_data.R"))
##############################################################
met_dd$dummy  =  factor(1)
dd = t(dd_filt)
###############################################################
tt=system.time(
  mod    <-  mms(y =  dd, fixed = ~ group,
                 random = ~ 1|dummy,
                 zi_fixed = ~1, 
                 data = met_dd, method = "zinb")
)

(mod$fit[[1]]$zi.fit)

saveRDS(mod, file=paste0(path,"autism/data/zinbmm_mod.rds"))
saveRDS(tt, file=paste0(path,"autism/data/zinbmm_runtime",".rds"))
