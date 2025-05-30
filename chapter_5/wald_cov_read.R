setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
source("reproducible/simp_single/Load_Packages.R")
source("reproducible/simp_single/func.R")
source("reproducible/simp_single/param.R")
source("reproducible/simp_single/initial_param0.R")
source("load_glmmTMB.R")

directory_path <-   "~/scratch/dataset/glmmTMB_mod/"
#################################################################
cc        =   commandArgs(trailingOnly  = TRUE)
i         =   as.integer(cc[1])

mod  =   readRDS(paste0(directory_path,"mod",i,".rds"))
  dd   =   wald_confint(mod, ntaxa, conf_level = .95, in_sd = 1)
dd$sim =   rep(i, nrow(dd))
dd	=  as.data.frame(dd)
saveRDS(dd, file = paste0("~/scratch/dataset/glmmTMB_wald/mod",i, ".rds"))

##########################################################################
if(FALSE){
files          <-   list.files(directory_path, full.names = TRUE)

ncores = Sys.getenv("SLURM_CPUS_PER_TASK") 
registerDoParallel(cores=ncores)


res = foreach(i = files, .combine ="rbind") %dopar% {
  mod  =   readRDS(i)
  dd   =   wald_confint(mod, ntaxa, conf_level = .95)
dd$sim =   rep(j, nrow(dd))
  dd 
}

dd      =  as.data.frame(res)
saveRDS(dd, file = paste0(getwd(), "/reproducible/simp_single/data/wald_cov.rds"))

dd$sim  =  rep(1:length(files),each=ntaxa)

saveRDS(dd, file = paste0(getwd(), "/reproducible/simp_single/data/wald_cov.rds"))
}
