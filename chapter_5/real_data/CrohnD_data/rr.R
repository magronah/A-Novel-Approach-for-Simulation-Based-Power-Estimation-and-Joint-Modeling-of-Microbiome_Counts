library(RhpcBLASctl)
library(glmmTMB)
###########################################################
path   =   paste0(getwd(),"/real_data/CrohnD_data/")
source(paste0(path,"prep_data.R"))
################################################################
par_ctrl <- glmmTMBControl(
  parallel = list(n = 10, autopar = TRUE)
)

options(glmmTMB_openmp_debug = TRUE)
blas_set_num_threads(1)
################################################################
form =  count ~ 1 + age + us(1 + group|taxon) +
               rr(0 + taxon | subject,2) +  offset(normalizer)
################################################################
tt = system.time(
  fit  <-  glmmTMB(form, data = df,
                    family  = nbinom2, 
                    prior   = gprior,
                    REML    = FALSE,
                    control = par_ctrl
  )
)
################################################################
file_path  =  paste0(path,"results/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}
###########################################################
saveRDS(tt, file=paste0(file_path,"rr_runtime.rds"))
saveRDS(fit, file=paste0(file_path,"rr_mod.rds"))


