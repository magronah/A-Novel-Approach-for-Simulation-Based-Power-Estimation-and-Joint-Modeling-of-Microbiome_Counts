library(RhpcBLASctl)
library(glmmTMB)
###########################################################
path   =   paste0(getwd(),"/real_data/soil_data/")
source(paste0(path,"prep_data.R"))
################################################################
par_ctrl <- glmmTMBControl(
  parallel = list(n = 10, autopar = TRUE)
)

options(glmmTMB_openmp_debug = TRUE)
blas_set_num_threads(1)
################################################################
form =  count ~ 1 + us(1 + group|taxon) + rr(0 + taxon | site,2) +  offset(normalizer)
################################################################
tt1 = system.time(
  fit1 <- glmmTMB(formula = form, 
                  data = df, 
                  family = nbinom2, 
                  ziformula = ~1, 
                  prior = gprior, 
                  REML = FALSE, 
                  control = par_ctrl)
)

tt2 = system.time(
  fit2 <- glmmTMB(formula = form, 
                  data = df, 
                  family = nbinom2, 
                  ziformula = ~1 + (1 | taxon),
                  prior = gprior, 
                  REML = FALSE, 
                  control = par_ctrl)
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
saveRDS(tt1, file=paste0(file_path,"rrzi_runtime.rds"))
saveRDS(fit1, file=paste0(file_path,"rrzi_mod.rds"))

saveRDS(tt2, file=paste0(file_path,"rrzi_each_runtime.rds"))
saveRDS(fit2, file=paste0(file_path,"rrzi_each_mod.rds"))
