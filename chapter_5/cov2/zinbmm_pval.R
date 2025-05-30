setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")

library(foreach)
library(huge)
library(glmmTMB)
library(Matrix)
source("func2.R")
source("initial_param0.R")
############################################################
path = paste0("~/scratch/dataset/RR/coverage/",200,"_",600,"/zinbmm2")
path

files <-   list.files(path, full.names = TRUE)

res = foreach(i = files,.combine ="c", .packages = "glmmTMB") %do% {
  mod  =   readRDS(i)
  unlist(lapply(mod,function(x) {summary(x)$coefficients$cond[, "Pr(>|z|)"][["grouptreat"]]}))
}


common_names <- Reduce(intersect, lapply(res, names))
print(common_names)
filtered_res <- lapply(res, function(x) x[common_names])
res      <- do.call(cbind, filtered_res)

dd    =   as.data.frame(res)
#rownames(dd)  =  paste0("taxon",1:nrow(dd))
colnames(dd)  =  paste0("sim",1:ncol(dd))
saveRDS(dd, file = paste0("cov2/",nsubj,"_",ntaxa,"/zinbmm_pval.rds"))
