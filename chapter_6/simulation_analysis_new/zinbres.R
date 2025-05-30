setwd("/home/agronahm/projects/def-bolker/agronahm/longitudinal_RR/")
library(NBZIMM)
library(foreach)
library(dplyr)
library(purrr)
library(tibble)
library(huge)
library(glmmTMB)
library(Matrix)
source("func.R")
############################################################
param_index   =   2
source("simulation_analysis_new/initial_param0.R")
############################################################
path = paste0("~/scratch/long/new_sim/",nsubj,"_",ntaxa,"_",ntime,"/zinb")
path
############################################################
files =   list.files(path, full.names = TRUE)

res   = foreach(i = files,.errorhandling = "remove",
                .packages = "NBZIMM") %do% {
  mod   =   readRDS(i)
  pp    =   fixed(mod)$dist
  dd    =  pp[(pp$variables) == "grouptreatment:time",][c("Estimate",
                                                          "Std.Error",
                                                          "padj")] %>%
            rownames_to_column("taxon")
  dd   =   dd %>%
    mutate(taxon = sub("--grouptreatment:time", "", taxon))
}
############################################################
res1 <- discard(res, ~ inherits(.x, "simpleError"))
common_taxa <- reduce(map(res1, ~ .x$taxon), intersect)
res_filtered <- map(res1, ~ filter(.x, taxon %in% common_taxa))
######################################################################
ll = length(res1)
padj <- (res_filtered
        %>% setNames(paste0("padjust", 1:ll))
        %>% purrr::map_dfr(pull, padj)
        %>% data.frame()
)
######################################################################
est <- (res_filtered
         %>% setNames(paste0("Estimate", 1:ll))
         %>% purrr::map_dfr(pull, Estimate)
         %>% data.frame()
)
#################################################################
sdr <- (res_filtered
        %>% setNames(paste0("sdr", 1:ll))
        %>% purrr::map_dfr(pull, Std.Error)
        %>% data.frame()
)
rownames(padj) = rownames(sdr) = rownames(est) = res_filtered[[1]]$taxon
#################################################################
df	=   lst(padj, est, sdr)
saveRDS(df, file = paste0("simulation_analysis_new/sim_data/",
                          nsubj,"_",ntaxa,"_",ntime,"/results/zinb.rds"))
