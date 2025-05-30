library(RhpcBLASctl)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(gllvm)
###########################################################
path1   =   paste0(getwd(),"/real_data/")
source(paste0(path1,"/fun.R"))
###########################################################
data("microbialdata")
data <- (microbialdata$Y); met_dd <- microbialdata$Xenv
dd  <- data[, apply(data, 2, function(col) any(col == 0))]
###########################################################
meta_dd    =   met_dd %>% 
  select(Site, Soiltype) %>% 
  rownames_to_column() %>% 
  setNames(c("sample","site","group")) 
#' group variable is the soiltype
#' subjects are the soil samples
##############################################################
ntax = ncol(dd)
dd_long    =   df_long(dd, otu_names = "OTU", 
                       subject_name = "sample", 
                       ntaxa = ntax)

dd_long    =   left_join(dd_long, meta_dd, by ="sample")  
###########################################################
tt = system.time(
  pp         <-    deseqfun(dd, meta_dd, alpha_level=0.1,
                            design   = ~group + site,
                            ref_name = "B",
                            ntaxa = ntax)
)


countdata     =    pp$data$countdata
meta_dd       =    pp$data$meta_data
normalize_fac =    sizeFactors(pp$object)
meta_dd$normalizer  =   normalize_fac
###########################################################
normalizer =   data.frame(normalizer = normalize_fac) %>% 
  rownames_to_column("sample")

df    =   left_join(dd_long, normalizer, by ="sample")  
################################################################
gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")
 
mean_count   =  rowMeans(countdata)
################################################################
file_dir  = "soil_data/results/"
saveRDS(tt, file =  paste0(path1, file_dir, "deseq_runtime.rds"))

saveRDS(pp$result, file =  paste0(path1, file_dir, "deseq_mod_est.rds"))
saveRDS(meta_dd, file =  paste0(path1, file_dir, "metadata.rds"))
saveRDS(mean_count, file =  paste0(path1, file_dir, "mean_count.rds"))
saveRDS(countdata, file =  paste0(path1, file_dir, "countdata.rds"))
saveRDS(ntax, file =  paste0(path1, file_dir, "ntaxa.rds"))

 