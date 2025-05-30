library(phyloseq)
library(dplyr)
library(fido)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(microbiome)
library(vegan)
library(pheatmap)
library(here)
##############################################################
path1   =   paste0("real_data/")
source(paste0(path1,"fun.R"))
##############################################################
data("atlas1006")

countdata  <-  data.frame(abundances(atlas1006)) 
colnames(countdata) <-  paste0("Sample-",1:ncol(countdata))
rownames(countdata) <-  paste0("taxon",1:nrow(countdata))

#' The full dataset has 130 taxa and 1151 subjects
#' We analyse only the data at the initial time point. 
#' The data at the initial time point has
#' uneven time intervals

meta_data  <-  sample_data(atlas1006)
meta_dd    <-  meta_data %>%
               data.frame()  %>%
               dplyr::select(c("subject","bmi_group","age")) %>% 
               dplyr::arrange(subject) %>%
               dplyr::distinct(subject, .keep_all = TRUE) %>% 
               setNames(c("subject","group", "age")) %>%
               na.omit() 

meta_dd$group <- factor(meta_dd$group, 
                        levels = c("lean", "underweight", "overweight", 
                                   "obese", "severeobese", "morbidobese"))
meta_dd$subject  =  rownames(meta_dd)
data  <-   countdata %>%  
              dplyr::select(rownames(meta_dd))
data  <-   data.frame(t(data))
###############################################################################
dd1  <- data[, apply(data, 2, function(col) any(col == 0))]
dd2  <- dd1[, apply(dd1, 2, function(col) sum(col) != 0)] 

#rownames(meta_dd)  = paste0("subject", 1:nrow(meta_dd))
#rownames(dd2)  = paste0("subject", 1:nrow(dd2))

ntax    =  ncol(dd2)
##############################################################
dd_long    =   df_long(dd2, otu_names = "taxon", 
                       subject_name = "subject", 
                       ntaxa=ntax)

dd_long    =   left_join(dd_long, meta_dd, by ="subject")  
###########################################################
dd         =    t(dd2)
tt = system.time(
pp         <-    deseqfun(dd, meta_dd, alpha_level=0.1,
                         design   = ~group + age,
                         ref_name = "lean",
                         ntaxa = ntax)
)
###########################################################
countdata     =    pp$data$countdata
meta_dd       =    pp$data$meta_data
normalize_fac =    sizeFactors(pp$object)
meta_dd$normalizer  =   normalize_fac
###########################################################
normalizer =   data.frame(normalizer = normalize_fac) %>% 
  rownames_to_column("subject")

df    =   left_join(dd_long, normalizer, by ="subject")  
################################################################
gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")


mean_count   =  rowMeans(countdata)
################################################################
file_dir  = "atlass_data/results/"

saveRDS(tt, file =  paste0(path1, file_dir, "deseq_runtime.rds"))

saveRDS(pp$result, file =  paste0(path1, file_dir, "deseq_mod_est.rds"))
saveRDS(meta_dd, file =  paste0(path1, file_dir, "metadata.rds"))
saveRDS(mean_count, file =  paste0(path1, file_dir, "mean_count.rds"))
saveRDS(countdata, file =  paste0(path1, file_dir, "countdata.rds"))
saveRDS(ntax, file =  paste0(path1, file_dir, "ntaxa.rds"))

