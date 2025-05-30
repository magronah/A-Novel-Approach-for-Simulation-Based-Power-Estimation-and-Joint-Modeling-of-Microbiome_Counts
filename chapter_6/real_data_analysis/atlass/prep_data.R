library(dplyr)
library(tidyverse)
library(microbiome)
library(NBZIMM)
data("atlas1006")
path = "real_data_analysis/atlass/"
##############################################################
countdata  <-  data.frame(abundances(atlas1006)) 
rownames(countdata) <-  paste0("taxon",1:nrow(countdata))

include_taxa    <-    readRDS(paste0(path,"include_taxa.rds"))
meta_data       <-    sample_data(atlas1006)
##############################################################
meta_dd    <-  meta_data %>%
  data.frame()  %>%
  dplyr::select("subject","bmi_group","age","sample","time") %>% 
  na.omit()

#removed missing data
meta_dd$sample <- gsub("-", ".", meta_dd$sample)
rownames(meta_dd) <- gsub("-", ".", rownames(meta_dd))

count_dd  <-  countdata %>%
  dplyr::select(meta_dd$sample)  

ddd  =  as.data.frame(t(count_dd))
ddd  =  ddd %>% 
   dplyr::select(all_of(include_taxa))

meta_dd$normalizer   =   log(colSums(count_dd))
##############################################################
meta_dd$bmi_group <- factor(meta_dd$bmi_group, 
                            levels = c("lean", "underweight", "overweight", 
                                       "obese", "severeobese", "morbidobese"))
mean_count   =  colMeans(ddd)
ntax     =  ncol(ddd)
################################################################
df = (ddd
      |> rownames_to_column("sample")
      |> pivot_longer(-sample, names_to="taxon", values_to = "count") 
      |> as.data.frame()
)
################################################################
ddf         =    left_join(df, meta_dd, by ="sample") 
ddf$nugget  =    1:nrow(ddf)
##########################################################
ddf$taxon     =   as.factor(ddf$taxon)
##########################################################
saveRDS(meta_dd, file =  paste0(path, "results/metadata.rds"))
saveRDS(ddf, file =  paste0(path, "results/dd_long.rds"))
saveRDS(mean_count, file =  paste0(path, "results/mean_count.rds"))
saveRDS(ddd, file =  paste0(path, "results/countdata.rds"))
saveRDS(ntax, file =  paste0(path, "results/ntaxa.rds"))

