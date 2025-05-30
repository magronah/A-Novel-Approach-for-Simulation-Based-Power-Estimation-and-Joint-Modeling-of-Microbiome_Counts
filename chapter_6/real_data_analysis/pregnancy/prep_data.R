library(glmmTMB)
library(RhpcBLASctl)
library(NBZIMM)
library(tidyverse)
View(otu_data)
#########################################################
otu_data     =  Romero$OTU
meta_dd      =  Romero$SampleData
#########################################################
include_taxa    =   readRDS("real_data_analysis/pregnancy/results/include_taxa.rds")
          ddd    =   otu_data[,include_taxa]
#########################################################
meta_dd$normalizer   =   log(meta_dd$Total.Read.Counts)
     meta_df         =   meta_dd %>% 
                               select("pregnant", "GA_Days", 
                               "normalizer","Subect_ID","Sample_ID")  %>% 
                               data.frame() %>% 
                               na.omit()
   
dd    =    ddd[rownames(ddd) %in% rownames(meta_df),]
mean_count   =  colMeans(dd)
ntax         =  ncol(dd)          
################################################################
df = (dd
      |> rownames_to_column("Sample")
      |> pivot_longer(-Sample, names_to="taxon", values_to = "count") 
      |> as.data.frame()
)

met_dd = (meta_dd
          |> rownames_to_column("Sample") 
          |> as.data.frame()
)
################################################################
dd_long         =    left_join(df, met_dd, by ="Sample") 
dd_long         =    as.data.frame(dd_long)
dd_long$nugget  =    1:nrow(dd_long)

pvars <- c("pregnant", "GA_Days", "count",
         "normalizer", "taxon",
         "nugget", "Subect_ID")
ddf         =   dd_long %>% 
                 select(any_of(pvars))  %>% 
                 as.data.frame() %>%
                 na.omit()
##########################################################
ddf$taxon     =   as.factor(ddf$taxon)
ddf$pregnant  =   as.factor(ddf$pregnant)
# dim(dd)
############################################################
file_dir  = "real_data_analysis/pregnancy/results/"
saveRDS(meta_df, file =  paste0(file_dir, "metadata.rds"))
saveRDS(ddf, file =  paste0(file_dir, "dd_long.rds"))
saveRDS(mean_count, file =  paste0(file_dir, "mean_count.rds"))
saveRDS(dd, file =  paste0(file_dir, "countdata.rds"))
saveRDS(ntax, file =  paste0(file_dir, "ntaxa.rds"))
