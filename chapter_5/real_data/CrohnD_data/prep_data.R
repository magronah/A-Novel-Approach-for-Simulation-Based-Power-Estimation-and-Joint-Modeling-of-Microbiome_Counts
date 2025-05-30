library(phyloseq)
library(dplyr)
library(fido)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(here)
##############################################################
path1   =   paste0("real_data/")
source(paste0(path1,"/fun.R"))
##############################################################
set.seed(899)
data(RISK_CCFA)

# making into a phyloseq object
CCFA_phylo <- phyloseq(otu_table(as.matrix(RISK_CCFA_otu), taxa_are_rows = TRUE),
                       sample_data(RISK_CCFA_sam), 
                       tax_table(as.matrix(RISK_CCFA_tax)))

dim(otu_table(CCFA_phylo,taxa_are_rows = TRUE)) # 9511 taxa 1359 samples
vv=(sample_data(RISK_CCFA_sam)) # 1359   73
#View(tax_table(CCFA_phylo))
unique(vv$diagnosis)
unique(vv$diagnosis)
# drop low abundant taxa and samples
dim(otu_table(dat,taxa_are_rows = TRUE))
dat <- CCFA_phylo %>% 
  subset_samples(disease_stat!="missing", 
                 immunosup!="missing") %>% 
  subset_samples(diagnosis %in% c("no", "CD")) %>% 
  subset_samples(steroids=="false") %>% 
  subset_samples(antibiotics=="false") %>% 
  subset_samples(biologics=="false")  %>% 
  subset_samples(biopsy_location=="Terminal ileum") %>% 
  tax_glom("Family") %>% 
  prune_samples(sample_sums(.) >= 5000,.) %>%
  filter_taxa(function(x) sum(x > 3) > 0.10*length(x), TRUE)

dim(otu_table(dat,taxa_are_rows = TRUE))

##' we consider various 
########################################################################
sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>% 
  mutate(age = as.numeric(as.character(age)),
         diagnosis = relevel(factor(diagnosis, ordered = FALSE), ref="no"), 
         disease_stat = relevel(factor(disease_stat, ordered = FALSE), ref="non-inflamed"))
########################################################
data    =   t(otu_table(dat))  %>% 
          as.data.frame()  %>% 
          setNames(paste0("taxon",1:ncol(t(otu_table(dat))))) 

dim(data)
dd  <-  data[, apply(data, 2, function(col) any(col == 0))]
dd  <-  dd[,-3]  #zinbmm cannot fit this as it says there is an NA
                 #I have not find NA though 

meta_dd  = sample_dat %>%  
          dplyr::select(diagnosis,age) %>%
          rownames_to_column()   %>%
          setNames(c("subject","group", "age")) 

ntax    =  ncol(dd)
##############################################################
dd_long    =   df_long(dd, otu_names = "taxon", 
                       subject_name = "subject", 
                       ntaxa=ntax)

dd_long    =   left_join(dd_long, meta_dd, by ="subject")  
###########################################################
tt  = system.time(
  pp        <-    deseqfun(dd, meta_dd, alpha_level=0.1,
                           design   = ~group + age,
                           ref_name = "CD",
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
file_dir  = "CrohnD_data/results/"
saveRDS(tt, file =  paste0(path1, file_dir, "deseq_runtime.rds"))

saveRDS(pp$result, file =  paste0(path1, file_dir, "deseq_mod_est.rds"))
saveRDS(meta_dd, file =  paste0(path1, file_dir, "metadata.rds"))
saveRDS(mean_count, file =  paste0(path1, file_dir, "mean_count.rds"))
saveRDS(countdata, file =  paste0(path1, file_dir, "countdata.rds"))
saveRDS(ntax, file =  paste0(path1, file_dir, "ntaxa.rds"))



