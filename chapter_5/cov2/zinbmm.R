setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")

library(NBZIMM)
library(DESeq2)
library(Matrix)
library(huge)

source("func2.R")
source("initial_param0.R")
path = paste0("cov2/",nsubj,"_",ntaxa,"/")
path
####################################################################
data	=   readRDS(paste0(path,"sim_data/","zinbmm_otu_meta_list_withzi_taxa.rds"))
################################################################
cc	=   commandArgs(trailingOnly  = TRUE)
i	=   as.integer(cc[1])
dd	=   data[[i]]
################################################################
countdata  =   dd$countdata
met_dd     =   dd$met_data
met_dd$dummy = factor(1)

otu_count    =   t(countdata)
  dds        =   DESeqDataSetFromMatrix(otu_count,met_dd, ~group)
  dds        =   DESeq(dds,sfType ="poscounts",minReplicatesForReplace=Inf) 
  normalizer =   sizeFactors(dds) 


mod        =   mms(y = countdata, fixed = ~group + offset(normalizer),
                      random =  ~ 1|dummy,
                      zi_fixed = ~1,
                      data = met_dd, method = "zinb")

file_path  =  paste0("~/scratch/coverage","/",nsubj,"_",ntaxa,"/","zinbmm/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}



saveRDS(mod, file=paste0(file_path,"mod",i,".rds"))

