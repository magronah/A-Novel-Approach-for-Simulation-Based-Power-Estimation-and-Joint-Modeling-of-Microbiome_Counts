library(tibble)
library(DESeq2)
library(Matrix)
library(huge)
library(here)

source("func2.R")
source("initial_param0.R")
####################################################################
path = paste0(nsubj,"_",ntaxa,"/")
path
####################################################################
data      =   readRDS(paste0(path,"otu_meta_list_withzi_taxa.rds"))
################################################################
cc      =   commandArgs(trailingOnly  = TRUE)
i       =   as.integer(cc[1]) 

for(i  in  1:5){
   dd      =   data[[i]]
################################################################
countdata  =   dd$countdata
met_dd     =   dd$met_data

res        =   deseqfun(countdata,met_dd,ref_name="control",ntaxa=ntaxa)
file_path  =   paste0("~/scratch/data/",nsubj,"_",ntaxa,"/","deseq/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}

saveRDS(res, file=paste0(file_path,"mod",i,".rds"))
}
