library(NBZIMM)
library(DESeq2)
library(Matrix)
library(huge)
library(here)

source("func2.R")
source("initial_param0.R")

path = paste0(nsubj,"_",ntaxa,"/")
path
####################################################################
data	  =   readRDS(paste0(path,"otu_meta_list_withzi_taxa.rds"))
################################################################
cc	=   commandArgs(trailingOnly  = TRUE)
i	  =   as.integer(cc[1])
for(i in 1:5){
dd	=   data[[i]]
################################################################
countdata  =   dd$countdata
met_dd     =   dd$met_data
met_dd$dummy = factor(1)

otu_count   =   t(countdata)
dds        =   DESeqDataSetFromMatrix(otu_count,met_dd, ~group)
dds        =   DESeq(dds,sfType ="poscounts",minReplicatesForReplace=Inf)
normalizer =   sizeFactors(dds)

mod        =   mms(y = countdata, fixed = ~group + offset(normalizer),
                   random =  ~ 1|dummy,
                   data = met_dd, method = "nb")

file_path  =  paste0("~/scratch/data/",nsubj,"_",ntaxa,"/","nbmm/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}



saveRDS(mod, file=paste0(file_path,"mod",i,".rds"))
}

