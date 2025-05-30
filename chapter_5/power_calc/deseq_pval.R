library(tibble)
library(DESeq2)
library(Matrix)
library(huge)
library(foreach)

source("/project/6006158/agronahm/reduced_rank_mm/func2.R")
source("/project/6006158/agronahm/reduced_rank_mm/initial_param0.R")
############################################################
# Define sets of nsubj and ntaxa
parameter_sets <- list(
  #list(nsubj = 100, ntaxa = 300),
  #list(nsubj = 150, ntaxa = 500),
  list(nsubj = 200, ntaxa = 600)
)

# Loop through parameter sets
for (params in parameter_sets) {
  
  nsubj <- params$nsubj
  ntaxa <- params$ntaxa
  
  path = paste0("~/scratch/dataset/RR","/",nsubj,"_",ntaxa,"/deseq")
  
  files   =   list.files(path, full.names = TRUE)
  res = foreach(i = files, .combine = "cbind") %do% {
    res        =    readRDS(i)
    res$padj
  }
  
  dd            =    as.data.frame(res)
  colnames(dd)  =    paste0("nsim", 1:ncol(dd))
  rownames(dd)  =    paste0("taxon",1:ntaxa)
  
  # Create output directory if it doesn't exist
  path2 <- paste0("/project/6006158/agronahm/reduced_rank_mm/", nsubj, "_", ntaxa, "/")

  file_path <- paste0(path2, "pvalues/")
  if (!dir.exists(file_path)) {
    dir.create(file_path, recursive = TRUE)
    cat("Folder created at:", file_path, "\n")
  } else {
    cat("Folder already exists at:", file_path, "\n")
  }
  
  saveRDS(dd, file = paste0(file_path,"deseq.rds"))
  
}
