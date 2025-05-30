library(foreach)
library(huge)
library(Matrix)
source("func2.R")
source("initial_param0.R")
############################################################
path = paste0("~/scratch/dataset/RR/coverage/",nsubj,"_",ntaxa,"/deseq")
path

files          <-   list.files(path, full.names = TRUE)

res = foreach(i = files, .combine = "cbind") %do% {
     r        =    readRDS(i)
     r$padj
}

dd            =    as.data.frame(res)
colnames(dd)  =    paste0("nsim", 1:ncol(dd))
rownames(dd)  =    paste0("taxon",1:nrow(dd))

saveRDS(dd, file = paste0("cov2/",nsubj,"_",ntaxa,"/deseq_pval.rds"))



