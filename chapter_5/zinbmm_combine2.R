setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")

library(foreach)
library(NBZIMM)

ntaxa   =  100

directory_path <-   "~/scratch/dataset/zinbmm_res/"
mod_path <-   "~/scratch/dataset/zinbmm_mod/"
files          <-   list.files(directory_path, full.names = TRUE)

names_list   =   list()
for(i in 1:100){
  mod              =    readRDS(paste0(mod_path,"mod",i,".rds"))
  names_list[[i]]  =    names(mod$fit)
}
names(names_list)  =  paste0("mod", 1:ntaxa)


common_names <- Reduce(intersect, names_list)
print(length(common_names))

res   = foreach(i = files, .combine="rbind") %do% {
    tryCatch({
    est   =  readRDS(i)
    est1  =  est[names(est) %in% common_names]
    stopifnot(names(est1) == common_names)
    est1
   }, error = function(e) {
      message(paste("Error in file:", i))
      return(NA)  # Return NA or another value that makes sense in your context
    })
  }

res            =    as.data.frame(res)
colnames(res)  =    common_names
rownames(res)  =    paste0("taxon",1:ntaxa)
saveRDS(res, file = paste0(getwd(), "/reproducible/simp_single/data/zinbmm_parametric.rds"))



