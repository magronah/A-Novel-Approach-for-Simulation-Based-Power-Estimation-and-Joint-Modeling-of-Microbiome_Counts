setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")
library(MASS)
library(nlme)
library(NBZIMM)
library(foreach)
library(huge)
library(glmmTMB)
library(Matrix)
source("func2.R")
source("initial_param0.R")
############################################################
path = paste0("~/scratch/coverage/",nsubj,"_",ntaxa,"/zinbmm")
path


files =  list.files(path, full.names = TRUE)

sub    =  500
list1  =  list2  =  list()

for(i in files[1:sub]){
     mod   =   readRDS(i)
     pp    =   fixed(mod)$dist
     ppp   =   pp[(pp$variables) == "grouptreat",][["Estimate"]]
     names(ppp)  =  mod$response
     list1[[i]]  =  ppp
     list2[[i]]  =  unlist(lapply(mod$fit, function(x){
                                        (summary(x))$tTable["grouptreat", "p-value"]}))

# Print the result
# unlist(lapply(mod,function(x) {summary(x)$coefficients$cond[, "Pr(>|z|)"][["grouptreat"]]}))
}

common_names1 <- Reduce(intersect, lapply(list1, names))
filtered_res1  <- lapply(list1, function(x) x[common_names1])
res1	 <- do.call(cbind, filtered_res1)

dd1    =   as.data.frame(res1)
print(common_names1)
#rownames(dd)  =  paste0("taxon",1:ntaxa)
colnames(dd1)  =    paste0("nsim",1:ncol(dd1))
saveRDS(dd1, file = paste0("cov2/",nsubj,"_",ntaxa,"/zinbmm.rds"))

###########################################################
common_names2 <- Reduce(intersect, lapply(list2, names))
filtered_res2  <- lapply(list2, function(x) x[common_names2])
res2     <- do.call(cbind, filtered_res2)

print(common_names2)
dd2    =   as.data.frame(res2)
colnames(dd2)  =  paste0("sim",1:ncol(dd2))
saveRDS(dd2, file = paste0("cov2/",nsubj,"_",ntaxa,"/pvalues/zinbmm.rds"))

stopifnot(common_names1 == common_names2)





if(FALSE){
setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")

library(MASS)
library(nlme)
library(NBZIMM)
library(foreach)
library(huge)
library(glmmTMB)
library(Matrix)
source("func2.R")
source("initial_param0.R")
############################################################
path = paste0("~/scratch/dataset/RR/coverage/",nsubj,"_",ntaxa,"/zinbmm")
path


files =   list.files(path, full.names = TRUE)

res   = foreach(i = files,.combine = "cbind",.packages = "NBZIMM") %do% {
     mod   =   readRDS(i)
     pp    =   fixed(mod)$dist
    ppp    =   pp[(pp$variables) == "grouptreat",][["Estimate"]]
    names(ppp)  =  mod$response
     ppp
}

common_names <- Reduce(intersect, lapply(res, names))
filtered_res <- lapply(res, function(x) x[common_names])
res	 <- do.call(cbind, filtered_res)


dd    =   as.data.frame(res)
#rownames(dd)  =	 paste0("taxon",1:ntaxa)
colnames(dd)  =    paste0("nsim",1:ncol(dd))
saveRDS(dd, file = paste0("cov2/",nsubj,"_",ntaxa,"/zinbmm.rds"))
}
