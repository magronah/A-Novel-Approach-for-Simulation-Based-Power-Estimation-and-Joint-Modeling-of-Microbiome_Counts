setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")

library(foreach)
library(doParallel)
library(NBZIMM)
#########################################################
## fixed doesnt run on the save model untill it is run once
data(Romero)
names(Romero)

otu = Romero$OTU; dim(otu)
sam = Romero$SampleData; dim(sam)
colnames(sam)

N = sam[, "Total.Read.Counts"]        
Days = sam$GA_Days; Days = scale(Days)
Age = sam$Age; Age = scale(Age)
Race = sam$Race
preg = sam$pregnant; table(preg)
subject = sam[, "Subect_ID"]; table(subject)

y = otu[, 1]

f1 = glmm.nb(y ~ Days + Age + Race + preg, random = ~ 1 | subject)
summary(f1)
data.frame(fixed(f1))["Days","dist.pvalue"]
#######################################################
directory_path <-   "~/scratch/dataset/new_sim/nbmm_mod/"

cc    =   commandArgs(trailingOnly  = TRUE)
i     =   as.integer(cc[1])

mod   =   readRDS(paste0(directory_path,"mod",i,".rds"))
df    =   fixed(mod)[[1]]
grouptreat <- grep("grouptreat", rownames(df), value = TRUE)

res   =  df[grouptreat, ]

# The first model:
f2 = mms(y = Romero$OTU, fixed = ~  Days + Age + Race + preg + offset(log(N)), 
         random = ~ 1 | subject, data = Romero$SampleData,
         min.p = 0.2, method = "nb")

fixed(f2)
#nsim  =   length(mod$fit)


#res   = foreach(j = 1:nsim, .combine="c", .packages = "NBZIMM") %do% {
#    tryCatch({
#       dd   =   fixed(mod)[[1]]
#       dd
#    }, error = function(e) {
#      message(paste("Error in file:", i, "at simulation:", j, "-", e))
#      return(NA)
#    })
#  }
ff = readRDS("~/Downloads/f2.rds")
fixed(ff)

saveRDS(f2,file= "~/Downloads/f2.rds")
#names(res)  =  names(mod$fit)
saveRDS(res, file = paste0("~/scratch/dataset/new_sim/nbmm_pval/pval",i,".rds" ))
