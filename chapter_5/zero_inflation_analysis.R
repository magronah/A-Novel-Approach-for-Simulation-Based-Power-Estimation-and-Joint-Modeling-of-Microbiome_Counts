setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
source("reproducible/new_sim/func.R")
source("reproducible/new_sim/param.R")
source("reproducible/new_sim/Load_Packages.R")
# source("reproducible/new_sim/initial_param0.R")
path = paste0(getwd(),"/reproducible/new_sim/data/")
####################################################################
dd          =   readRDS(paste0(path,"otu_meta_list_nozi.rds"))
zi_dd1       =   readRDS(paste0(path,"otu_meta_list_withzi.rds"))
zi_taxa_dd	=   readRDS(paste0(path,"otu_meta_list_withzi_taxa.rds"))

zi_dd       =   readRDS(paste0(path,"otu_meta_list_withzi_test.rds"))
zi_taxa_dd	=   readRDS(paste0(path,"otu_meta_list_withzi_taxa_test.rds"))

i   =   1
View((zi_dd[[i]]$countdata))
View(zres_dd2[[i]]$countdata)

plot(rowSums(zi_dd[[i]]$countdata), rowSums(zi_dd1[[i]]$countdata))

plot(rowSums(ztaxa_res_dd2[[i]]$countdata), rowSums(zres_dd2[[i]]$countdata))
####################################################################
zi_dd_count      <- lapply(zi_dd, function(x) x$countdata)
zi_taxa_dd_count <- lapply(zi_taxa_dd, function(x) x$countdata)

zi_dd_count_merge       <- do.call(rbind, zi_dd_count)
zi_taxa_dd_count_merge  <- do.call(rbind, zi_taxa_dd_count)
dim(zi_taxa_dd_count_merge)

pp1   = apply(zi_dd_count_merge,2, function(x) sum(x == 0) / length(x))
pp2   = apply(zi_taxa_dd_count_merge,2, function(x) sum(x == 0) / length(x))

ddf   =  data.frame(prop = c(pp1, pp2), 
                    type = rep(c("pp1","pp2"), each = length(pp1)))

ggplot(ddf, aes(x = prop, group=type, color = type)) +
  geom_density() 
