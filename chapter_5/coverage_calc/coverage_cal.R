setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
source("reproducible/new_sim/func.R")
source("reproducible/new_sim/Load_Packages.R")
###############################################################
path = paste0(getwd(),"/reproducible/new_sim/coverage_calc/data/")
true_param  =   readRDS(paste0(path, "true_param.rds"))
mean_count  =   readRDS(paste0(path,"mean_counts.rds"))

rr_est      =   readRDS(paste0(path,"rr_cov_parametric.rds"))
rr_est      =   lapply(rr_est, function(x){as.data.frame(x)} )
################################################################
dd_list     =   lapply(rr_est, function(x){dd_long(x,true_param,label="with_rr")})
para_conf   =   lapply(dd_list, function(x){para_confint(x,true_param, alpha = 0.05)})
cov_dd      =   as.data.frame(lapply(para_conf, function(x){coverage_cal(x)$coverage}))
true_param$coverage =  rowSums(cov_dd)

ggplot(true_param, aes(true_param, coverage)) +
  geom_point() +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red") + 
  geom_hline(yintercept = 90, linetype = "dashed", color = "red") + 
  theme_bw()


true_param$mean_val  = rowMeans(mean_count)

ggplot(true_param, aes(mean_val, coverage)) +
  geom_point() +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red") + 
  geom_hline(yintercept = 90, linetype = "dashed", color = "red") + 
  theme_bw()


175/2
View(true_param)



nbmm0      =   readRDS(paste0(path0,"nbmm_parametric.rds"))

rr1        =   readRDS(paste0(path1,"rr_cov_parametric.rds"))
deseq1     =   readRDS(paste0(path1,"deseq_cov_parametric.rds"))
nbmm1      =   readRDS(paste0(path1,"nbmm_cov_parametric.rds"))
################################################################
rr       =   cbind(rr0, rr1)
deseq    =   cbind(deseq0, deseq1)
nbmm     =   cbind(nbmm0, nbmm1)
colnames(rr)  =  colnames(deseq)  =   colnames(nbmm)  =  paste0("sim",1:ncol(rr))
################################################################
rr_list   =   deseq_list  =  nbmm_list  =  list()  
num_columns =   ncol(rr)
slice_size  =   200

for (i in seq(1, num_columns, by = slice_size)) {
  end <- min(i + slice_size - 1, num_columns)  
  rr_list[[length(rr_list) + 1]] <- rr[, i:end]
  deseq_list[[length(deseq_list) + 1]] <- deseq[, i:end]
  nbmm_list[[length(nbmm_list) + 1]] <- nbmm[, i:end]
}

res_rr  =  foreach(i = 1:length(rr_list), .combine="rbind") %do% {
  rr_dd =  dd_long(rr_list[[i]],true_param,label="with_rr")
  rr_para_conf      =   para_confint(rr_dd, true_param, alpha = 0.05)
  coverage_cal(rr_para_conf)
}

res1   =  res_rr %>%
            group_by(param_name) %>%
           summarise(sum(coverage)/100)


View(res1)



res[res$param_name %in% "taxon2", ]
dim(res)
deseq_dd[[i]]     =   dd_long(deseq_list[[i]],true_param,label="deseq")
nbmm_dd[[i]]      =   dd_long(nbmm_list[[i]],true_param,label="nbmm")

rr_para_conf$coverage    = coverage_cal(rr_para_conf)$coverage
pval_us    = coverage_cal(rr_para_conf)$coverage
pval_deseq = coverage_cal(deseq_para_conf)$coverage
pval_nbmm  = coverage_cal(nbmm_para_conf)$coverage

rr_dd = deseq_dd = nbmm_dd = list()
##################################################################
for(i in 1:length(rr_list)){
  rr_dd[[i]]       =   dd_long(rr_list[[i]],true_param,label="with_rr")
  deseq_dd[[i]]     =   dd_long(deseq_list[[i]],true_param,label="deseq")
  nbmm_dd[[i]]      =   dd_long(nbmm_list[[i]],true_param,label="nbmm")
  ########################################################################### 
  rr_para_conf      =   para_confint(rr_dd, true_param, alpha = 0.05)
  deseq_para_conf   =   para_confint(deseq_dd, true_param, alpha = 0.05)
  nbmm_para_conf   =    para_confint(nbmm_dd, true_param, alpha = 0.05)
  ################################################################################
}

dim(rr_para_conf)







