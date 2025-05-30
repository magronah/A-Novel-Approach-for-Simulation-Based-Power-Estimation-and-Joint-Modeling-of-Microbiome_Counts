#setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
library(scam)
library(ggrastr)
library(metR)
library(latex2exp)
library(tibble)
library(ggplot2)
library(patchwork)
source("reproducible/new_sim/func.R")
path      =   paste0(getwd(),"/reproducible/new_sim/data/")
fig_path  =   paste0(getwd(),"/reproducible/new_sim/fig/")
#####################################################################
mean_count     =   readRDS(paste0(path,"/mean_count_dd.rds"))
true_param     =   readRDS(paste0(path,"/true_param_new.rds"))
pvalue_rr      =   readRDS(paste0(path,"/pval_rr.rds"))
pvalue_us      =   readRDS(paste0(path,"/pval_us.rds"))
pvalue_deseq   =   readRDS(paste0(path,"/pval_deseq2.rds"))
pvalue_nbmm    =   readRDS(paste0(path,"/pval_nbmm.rds"))
pvalue_zinbmm  =   readRDS(paste0(path,"/pval_zinbmm.rds"))
##how to get the p values for the zeroes inflation model fit
##################################################################
lmean_abund         =   rowMeans(mean_count)
lfoldchange         =   true_param$true_param
names(lmean_abund)  =   rownames(lfoldchange)
##################################################################
### pvalues
pval_rr       =   pvalue_rr$pvalue
pval_us       =   pvalue_us$pvalue
pval_deseq    =   as.numeric(rowMeans(pvalue_deseq))
pval_nbmm     =   as.numeric(rowMeans(pvalue_nbmm))
pval_zinbmm   =   as.numeric(rowMeans(pvalue_zinbmm))
##################################################################
        ### fit GAM
        fit_rr      =   gam_fit(pval_rr,lfoldchange,lmean_abund)
combined_dd_rr      =   fit_rr$combined_data
power_estimate_rr   =   fit_rr$power_estimate

         fit_us     =   gam_fit(pval_us,lfoldchange,lmean_abund)
combined_dd_us      =   fit_us$combined_data
power_estimate_us   =   fit_us$power_estimate

      fit_deseq        =   gam_fit(pval_deseq,lfoldchange,lmean_abund)
combined_dd_deseq      =   fit_deseq$combined_data
power_estimate_deseq   =   fit_deseq$power_estimate


        fit_nbmm      =   gam_fit(pval_nbmm,lfoldchange,lmean_abund)
combined_dd_nbmm      =   fit_nbmm$combined_data
power_estimate_nbmm   =   fit_nbmm$power_estimate

      fit_zinbmm        =   gam_fit(pval_zinbmm,lfoldchange,lmean_abund)
combined_dd_zinbmm      =   fit_zinbmm$combined_data
power_estimate_zinbmm   =   fit_zinbmm$power_estimate
###########################################################################
## contour plot
cont_breaks = seq(0,1,0.1); size= 25
pp_rr        =  contour_plot_fun(combined_dd_rr,power_estimate_rr,
                              cont_breaks, title_name = "reduced rank",
                              cust_theme = size)


pp_us     =  contour_plot_fun(combined_dd_us,power_estimate_us,
                              cont_breaks, title_name = "without reduced rank term",
                              cust_theme = size)

pp_deseq  =  contour_plot_fun(combined_dd_deseq,power_estimate_deseq,
                              cont_breaks, title_name = "deseq",
                              cust_theme = size)

pp_nbmm   =  contour_plot_fun(combined_dd_nbmm,power_estimate_nbmm,
                              cont_breaks, title_name = "nbmm",
                              cust_theme = size)

pp_zinbmm =  contour_plot_fun(combined_dd_zinbmm,power_estimate_zinbmm,
                              cont_breaks, title_name = "zinbmm",
                              cust_theme = size)


ppp=(pp_rr|pp_us|pp_deseq)/(pp_nbmm|pp_zinbmm|plot_spacer()) +  plot_layout(guides = "collect")

path2 = "~/Documents/PhD_Thesis/Latex_Files/RR_Model_Single_Time/doc/fig/"
ggsave(paste0(path2,"power_contour.png"), plot = ppp, width = 20, height = 12, dpi = 500)

###############################################################################
mod_obj_list  =   list(    rr  =  fit_rr$fit_2d, 
                   without_rr  =  fit_us$fit_2d, 
                       deseq   =  fit_deseq$fit_2d, 
                       nbmm    =  fit_nbmm$fit_2d,
                       zinbmm  =  fit_zinbmm$fit_2d
                   )
colour    =   c("blue","red", "green", "orange", "purple")
################################################################
mean_abun  =  c(2^3,2^4,2^5)
abs_lfc    =  seq(0.1,3,0.2)
newdata_list = pow_dd_list  = list()
for(i in 1:length(mean_abun)){
  newdata_list[[i]]  =     data.frame(lmean_abund =  rep(mean_abun[[i]], 
                                                         length(abs_lfc)),
                                      abs_lfc     =  abs_lfc)
}

names(newdata_list) = paste0("mean_abund=", mean_abun)
for(i in 1:length(mean_abun)){
  pow_dd_list[[i]]   =  power_predict_dd(mod_obj_list, newdata_list[[i]])
}

        power_dd      =   do.call(rbind,pow_dd_list)
power_dd$lmean_abund  =  factor(power_dd$lmean_abund)                                  

levels(power_dd$lmean_abund) <-  c(paste0("mean count = ", unique(power_dd$lmean_abund)))
mod_name  =   unique(power_dd$model); nnr =23


gg1 = ggplot(power_dd,aes(x = abs_lfc, y = power,group=model, 
                color = factor(model, levels = unique(model))))+
  geom_point() +
  geom_line(size = 1) + 
  scale_color_manual(values = setNames(colour, mod_name)) +  
  labs(color = "model", x= TeX("|$\\log$(fold change)|") )  +
  facet_wrap(~lmean_abund) +
  custom_theme(nnr)

ggsave(paste0(fig_path,"power_predict.png"), plot = gg1, width = 15, height = 6, dpi = 500)

path2 = "~/Documents/PhD_Thesis/Latex_Files/RR_Model_Single_Time/doc/fig/"
ggsave(paste0(path2,"power_predict.png"), plot = gg1, width = 15, height = 6, dpi = 500)
