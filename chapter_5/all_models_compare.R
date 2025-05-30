library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
source("reproducible/new_sim/func.R")
path        =   paste0(getwd(), "/reproducible/new_sim/")
##################################################################
true_param  =  readRDS(paste0(path, "10_50/true_param.rds"))
rr_est      =  readRDS(paste0(path, "10_50/parametric_rr.rds"))
nbmm_est    =  readRDS(paste0(path, "10_50/parametric_nbmm.rds"))
zinbmm_est  =  readRDS(paste0(path, "10_50/parametric_zinbmm.rds"))
deseq_est   =  readRDS(paste0(path, "10_50/parametric_deseq.rds"))
us_est      =  readRDS(paste0(path, "10_50/parametric_us.rds"))
########################################################################
true_param  =  readRDS(paste0(path, "true_param_new.rds"))
rr_est      =  readRDS(paste0(path, "rr_parametric_new.rds"))
nbmm_est    =  readRDS(paste0(path, "nbmm_parametric.rds"))
zinbmm_est  =  readRDS(paste0(path, "zinbmm_parametric.rds"))
deseq_est   =  readRDS(paste0(path, "deseq_parametric2_new.rds"))
us_est      =  readRDS(paste0(path, "us_only_parametric.rds"))
############################################################################
##convert to long format
rr_dd       =   dd_long(rr_est,true_param,label="rr")
nbmm_dd     =   dd_long(nbmm_est,true_param,label="nbmm")
zinbmm_dd   =   dd_long(zinbmm_est,true_param,label="zinbmm")
deseq_dd    =   dd_long(deseq_est,true_param,label="deseq")
us_dd       =   dd_long(us_est,true_param,label="without_rr")
#########################################################################################
rr_para_conf     =   para_confint(rr_dd, true_param, alpha = 0.05)
nbmm_para_conf   =   para_confint(nbmm_dd, true_param, alpha = 0.05)
zinbmm_para_conf =   para_confint(zinbmm_dd, true_param, alpha = 0.05)
deseq_para_conf  =   para_confint(deseq_dd, true_param, alpha = 0.05)
us_para_conf     =   para_confint(us_dd, true_param, alpha = 0.05)
#########################################################################################
## Combine average estimates with true estimates
dd       =   rbind(rr_para_conf,deseq_para_conf,nbmm_para_conf,
                   zinbmm_para_conf, us_para_conf)
dd       =   dd %>% arrange(true_param)

##### comparing average estimate with true parameter value
n  = 13
ggplot(dd, aes(x = param_name, y = average_estimate, group = model, 
                    color= factor(model, levels = c("rr", "without_rr","deseq","nbmm","zinbmm","true_param")))) +
  geom_point() +
  geom_point(aes(x = param_name, y = true_param,color="true_param")) +
  geom_line() +
  scale_color_manual(values = c("rr" = "red", "without_rr" = "pink", "deseq" = "green", 
                                "nbmm" = "orange", "zinbmm" = "purple",
                                "true_param" = "blue")) +
  labs(color = "model") +   
  xlab("taxa") +
  ylab("average estimate")  +
  custom_theme(n) +
  theme(axis.text.x = element_blank()) 
#########################################################################################
## Combine mean squared error and bias
##### Compare their MSE
rr_error      =    error_cal(rr_est, true_param, model = "rr")
nbmm_error    =    error_cal(nbmm_est, true_param, model = "nbmm")
zinbmm_error  =    error_cal(zinbmm_est, true_param, model = "zinbmm")
deseq_error   =    error_cal(deseq_est, true_param, model = "deseq")
us_error      =    error_cal(us_est, true_param, model = "without_rr")


full_summary_dd_all         =    rbind(rr_error$full_summary_dd,  
                                       us_error$full_summary_dd,
                                       nbmm_error$full_summary_dd,
                                       zinbmm_error$full_summary_dd,
                                       deseq_error$full_summary_dd)

avg_mse      =   rbind(us_error$avg_mse,nbmm_error$avg_mse, 
                       deseq_error$avg_mse, zinbmm_error$avg_mse)

avg_bias      =   rbind(us_error$avg_bias, deseq_error$avg_bias, 
                        nbmm_error$avg_bias, zinbmm_error$avg_bias)

         nr    =    4
avg_bias_rr   =  rr_error$avg_bias %>% slice(rep(1:n(), each = nr))
avg_mse_rr    =  rr_error$avg_mse %>% slice(rep(1:n(), each = nr))

avg_bias$type =  avg_bias_rr$type = paste0("bias",1:nrow(avg_bias))
avg_mse$type  =  avg_mse_rr$type = paste0("mse",1:nrow(avg_mse))

bias_dd   =  rbind(avg_bias,avg_bias_rr)
mse_dd    =  rbind(avg_mse,avg_mse_rr)
################################################################################
ggplot(full_summary_dd_all, aes(x =  param_name, y = mse, 
                                group = model, color = model)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("rr" = "red", "without_rr" = "pink", "deseq" = "green", 
                                "nbmm" = "orange", 
                                "zinbmm" = "purple",
                                "true_param" = "blue")) +
  xlab("taxa") +
  ylab("mean squared error") +
  custom_theme(n) +
  theme(axis.text.x = element_blank()) 


n=15
ppm = ggplot(mse_dd, aes(x =  type, y = average_value,group = type,  
                color = factor(model, levels = c("rr", "without_rr","nbmm","deseq","zinbmm")))) +
  geom_point(size=2) +
  scale_color_manual(values = c("rr" = "blue", "without_rr" = "red", "deseq" = "green",
                                "nbmm" = "orange", "zinbmm" = "purple")) +
  labs(color = "model") +
  ylab("average mse across taxa") +
  custom_theme(n) +
  theme(axis.title.x = element_blank())

ggsave(paste0(fig_path,"avgmse.png"), plot = ppm, width = 8, height = 6, dpi = 500)
################################################################################
ppb = ggplot(bias_dd, aes(x =  type, y = average_value,group = type,  
                   color = factor(model, levels = c("rr", "without_rr","deseq","zinbmm","nbmm")))) +
  geom_point(size=2) +
  scale_color_manual(values = c("rr" = "blue", "without_rr" = "red", "deseq" = "green",
                                "nbmm" = "orange", "zinbmm" = "purple")) +
  labs(color = "model") +
  ylab("average bias across taxa") +
  custom_theme(n) +
  theme(axis.title.x = element_blank())

ggsave(paste0(fig_path,"avgbias.png"), plot = ppb, width = 8, height = 6, dpi = 500)

################################################################################
coverage_cal = function(dd){
  df <- dd %>%
    mutate(
      coverage = if_else(true_param > lwr & true_param < upr, 1, 0)
    )
  df
}

cov1 = coverage_cal(rr_para_conf)
mean(cov1$coverage)

cov2 = coverage_cal(us_para_conf)
mean(cov2$coverage)

cov3 = coverage_cal(deseq_para_conf)
mean(cov3$coverage)

cov4 = coverage_cal(nbmm_para_conf)
mean(cov4$coverage)

cov5 = coverage_cal(zinbmm_para_conf)
mean(cov5$coverage)

################################################################################
plt  = list()
numbers <- 1:100
group_size <- 10
number_list <- split(numbers, ceiling(seq_along(numbers) / group_size))
position_dodge_width <- 0.5

for(i in 1:length(number_list)){
  elem     =  number_list[[i]]  
  dd_sub   =  dd[dd$param_name %in% paste0("taxon",elem),]
  plt[[i]] =  ggplot(dd_sub, aes(x = param_name, y = true_param, group = model, color = model)) +
    geom_pointrange(aes(ymin = lwr, ymax = upr),
                    position = position_dodge(width = position_dodge_width)) +
    scale_color_manual(values = c("rr" = "red", "without_rr" = "pink", "deseq" = "green",
                                  "nbmm" = "orange", 
                                  "zinbmm" = "purple")) +
    geom_point(position = position_dodge(width = position_dodge_width), color="black",
               size = 2.5) +  
     theme_bw() +
    theme(axis.text.x = element_blank()) +
    xlab("taxa") +
    ylab("true parameter") 
}

(plt[[1]]|plt[[2]])/(plt[[3]]|plt[[4]])/(plt[[5]]|plt[[6]])/(plt[[7]]|plt[[8]])/(plt[[9]]|plt[[10]]) +  plot_layout(guides = "collect")



dd2    =    rbind(rr_dd, nbmm_dd)
dd2    =    dd2 %>% 
            arrange(true_param)

sub       =    list(1:5, 51:55, 96:100)
#taxa_sub  =   unique(dd2$param_name)[sub] 
#dd2_sub   =   dd2[(dd2$param_name) %in% taxa_sub,]



dd_sub_list  =  plt  = list()
numbers <- 1:length(sub)
group_size <- 10
number_list <- split(numbers, ceiling(seq_along(numbers) / group_size))
position_dodge_width <- 0.5

unique(dd_sub$param_name)
for(i in 1:3){
  elem     =  sub[[i]]  
  dd_sub   =  dd2[dd2$param_name %in% paste0("taxon",elem),]
  plt[[i]] =  ggplot(dd_sub, aes(x = estimate, group = model_type, color = model_type)) +
              geom_density() +
              geom_vline(aes(xintercept = true_param), color ="red") +
              facet_wrap(~param_name, scales = "free",   ncol = 5) +
               theme_bw() 
    # xlab("taxa name") +
    # ylab("true parameter") 
}

(plt[[1]])/(plt[[2]])/(plt[[3]]) +  plot_layout(guides = "collect")





ggplot(dd2_sub, aes(x = estimate, group = model_type, color = model_type)) +
  geom_density() +
  geom_vline(aes(xintercept = true_param), color ="red") +
  facet_wrap(~param_name, scales = "free") +
  theme_bw()

  geom_point(position = position_dodge(width = position_dodge_width), color="black",
             size = 2.5) +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

###################################################################
##### Compare their MSE
rr_error       =    error_cal(rr_est, true_param, model_type = "rr")
deseq_error    =    error_cal(deseq_est, true_param, model_type = "deseq")
nbmm_error     =    error_cal(nbmm_est, true_param, model_type = "nbmm")
zinbmm_error   =    error_cal(zinbmm_est, true_param,model_type = "zinbmm")

full_summary_dd_all         =    rbind(rr_error$full_summary_dd,    
                                      deseq_error$full_summary_dd)
                         # nbmm_error$dd, zinbmm_error$dd)

single_summary_dd_all         =    rbind(rr_error$single_sum_dd,    
                                       deseq_error$single_sum_dd)

ggplot(full_summary_dd_all, aes(x =  param_name, y = mse, 
                                group = model_type, color = model_type)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("rr" = "red", "deseq" = "orange", "zinbmm" = "green", "true_param" = "blue")) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  xlab("taxa") +
  ylab("mean squared error") 



ggplot(full_summary_dd_all, aes(x =  param_name, y = bias, group = model_type, color = model_type)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("rr" = "red", "deseq" = "orange", "zinbmm" = "green", "true_param" = "blue")) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  xlab("taxa") +
  ylab("mean squared error") 



single_summary_dd_all$value  = rep("value",nrow(single_summary_dd_all))

ggplot(single_summary_dd_all, aes(x =  value, y = mse,
                                  color =model_type)) +
  geom_point() +
  #scale_color_manual(values = c("rr" = "red", "deseq" = "orange", "zinbmm" = "green", "true_param" = "blue")) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  xlab("model") +
  ylab("mean squared error") 
###################################################################
##coverage 
coverage_cal = function(dd){
  df <- dd %>%
    mutate(
      coverage = if_else(true_param > lwr & true_param < upr, 1, 0)
    )
  df
}

cov = coverage_cal(rr_para_conf)
mean(cov$coverage)
cov2 = coverage_cal(deseq_para_conf)
mean(cov2$coverage)
##################
## The reduced rank follows the trend of the true parameter than the other model

# dd             =   rbind(rr_para_conf, deseq_para_conf,
#                          nbmm_para_conf, zinbmm_para_conf)

# Calculate percentage difference with rr as the reference
rr_bias <- single_summary_dd_all$bias[single_summary_dd_all$model_type == "rr"]
single_summary_dd_all$percentage_diff <- 100 * (single_summary_dd_all$bias - rr_bias) / abs(rr_bias)

rr_para_conf    =  rr_para_conf %>% 
  mutate(param_name = factor(param_name, 
                             levels = unique(param_name[order(true_param)])))

ggplot(rr_para_conf, aes(x = param_name, y = true_param, group = model_type, color = model_type)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr)) +
  theme_bw()

#

deseq_para_conf    =  deseq_para_conf %>% 
  mutate(param_name = factor(param_name, 
                             levels = unique(param_name[order(true_param)])))

ggplot(deseq_para_conf, aes(x = param_name, y = true_param, group = model_type, color = model_type)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr)) +
  theme_bw()
################################################################################
plot(density(as.numeric(deseq_est[1,]), probability = T))

#We now check if  I can fit t distribution to make the wald confidence equal to 
## to the parameter confidence interval 
#############################################
library(fitdistrplus)
fit = fitdistr(as.numeric(deseq_est[1,]), "t", df = 10)
m  =  (fit$estimate)

### power calculation
## simulate with new effects each time, no I cant do that? Keep the same effect
## effects
##' 
