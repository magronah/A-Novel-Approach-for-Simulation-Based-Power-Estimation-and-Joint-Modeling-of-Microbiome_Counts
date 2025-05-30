library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
source("reproducible/new_sim/func.R")
path        =   paste0(getwd(), "/reproducible/new_sim/data/")
########################################################################################
true_param      =  readRDS(paste0(path, "true_param_new.rds"))
rr_est          =  readRDS(paste0(path, "rr_parametric_new.rds"))
rr_zi_est       =  readRDS(paste0(path, "rr_zi_parametric.rds"))
rr_zitaxa_est   =  readRDS(paste0(path, "rr_zitaxa_parametric.rds"))
#########################################################################################
##convert to long format
rr_dd          =   dd_long(rr_est,true_param,label="rr_without_zi")
rr_zi_dd       =   dd_long(rr_zi_est,true_param,label="rr_single_zi")
rr_zitaxa_dd   =   dd_long(rr_zitaxa_est,true_param,label="rr_zi_per_taxa")
#########################################################################################
rr_para_conf        =   para_confint(rr_dd, true_param, alpha = 0.05)
rr_zi_para_conf     =   para_confint(rr_zi_dd, true_param, alpha = 0.05)
rr_zitaxa_para_conf =   para_confint(rr_zitaxa_dd, true_param, alpha = 0.05)
#########################################################################################
## Combine average estimates with true estimates
dd       =   rbind(rr_para_conf,rr_zi_para_conf, rr_zitaxa_para_conf)
dd       =   dd %>% arrange(true_param)

##### comparing average estimate with true parameter value
n  = 18
g1 = ggplot(dd, aes(x = param_name, y = average_estimate, group = model, 
                    color= factor(model, levels = c("rr_without_zi", "rr_single_zi",
                                                    "rr_zi_per_taxa","true_param")))) +
  geom_point() +
  geom_point(aes(x = param_name, y = true_param,color="true_param")) +
  geom_line() +
  scale_color_manual(values = c("rr_without_zi" = "red", "rr_single_zi" = "black",
                                "rr_zi_per_taxa" = "green", "true_param" = "blue")) +
  labs(color = "model") +   
  xlab("taxa") +
  ylab("average estimate")  +
  custom_theme(n) +
  theme(axis.text.x = element_blank()) 


#ggsave(paste0(fig_path,"rr_est.png"), plot = g1, width = 8, height = 6, dpi = 500)
###################################################################
##### Compare their MSE
rr_error          =   error_cal(rr_est, true_param, model = "rr_without_zi")
rr_zi_error       =   error_cal(rr_zi_est, true_param, model = "rr_single_zi")
rr_zitaxa_error   =   error_cal(rr_zitaxa_est, true_param, model = "rr_zi_per_taxa")

full_summary_dd_all         =    rbind(rr_error$full_summary_dd,    
                                       rr_zi_error$full_summary_dd,
                                       rr_zitaxa_error$full_summary_dd)

single_summary_dd_all       =    rbind(rr_error$single_sum_dd,    
                                       rr_zi_error$single_sum_dd,
                                       rr_zitaxa_error$single_sum_dd)
###################################################################
g2= ggplot(full_summary_dd_all, aes(x =  param_name, y = mse, 
                                    group = model, color = model)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("rr_without_zi" = "red", "rr_single_zi" = "black",
                                "rr_zi_per_taxa" = "green")) +
  xlab("taxa") +
  ylab("mean squared error") +
  custom_theme(n) +
  theme(axis.text.x = element_blank()) 

#ggsave(paste0(fig_path,"us_rr_mse.png"), plot = g2, width = 8, height = 6, dpi = 500)
###################################################################
g3 = ggplot(full_summary_dd_all, aes(x =  param_name, y = bias, group = model, color = model)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("rr_without_zi" = "red", "rr_single_zi" = "black",
                                "rr_zi_per_taxa" = "green")) +
  xlab("taxa") +
  ylab("bias") + 
  custom_theme(n) +
  theme(axis.text.x = element_blank()) 

#ggsave(paste0(fig_path,"us_rr_bias.png"), plot = g3, width = 8, height = 6, dpi = 500)
###################################################################
g4=ggplot(single_summary_dd_all, aes(x =  type, y = average_value,
                                     group = type,  color =model)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("rr_without_zi" = "red", "rr_single_zi" = "black",
                                "rr_zi_per_taxa" = "green")) +
  ylab("average value across taxa") +
  custom_theme(n) +
  theme(axis.title.x = element_blank())

#The model's performance for both zero inflations is the same. 

#ggsave(paste0(fig_path,"us_rr_avgmse.png"), plot = g4, width = 8, height = 6, dpi = 500)
###################################################################
##confidence intervals
plt  = list(); numbers <- 1:100; group_size <- 10

number_list <- split(numbers, ceiling(seq_along(numbers) / group_size))
position_dodge_width <- 0.5

for(i in 1:length(number_list)){
  elem     =  number_list[[i]]  
  dd_sub   =  dd[dd$param_name %in% paste0("taxon",elem),]
  plt[[i]] =  ggplot(dd_sub, aes(x = param_name, y = true_param, group = model, color = model)) +
    geom_pointrange(aes(ymin = lwr, ymax = upr),
                    position = position_dodge(width = position_dodge_width)) +
    geom_point(position = position_dodge(width = position_dodge_width), color="black",
               size = 2.5) +  
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank()) +  # Remove individual x labels
    xlab("taxa") +
    ylab("true parameter") 
}


# Combine the plots and add a single x-axis label at the bottom
combined_plot <- (plt[[1]] | plt[[2]]) / (plt[[3]] | plt[[4]]) / 
  (plt[[5]] | plt[[6]]) / (plt[[7]] | plt[[8]]) / 
  (plt[[9]] | plt[[10]]) +
  plot_layout(guides = "collect") +
  theme(plot.margin = unit(c(1, 1, 2, 1), "cm"))  # Adjust bottom margin

# Add the x-axis label "taxa" at the bottom
combined_plot + plot_annotation(caption = "taxa") & 
  theme(plot.caption = element_text(hjust = 0.5, vjust = 1, size = 12, face = "bold"))
(plt[[1]]|plt[[2]])/(plt[[3]]|plt[[4]])/(plt[[5]]|plt[[6]])/(plt[[7]]|plt[[8]])/(plt[[9]]|plt[[10]]) +  plot_layout(guides = "collect")

#################################################################################
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
cov2 = coverage_cal(us_para_conf)
mean(cov2$coverage)
###################################################################

