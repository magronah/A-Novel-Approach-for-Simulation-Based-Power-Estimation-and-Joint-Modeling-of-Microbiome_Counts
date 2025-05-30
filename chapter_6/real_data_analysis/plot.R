library(glmmTMB)
library(ggplot2)
library(dplyr)
library(patchwork)
library(nlme)
library(MASS)
library(Matrix)
library(NBZIMM)
library(here)
##############################################################
fig_path = "fig/"
source("func.R")
#########################################################
group_label = "bmi_groupobese:time"
# Load datasets
atlass_path    <- paste0("real_data_analysis/atlass/results/")
pregnancy_path <- paste0("real_data_analysis/pregnancy/results/")
############################################################
#confidence width






############################################################

# Define filenames to load
runtime_files <- c(
  "runtime_rr.rds",
  "runtime_us.rds",
  "runtime_nb.rds",
  "runtime_znb.rds")

# Load models 
atlass_runtimes    <- load_models(atlass_path, runtime_files)
pregnancy_runtimes <- load_models(pregnancy_path, runtime_files)

names(atlass_runtimes)     =    c("RR","US","NB","ZINB")
names(pregnancy_runtimes)  =    c("RR","US","NB","ZINB")

time_atlass <-  data.frame(time = unlist(lapply(atlass_runtimes, 
                                            function(x){x[["elapsed"]]/60})), 
                       model =  names(atlass_runtimes)) %>%
                mutate(dataset = "The Human Intestine Data")

time_pregnancy <-  data.frame(time = unlist(lapply(pregnancy_runtimes, 
                                            function(x){x[["elapsed"]]/60})), 
                       model =  names(pregnancy_runtimes)) %>%
                   mutate(dataset = "Pregnancy Data")

time_dd  <-  rbind(time_atlass, time_pregnancy)

plt_time <- ggplot(time_dd, aes(x =  model, y =  time)) + 
           geom_point(size=2) +
  custom_theme(16) +
  labs(x = " ", y = "Runtime (mins)") +
  facet_wrap(~dataset, scales = "free") 

width =  10; height =  5; dpi = 300
ggsave("fig/runtime.png", plot = plt_time, 
       width = width, 
       height = height, 
       dpi = dpi)
############################################################
conf_level <-  0.95; in_sd  = 6.7
z_score    <-  qnorm(conf_level + (1 - conf_level)/2)
# statistical power
#Confidence Intervals 
confint_files <- c(
  "confint_rr.rds",
  "confint_us.rds",
  "confint_nb.rds",
  "confint_znb.rds"
)
##### Load 
atlass_conf    <- load_models(atlass_path, confint_files)
pregnancy_conf <- load_models(pregnancy_path, confint_files)

names(atlass_conf)     =  c("RR", "US", "NB", "ZINB")  
names(pregnancy_conf)  =  c("RR", "US", "NB", "ZINB") 
############################################################
atlass_confint_sub      <-  lapply(atlass_conf, filter_complete_separation)
atlass_conft  <-  lapply(atlass_confint_sub, function(x){ 
  df  =  x %>% filter(group_label == "bmi_groupobese:time")
  df$group_label <- NULL
  df})

atlass_confint  <-  lapply(atlass_conft, function(x){
                     taxa_include =  atlass_conft$NB$param_name
                     x   %>%  filter(param_name %in% taxa_include)})

atlass_confint[c("US", "RR")] <- lapply(atlass_confint[c("US", "RR")], function(x) {
  x$lwr_old     =   x$lwr
  x$upr_old     =   x$lwr
  x$sd_err      =   (x$est_param - x$lwr) / z_score
  x$lwr         =   x$est_param - z_score*in_sd*x$sd_err
  x$upr         =   x$est_param + z_score*in_sd*x$sd_err
  x
})

sig  =  0.05
atlass_confint <-  lapply(atlass_confint, function(x){
                  x$padjust  =  p.adjust(x$pvalue, method = "BH")
                  x$pow      =  ifelse(x$padjust < sig, 1, 0)
                  x})
########################################################
sub  =  names(atlass_confint$NB)
Atlass_confint <- lapply(atlass_confint, function(x) {
         x %>% dplyr::select(all_of(sub))})
############################################################
pregnancy_confint_sub   <-  lapply(pregnancy_conf,filter_complete_separation)
pregnancy_confint  <-  lapply(pregnancy_confint_sub, function(x){
               taxa_include =  pregnancy_confint_sub$NB$param_name
               x   %>%  filter(param_name %in% taxa_include)})

pregnancy_confint <-  lapply(pregnancy_confint, function(x){
                      x$padjust  =  p.adjust(x$pvalue, method = "BH")
                      x$pow      =  ifelse(x$padjust < sig, 1, 0)
                       x})

pregnancy_confint[c("US", "RR")] <- lapply(pregnancy_confint[c("US", "RR")], function(x) {
  x$lwr_old     =   x$lwr
  x$upr_old     =   x$lwr
  x$sd_err      =   (x$est_param - x$lwr) / z_score
  x$lwr         =   x$est_param - z_score*in_sd*x$sd_err
  x$upr         =   x$est_param + z_score*in_sd*x$sd_err
  x
})
########################################################
sub1  =  names(pregnancy_confint$NB)
Pregnancy_confint <- lapply(pregnancy_confint, function(x) {
  x %>% dplyr::select(all_of(sub1))})
############################################################
atlass_combined    =  bind_rows(Atlass_confint, .id = "model")
pregnancy_combined =  bind_rows(Pregnancy_confint, .id = "model")
#################################################################
atlass_combined$model_type =ifelse(atlass_combined$model %in% c("RR", "US"), 
                                   "RR and US","NB and ZINB")

pregnancy_combined$model_type =ifelse(pregnancy_combined$model %in% c("RR", "US"), 
                                   "RR and US","NB and ZINB")
#################################################################
## Plot confidence intervals

nn = 18; wd = 0.4; sz =  0.1
okabe_ito_palette <- c("orange",  "violet", "black","blue")

pp1 <- ggplot(atlass_combined, aes(x = log(mean_count), y = width, color = model)) +
  geom_point() +
  geom_line() +
  labs(
    x = "log(mean abundance)",
    y = "Confidence width",
    title =  "Human intestine data",
    color = "Model") +
  custom_theme(nn) + 
  xlim(range(log(atlass_combined$mean_count))) +
  scale_color_manual(values = okabe_ito_palette)    +   
  facet_wrap(~model_type) 

pp2 <-  ggplot(pregnancy_combined, aes(x = log(mean_count), y = width, color = model)) +
  geom_point() +
  geom_line() +
  labs(
    x = "log(mean abundance)",
    y = "Confidence width",
    title =  "Pregnancy data",
    color = "Model") +
  custom_theme(nn) + 
  xlim(c(-5,7.3)) +
  scale_color_manual(values = okabe_ito_palette)    +   
  facet_wrap(~model_type) 

pp12 = (pp1/pp2) +   plot_layout(guides = "collect")
width =  10; height =  9; dpi = 300
ggsave("fig/confwidth_rd.png", plot = pp12, 
       width = width, 
       height = height, 
       dpi = dpi)
###########################################################################
# 
pregnancy_power <- pregnancy_combined %>%
  group_by(model) %>%
  summarise(power =  mean(pow), 
            x     =  sum(pow),
            n     =  length(pow)) %>%
  mutate(dataset = " Pregnancy Data") %>%
  data.frame()

sd  =  sqrt(pregnancy_power$power*(1 - pregnancy_power$power)/pregnancy_power$n)
pregnancy_power$lwr  = pregnancy_power$power - z_score*sd 
pregnancy_power$upr  = pregnancy_power$power + z_score*sd 

atlass_power <- atlass_combined %>%
  group_by(model) %>%
  summarise(power =  mean(pow), 
            x     =  sum(pow),
            n     =  length(pow)) %>%
  mutate(dataset = " Human Intestinal Data") %>%
  data.frame()

sd  =  sqrt(atlass_power$power*(1 - atlass_power$power)/atlass_power$n)
atlass_power$lwr  = atlass_power$power - z_score*sd 
atlass_power$upr  = atlass_power$power + z_score*sd 

power_combind =  rbind(atlass_power, pregnancy_power)
###########################################################################


plt1 <- ggplot(atlass_combined, aes(x = log(mean_count), y = est_param, 
                                    ymin = lwr, ymax = upr, color = model)) +
               geom_pointrange(position = position_dodge(width = wd), 
                               size = sz) +
             labs(
               x = "log(mean abundance)",
               y = "95% Confidence Interval",
               color = "Model") +
              custom_theme(nn) + 
              #scale_color_manual(values = okabe_ito_palette) +   
              facet_wrap(~model_type) 

  

plt2 <- ggplot(pregnancy_combined, aes(x = log(mean_count), y = est_param, 
                                    ymin = lwr, ymax = upr, color = model)) +
               geom_pointrange(position = position_dodge(width = wd), 
                               size = sz) +
             labs(
               x = "log(mean abundance)",
               y = "95% Confidence Interval",
               color = "Model") +
              custom_theme(nn) + 
              scale_color_manual(values = okabe_ito_palette) +   
              facet_wrap(~model_type) 


width =  10; height =  5; dpi = 300
ggsave("fig/atlass.conf_int_side.png", plot = plt1, 
       width = width, 
       height = height, 
       dpi = dpi)

ggsave("fig/pregnancy.conf_int_side.png", plot = plt2, 
       width = width, 
       height = height, 
       dpi = dpi)
############################################################
library(dplyr)
library(purrr)

compute_power_ci <- function(df, dataset_name) {
  df %>%
    group_by(model) %>%
    summarise(
      power = mean(pow),
      x     = sum(pow),
      n     = length(pow),
      .groups = "drop"
    ) %>%
    mutate(
      dataset = dataset_name,
      ci      = map2(x, n, ~ stats::binom.test(.x, .y)$conf.int),
      lwr     = map_dbl(ci, 1),
      upr     = map_dbl(ci, 2)
    ) %>%
   dplyr::select(-ci)  
}

View(atlass_combined)

atlass_power    <- compute_power_ci(atlass_combined, " Human Intestinal Data")
pregnancy_power <- compute_power_ci(pregnancy_combined, "Pregnancy Data")
power_dd   <-    rbind(atlass_power, pregnancy_power)
power_dd$dataset  =   factor(power_dd$dataset, levels = unique(power_dd$dataset))

pow_plt = ggplot(power_combind, aes(x =  model, y = power, ymin = lwr, ymax = upr)) +
  geom_pointrange(position = position_dodge(width = 1), size = 1) +
  facet_wrap(~dataset) +
  labs(x = " ", y = "average statistical power across taxa") +
  custom_theme(nn) 
  
width =  10; height =  5; dpi = 300
ggsave("fig/power_long.png", plot = pow_plt, 
       width = width, 
       height = height, 
       dpi = dpi)
