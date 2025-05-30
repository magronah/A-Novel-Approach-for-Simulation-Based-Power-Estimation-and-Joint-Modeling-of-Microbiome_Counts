setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
source("reproducible/simp_single/Load_Packages.R")
source("reproducible/simp_single/func.R")
source("reproducible/simp_single/param.R")
source("reproducible/simp_single/initial_param0.R")
source("load_glmmTMB.R")

path = paste0(getwd(), "/reproducible/simp_single/data/")
###########################################################################
true_param      =    readRDS(paste0(path,"true_param.rds"))
names(true_param) = c("param_name","true_param")
true_param_all  =    readRDS(paste0(path,"true_all_b.rds"))
###########################################################################
##sanity check
stopifnot(true_param_all[seq(2,(2*ntaxa),2)] == true_param$true_param)

ref_mod       =    readRDS(paste0(path,"glmmTMB_101.rds"))
rr_est        =    readRDS(paste0(path,"est_glmmTMB.rds"))

wald          =   wald_confint(ref_mod, ntaxa, conf_level = .95)
wald_dd       =   left_join(wald, true_param, by = "param_name")

dd_long1        =   dd_long(rr_est,label="rr")
dd_para_conf1   =   para_confint(dd_long1, true_param, alpha = 0.05)
dd_para_conf11  =   dd_para_conf1 %>% dplyr::rename(est_param = average_estimate)
##################################################################################
p1 =  ggplot(dd_para_conf11, aes(x = param_name, y = est_param)) +
          geom_pointrange(aes(ymin = lwr, ymax = upr),size =0.2) +
  labs(y ="average_param") +
  theme(
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    axis.text.x = element_blank()   # Optionally remove x-axis text labels
  )

p2 =  ggplot(dd_para_conf11, aes(x = param_name, y = true_param)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr), size =0.2) + theme(
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    axis.text.x = element_blank()   # Optionally remove x-axis text labels
  )

p1|p2
#################################################################################

dd_full      =    rbind(dd_para_conf11, wald_dd)

dd_sub_list  =  plt  = list()
numbers <- 1:50
group_size <- 10
number_list <- split(numbers, ceiling(seq_along(numbers) / group_size))
position_dodge_width <- 0.5

for(i in 1:length(number_list)){
  
  elem     =  number_list[[i]]  
  dd_sub   =  dd_full[dd_full$param_name %in% paste0("taxon",elem),]
  plt[[i]] =  ggplot(dd_sub, aes(x = param_name, y = true_param, group = model_type, color = model_type)) +
    geom_pointrange(aes(ymin = lwr, ymax = upr),
                    position = position_dodge(width = position_dodge_width)) +
    geom_point(position = position_dodge(width = position_dodge_width), color="black",
               size = 3) 
}


(plt[[1]]|plt[[2]])/(plt[[3]]|plt[[4]]) +  plot_layout(guides = "collect")
