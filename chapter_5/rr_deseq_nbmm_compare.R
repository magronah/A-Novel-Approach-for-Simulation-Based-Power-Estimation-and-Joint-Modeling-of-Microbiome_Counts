library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(here)
source("func2.R")
fig_path=   "fig/"
#####################################################
####Bias and RMSE calculation
strg    =   c("100_300/","150_500/", "200_600/")
titles  =   c("50 subjects per group and 300 taxa",
             "75 subjects per group and 500 taxa",
             "100 subjects per group and 600 taxa")
########################################################
n = 20; p1  =   p2   =  p3 =  p4  = p5  = p6 = p7 = list()
BIAS = VARIAN  =  list()
sz =1; wd  = 0.3; lw = 1
for(i in 1:length(strg)){
  ########################################################
  path  =   paste0(strg[[i]]) 
  dd    =   load_data(path) 

  conf =  dd$confint
  ddd  = do.call(rbind,conf)
  ddd  =  ddd %>%
          filter(ddd$model != "deseq_NS" )
  ########################################################
 p1[[i]] = ggplot(ddd, aes(x = param_name, y = average_estimate, color = model, group = model)) +
    geom_point(size = 1, alpha = 0.5) +  
    geom_line() +   
    geom_point(aes(x = param_name, y = true_param, color = "true group effect"), 
               size = 1.5, shape = 17) +   
    ylim(-5, 5)  +# Set y-axis limits
    # scale_color_manual(
    #   values = c("blue", "red", "green", "black", "purple", "orange", "cyan"),  
    #   name = "Model"
    # ) +  
    labs(title = titles[[i]],
         x = "Taxa",   
         y = "Average group effect estimates",
         color = "Model") +  
    theme_bw(base_size = n) +
    theme(
      plot.title = element_text(hjust = 0.5, size = n),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),  
      axis.ticks.x = element_blank(),
      text = element_text(size = n, family = "Roboto")
    )  
  #####################################################
  mse   =   err_extract(dd, "avg_mse")
  bias  =   err_extract(dd, "avg_bias")
  var   =   err_extract(dd, "avg_var")
  
  p2[[i]] = ggplot(mse, aes(model, average_value)) +
    geom_point(size = 3*sz) + 
    geom_errorbar(aes(ymin = lwr, ymax = upr), size=sz,  width = wd) +
    geom_hline(yintercept = mse["RRzi",1],  linewidth = lw,linetype = "dashed", color = "red") +
    custom_theme(n) +
    labs(title= titles[[i]],
         y = "Average RMSE across taxa",
         x = " "
         )
  
  BIAS[[i]] = bias
  p3[[i]] = ggplot(bias, aes(model, average_value)) +
    geom_point(size = 3*sz) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), size=sz, width = wd) +
    geom_hline(yintercept = bias["RRzi",1],  linewidth = lw,linetype = "dashed", color = "red") +
    custom_theme(n) +
    labs(title = titles[[i]],
         y = "Average bias across taxa"
           #"Comparison of average bias across taxa",x = " ",y = "Average Bias"
           )
  
  
  VARIAN[[i]] = var
  var$type   =   if_else(var$model %in% c("NB","ZNB"), "NB and ZNB", 
                         "RR, RRzi US, USzi, DE and DE_noSh")
  
  p4[[i]] = ggplot(var, aes(model, average_value)) +
    geom_point(size = 3*sz) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), size=sz, width = wd) +
    geom_hline(yintercept = (var["RRzi",1]), linewidth = lw, linetype = "dashed", color = "red") +
    custom_theme(n) +
    labs(title = titles[[i]],
         x = " ",
         y = "Average variance"
         # y = "Average variance of error across taxa"
         #"Comparison of average bias across taxa",x = " ",y = "Average Bias"
    ) +
    facet_wrap(~type, scale ="free")
  p4
  #####################################################
  ####Coverage Calculation 
  pp       =   do.call(rbind,dd[["confint"]]) %>%  
    arrange(true_param)  %>%  
    as.data.frame(row.names = NULL)
  
  # Compute the mean of lwr and upr by model type
  mean_values <- aggregate(cbind(lwr, upr, CI_width,true_param) ~ model, 
                           data = pp, FUN = mean)
  
  # Create the plot
  p5[[i]]= ggplot(mean_values, aes(x = model, ymin = lwr, ymax = upr, y = (lwr + upr) / 2)) +
    geom_pointrange(size = 0.8) +
    geom_hline(aes(yintercept = true_param), linetype = "dashed", color = "red", size = 1) +
    labs(
      title = "Confidence Intervals for Models",
      x = "Model",
      y = "Estimate",
      color = "Model"
    ) +
    custom_theme(n) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_flip()
  
  
  p6[[i]]  = ggplot(mean_values, aes(x = model, y = CI_width)) +
    geom_point() +
    labs(
      title = "Confidence Width for Models",
      x = "Model",
      y = "Confidence Width",
    ) +
    custom_theme(n)
  
  num_taxa  = nrow(dd$dd$DE)
  #####################################################
  # Sum coverage by model type
  coverage_dd <- aggregate(coverage ~ model, data = pp, sum)
  coverage_dd$coverage =  coverage_dd$coverage/num_taxa  
  
  
  p7[[i]]  = ggplot(coverage_dd, aes(x = model, y = coverage)) +
    geom_point() +
    custom_theme(n)
}
##################################################################
trend <- (p1[[1]]|p1[[2]]|p1[[3]]) + plot_layout(guides = "collect")  
size = 3; width =  28; height =  8; dpi = 300
ggsave("fig/trend2.png", plot = trend, 
       width = width, 
       height = height, 
       dpi = dpi)

rmse <- (p2[[1]]|p2[[2]]|p2[[3]]) + plot_layout(guides = "collect")  
ggsave("fig/rmse.png", plot = rmse, 
       width = width, 
       height = height, 
       dpi = dpi)

bias <- (p3[[1]]|p3[[2]]|p3[[3]]) + plot_layout(guides = "collect")  
ggsave("fig/bias2.png", plot = bias, 
       width = width, 
       height = height, 
       dpi = dpi)
size = 3; width =  15; height =  12; dpi = 300

var_p <- (p4[[1]]/p4[[2]]/p4[[3]]) + plot_layout(guides = "collect")  
ggsave("fig/var_plt_3.png", plot = var_p, 
       width = width, 
       height = height, 
       dpi = dpi)
##################################################################
filenames  <-  c("us.rds", "uszi.rds", "rr.rds")
path       =   paste0(strg[[3]])
dd         =   load_models(paste0(path, "confint/"),filenames)
names(dd)  =   c("us", "uszi", "rr")     
##################################################################
true_param    =   readRDS(paste0(path, "true_param.rds"))
##################################################################
alpha  <-  0.05
ddd <- lapply(dd, function(x){
  ddf <- lapply(x, function(y){
    y$p.adj    =   p.adjust(y$pvalue, "BH")
       df      =   left_join(true_param,y, by = "param_name")
    df$cov     =   ifelse(df$lwr  < df$true_param & df$true_param < df$upr,1, 0)
    df$pval_reject  =   ifelse(df$p.adj  < alpha,1, 0)
    df
  })
  names(ddf)  = names(x)
  ddf
})
##################################################################
library(dplyr)
library(purrr)
cov_dd  <-   lapply(ddd, function(x){
  p = x %>%  map_dfc(~ pull(.x, cov))
  mean(rowMeans(p))
})

covv = data.frame(coverage =do.call("rbind", cov_dd)) %>%
       rownames_to_column("model")
##################################################################
size = 3; nn= 14; width =  7; height =  5; dpi = 300
plt0 = ggplot(covv, aes(model, coverage)) + 
         geom_point(size = size)  +  
         custom_theme(nn) +
        labs(x = "")

ggsave("fig/coverage.png", plot = plt0, 
       width = width, 
       height = height, 
       dpi = dpi)
##################################################################
pval_dd  <-   lapply(ddd, function(x){
  p = x %>%  map_dfc(~ pull(.x, pval_reject))
  mean(colMeans(p))
})

pval_dd_ = data.frame(pow = do.call("rbind", pval_dd))  %>%
  rownames_to_column("model")
 
plt1 = ggplot(pval_dd_, aes(model, pow)) + 
       geom_point(size = size)+
       custom_theme(nn)  +
       labs(x = "", y ="power")

ggsave("fig/power.png", plot = plt1, 
       width = width, 
       height = height, 
       dpi = dpi)

plt0 + plt1
####################################################################

                  
########################################################
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(here)
source("func2.R")
fig_path=   "fig/"
#####################################################
####Bias and RMSE calculation
strg    =   c("100_300/","150_500/", "200_600/")
titles  =   c("50 subjects per group and 300 taxa",
             "75 subjects per group and 500 taxa",
             "100 subjects per group and 600 taxa")
########################################################
# Define filenames to load
filenames <- c(
  "rr_aicc.rds",
  "rrzi_aicc.rds",
  "us_aicc.rds",
  "uszi_aicc.rds",
  "nbmm_aicc.rds",
  "zinbmm_aicc.rds",
  "deseq_aicc.rds"  
)



# Load models for autism data
aicc <- load_models(strg[2], filenames)
# Assigning names 
aicc_names <- c("RR","RRzi","US","USzi","NB","ZNB","DE")
names(aicc)    =  aicc_names
# names(crohn_models)    =   names(soil_models)    =    mod_names
# mod_list     =   lst(autism_models,atlass_models,crohn_models,soil_models)

sim_list <- lapply(paste0("sim", 1:10), function(sim) {
  sapply(aicc, function(x) which(x == sim))
})

names(sim_list) <- paste0("sim", 1:10)


v   =   as.data.frame(do.call(rbind,aicc))

for(i in 1:10){
  vv  =   v %>% dplyr::select(paste0("sim",i)) %>%
         rownames_to_column("model")
  vv[[paste0("sim", i)]] <- as.numeric(vv[[paste0("sim", i)]])  
  vv$aicc_diff <- vv[[paste0("sim", i)]] - min(vv[[paste0("sim", i)]], na.rm = TRUE)
  vvv  =vv %>% arrange(aicc_diff)

  print(vvv)
}


View(v)
names(sim_list) <- paste0("sim", 1:10)

sim_dd  <-   t(data.frame(lapply(aicc, function(x){
              x[["sim1"]]
})))


class(aicc)
lapply(aicc, names)
  
                  
height =  5; width = 17
##################################################################
RMSE = (p2[[1]]|p2[[2]]|p2[[3]]) +  plot_layout(guides = "collect") 

ggsave("fig/mse.png", plot = RMSE, width = width, 
        height = height, dpi = 300)
##################################################################
BIAS = (p3[[1]]|p3[[2]]|p3[[3]]) +  plot_layout(guides = "collect") 

ggsave("fig/bias.png", plot = BIAS, width = width,
       height = height, dpi = 300)
##################################################################
VAR = (p4[[1]]|p4[[2]]|p4[[3]]) +  plot_layout(guides = "collect") 

ggsave("fig/var.png", plot = VAR, width = width,
       height = height, dpi = 300)
##################################################################
(p4[[1]]|p4[[2]]|p4[[3]]) +  plot_layout(guides = "collect") 

(p5[[1]]|p5[[2]]|p5[[3]]) +  plot_layout(guides = "collect") 

(p6[[1]]|p6[[2]]|p6[[3]]) +  plot_layout(guides = "collect") 
##################################################################
####Confidence Interval Width
#mm   =                      
ggplot(mm, aes(x = model, y = true_param)) +
  # Add confidence interval lines
  geom_linerange(aes(ymin = lwr, ymax = upr), color = "blue", size = 1.5) +
  # Add true parameter line for reference
  geom_hline(aes(yintercept = true_param), linetype = "dashed", color = "red") +
  # Add coverage as points, colored based on coverage value
  geom_point(aes(color = coverage), size = 4) +
  # Set colors for coverage
  scale_color_gradient(low = "red", high = "green") +
  # Rotate x-axis labels for better readability
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # Add labels and title
  labs(
    title = "Model Confidence Intervals and Coverage",
    x = "Model",
    y = "True Parameter"
  ) +
  theme_minimal()






#####################################################
path1   =   paste0(getwd(),"/50_200_previous_500sim/")
path1   =   paste0(getwd(),"/50_200/")
path2   =   paste0(getwd(),"/100_300/")
path3   =   paste0(getwd(),"/150_500/")
#####################################################
dd1     =   load_data(path1) 
dd2     =   load_data(path2) 
dd3     =   load_data(path3) 
#####################################################
pp1     =   do.call(rbind,dd1[["confint"]]) %>% 
              arrange(true_param)
pp1$width =  pp1$upr - pp1$lwr

pp2     =   do.call(rbind,dd2[["confint"]]) %>% 
              arrange(true_param)
pp3     =   do.call(rbind,dd3[["confint"]]) %>% 
              arrange(true_param)

names(pp1)
mod <- c("rr","nbmm")
pp11  = pp1  %>%
        mutate(lwr2 = lwr - average_estimate,
               upr2 = upr - average_estimate)

ggplot(pp11 %>% filter(model %in% mod) , aes(y= param_name, color = model)) +
  #geom_text(aes(x=upr,label=upr,color=model),size=2, hjust=-0.4,vjust=-0.9, show.legend = FALSE)+
  #geom_text(aes(x=lwr,label=lwr,color=model),size=2, hjust=1.2,vjust=-0.9, show.legend = FALSE)+
  geom_errorbarh(aes(xmin=lwr2,xmax=upr2,color=model),height=4)+
  facet_wrap(~model) +
  #geom_point() +
  #geom_line() +
  theme_bw()
###################################################################
# Create plots in a loop
plot_list <- list()

for (i in seq(1, 200, by = 10)) {
  taxon_range <- paste0("taxon", i:(i + 9))
  plot_list[[length(plot_list) + 1]] <- ggplot(pp1 %>% filter(param_name %in% taxon_range), 
                                               aes(x = true_param, y = param_name, color = model)) +
    geom_point(position = position_dodge(width = 0.5)) + 
    coord_flip() +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), position = position_dodge(width = 0.5), height = 0.2) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(x = "Estimate", title = paste("Parameters", i, "-", i + 9))
}

(plot_list[[1]] | plot_list[[2]]) / 
(plot_list[[3]] | plot_list[[4]]) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
############################################################
bias1    =    err_extract(dd1, "avg_bias")
bias11   =    reorganise_dd(bias1, "bias")
bias11$label  = factor("50 subjects and 200 ntaxa")

bias2   =   err_extract(dd2, "avg_bias")
bias22   =   reorganise_dd(bias2, "bias")
bias22$label  = factor("100 subjects and 300 ntaxa")

bias3   =   err_extract(dd3, "avg_bias")
bias33   =   reorganise_dd(bias3, "bias")
bias33$label  = factor("150 subjects and 500 ntaxa")
###################################################
mse1    =   err_extract(dd1, "avg_mse")

mse11   =   reorganise_dd(mse1, "mse")
mse11$label  = factor("50 subjects and 200 ntaxa")

mse2   =   err_extract(dd2, "avg_mse")
mse22   =   reorganise_dd(mse2, "mse")
mse22$label  = factor("100 subjects and 300 ntaxa")

mse3   =   err_extract(dd3, "avg_mse")
mse33   =   reorganise_dd(mse3, "mse")
mse33$label  = factor("150 subjects and 500 ntaxa")
###################################################
mse_dd  =  rbind(mse11,mse22,mse33)
bias_dd =  rbind(bias11,bias22,bias33)
##################################################
n = 11
oka_col = c(
  "#0000FF",
  "#556B2F", 
  "#E23D28", 
  "#E69F00", 
  "#000000",
  "#56B4E9", 
  "#D55E00"  # Vermilion
  #"#56B4E9", # Sky Blue
  #"#009E73", # Bluish Green
  #"#F0E442", # Yellow
  #"#0072B2", # Blue
)

mse222 = mse2[mse2$model %in% c("rrzi", "deseq", "zinbmm", "nbmm"), ]
rrzi_value <- mse222$average_value[mse222$model == "rrzi"]

mse222$percentage_change <- ((mse222$average_value - rrzi_value) / rrzi_value) * 100

mse222$model <- factor(mse222$model, levels = c("rrzi", "deseq", "zinbmm", "nbmm"))

g1 <- ggplot(mse222, aes(x = model, y = average_value, color = model)) +
  geom_point(size = 3) +  # Increase point size
  geom_hline(yintercept = rrzi_value, linetype = "dashed", color = "black", 
             size= 1.2) +  # Increase line thickness
  scale_color_manual(values = oka_col) +
  ylab("average value across taxa") +
  custom_theme(n) +
  labs(color = "model") +
  theme(axis.title.x = element_blank())

g1




g1=ggplot(mse33, aes(x =  type, y = average_value, 
                         color = factor(model, levels = unique(model)))) +
  geom_point() +
  scale_color_manual(values = oka_col) +
  ylab("average value across taxa") +
  custom_theme(n) +
  labs(color = "model")+
  #facet_wrap(~label, scales = "free") +
  theme(axis.title.x = element_blank())
g1

dim(dd1$dd$rr)
x11()



g2 =  ggplot(bias_dd, aes(x =  type, y = average_value,  color =model)) +
  geom_point() +
  scale_color_manual(values = oka_col) +
  ylab("average value across taxa") +
  custom_theme(n) +
  labs(color = "model")+
  facet_wrap(~label, scales = "free") +
  theme(axis.title.x = element_blank())
g2


#########################################################################

g3 =  ggplot(mse3, aes(x =  type, y = average_value,  color =model)) +
  geom_point() +
  #scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange")) +
  ylab("average value across taxa") +
  custom_theme(n) +
  theme(axis.title.x = element_blank())
###############################################################
g4 =  ggplot(bias1, aes(x =  type, y = average_value,  color =model)) +
      geom_point() +
    #scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange")) +
      ylab("average value across taxa") +
      custom_theme(n) +
      theme(axis.title.x = element_blank())

g5 =  ggplot(bias2, aes(x =  type, y = average_value,  color =model)) +
  geom_point() +
  #scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange")) +
  ylab("average value across taxa") +
  custom_theme(n) +
  theme(axis.title.x = element_blank())


g6 =  ggplot(bias3, aes(x =  type, y = average_value,  color =model)) +
  geom_point() +
  #scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange")) +
  ylab("average value across taxa") +
  custom_theme(n) +
  theme(axis.title.x = element_blank())


(g1|g2|g3)/(g4|g5|g6) +   plot_layout(guides = "collect") 
  theme(legend.position = "bottom") # Position the legend


g4

View(dd1$error$rr)

 
g2= ggplot(err1, aes(x =  param_name, y = mse, 
                                    group = model, color = model)) +
  geom_point() +
  geom_line() +
  #scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange","true_param" = "blue")) +
  xlab("taxa") +
  ylab("mean squared error") +
  #custom_theme(n) +
  theme(axis.text.x = element_blank()) 

g2
##plots


##### comparing average estimate with true parameter value
n  = 13
g1 = ggplot(pp1, aes(x = param_name, y = average_estimate, group = model, 
                    color= factor(model, levels = unique(model) ))) +
  geom_point() +
  geom_point(aes(x = param_name, y = true_param,color="true_param")) +
  geom_line() +
  #scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange", "true_param" = "blue")) +
  labs(color = "model") +   
  xlab("taxa") +
  ylab("average estimate")  +
  #custom_theme(n) +
  theme_bw() +
  theme(axis.text.x = element_blank()) 

ggsave(paste0(fig_path,"deseq_rr_est.png"), plot = g1, width = 8, height = 6, dpi = 500)
###################################################################
names(dd1$error$rr)
##### Compare their MSE
rr_error    =    error_cal(rr, true_param, model = "with_rr")
deseq_error    =    error_cal(deseq_est, true_param, model = "without_rr")

full_summary_dd_all         =    rbind(rr_error$full_summary_dd,    
                                       deseq_error$full_summary_dd)

single_summary_dd_all       =      rbind(rr_error$single_sum_dd,    
                                       deseq_error$single_sum_dd)
###################################################################
g2= ggplot(full_summary_dd_all, aes(x =  param_name, y = mse, 
                                    group = model, color = model)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange","true_param" = "blue")) +
  xlab("taxa") +
  ylab("mean squared error") +
  custom_theme(n) +
  theme(axis.text.x = element_blank()) 

ggsave(paste0(fig_path,"deseq_rr_mse.png"), plot = g2, width = 8, height = 6, dpi = 500)
###################################################################
g3 = ggplot(err1, aes(x =  param_name, y = bias, group = model, color = model)) +
  geom_point() +
  geom_line() +
 # scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange" )) +
  xlab("taxa") +
  ylab("bias") + 
  #custom_theme(n) +
  theme(axis.text.x = element_blank()) 
g3
ggsave(paste0(fig_path,"deseq_rr_bias.png"), plot = g3, width = 8, height = 6, dpi = 500)
###################################################################
g4=ggplot(single_summary_dd_all, aes(x =  type, y = average_value,
                                     group = type,  color =model)) +
  geom_point() +
  scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange")) +
  ylab("average value across taxa") +
  custom_theme(n) +
  theme(axis.title.x = element_blank())

ggsave(paste0(fig_path,"deseq_rr_avgmse.png"), plot = g4, width = 8, height = 6, dpi = 500)
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
cov2 = coverage_cal(deseq_para_conf)
mean(cov2$coverage)
###################################################################
dd_list  = list()
for(i in 1:length(strg)){
  path    =   paste0(getwd(),strg[[i]])
  dd      =   load_data(path) 
  dd$dd$true_param
  
  conf    =   dd$confint
  ddd     =  do.call(rbind,conf)
  ddd$type  =  rep(i, nrow(ddd))
  dd_list[[i]]  =  dd
  #%>%
  # filter(ddd$model != "deseq_NS" )
}

df       =   data.frame(do.call(rbind, dd_list))
df$type  =   factor(df$type)
df$param_name  =   factor(df$param_name)
str(df)
df$type  =   as.factor(df$type, levels = unique(df$type))

ddf     =  df %>% 
  filter(type == plt_title[[2]])

path    =   paste0(getwd(),strg[[2]])
dd  =   load_data(path) 
conf =  dd$confint
df       =   do.call(rbind, conf)

class(df)

df =  rbind(dd_list[[1]],dd_list[[2]],dd_list[[3]])
ggplot(df, aes(x = param_name, y = average_estimate, color = model, group = model)) +
  geom_point(size = 1, alpha = 0.5) +  
  geom_line()  +
  facet_wrap(~type)

#geom_point(aes(x = param_name, y = true_param, color = "True Parameter"), 
#           size = 1, shape = 17) +   
#scale_color_manual(
#  values = c("blue", "red", "green", "black", "purple", "orange", "cyan"),  
#  name = "Model"
#) + 
# facet_wrap(~type)


str(df)
dd_list
View(dd_list[[2]])


length(unique(dd_list[[1]]$param_name))
df$model  = as.factor(df$model)


df = do.call(rbind, ddd_list)

df$ti =  factor(rep(titles,c(1800,3000,3600)), 
                levels = titles)
ggplot(df, aes(x = param_name, y = average_estimate, color = model, group = model)) +
  geom_point(size = 1, alpha = 0.5) +  
  geom_line() +   
  geom_point(aes(x = param_name, y = true_param, color = "True Parameter"), 
             size = 1, shape = 17) +   
  scale_color_manual(
    values = c("blue", "red", "green", "black", "purple", "orange", "cyan"),  
    name = "Model"
  ) + 
  facet_wrap(~ti)

labs(#title = titles[[i]],#"Comparison of averages of effect size estimates
  #for each model with true effect side",
  x = "Taxa",   
  y = "Average group effect size estimate across simulations",
  color = "Model") +   
  theme_bw(base_size = n) +  
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),
    text = element_text(size = n, family = "Roboto")
  ) 
####################################################################  

