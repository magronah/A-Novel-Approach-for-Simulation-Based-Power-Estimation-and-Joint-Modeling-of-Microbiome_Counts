library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(here)
source("func.R")
fig_path=   "fig/"
#####################################################
path0  =  "simulation_analysis_new/sim_data/"
strg    =   c("30_100_3/","50_200_4/")
titles  =   c("30 subjects, 100 taxa and 3 time points",
              "50 subjects, 200 taxa and 4 time points")
########################################################
n = 18; p1  =   p2   =  p3 =  p4  = p5  = p6 = p7 = list()
# sz =2; wd  = 0.8; lw = 2
titles  =   c("30 subjects, 100 taxa, 3 time points",
              "50 subjects, 200 taxa, 4 time points")
oka_col = c("#556B2F", "#E23D28","#0000FF","#E69F00","#000000")  

for(i in 1:length(strg)){

  path  =   paste0(path0,strg[[i]],"results/") 
  dd    =   load_data(path) 

  conf =  dd$confint
  ddd  =  do.call(rbind,conf)
  p1[[i]] <- make_plot(ddd, "true_param", "CI_width", 
                             "Confidence width", titles[[i]], 
                             show_true = FALSE, 
                             use_custom_theme = TRUE,
                           x_label  = "True effect")
                        
  p2[[i]] <- make_plot(ddd, "param_name", "average_estimate", 
                                   "Average effect estimate", titles[[i]],
                                   show_true = TRUE,
                       x_label  = "True effect")
  p2[[i]]  =    p2[[i]] + labs(x = "taxa") +    theme(axis.text.x = element_blank())

  
  
  err       =  combine_model_component(dd,"full_summary_dd")
  err$rmse  =  sqrt(err$mse)
  
  # p3[[i]] <- make_plot(err, "true_param", "bias", "Bias", titles[[i]], 
                       # x_label  = "True effect")
  
  #p3[[i]] <- p3[[i]]  +  geom_hline(yintercept = 0, linetype = "dashed", 
  #                                  linewidth  = 1,
  #                                  color = "black")
  err$type  =  ifelse(err$model %in% c("RR","US"), "RR and US", "NB and ZINB")
  err$sd_err =  sqrt(err$variance)
  
  p3[[i]]  = ggplot(err, aes(true_param, bias, color = model))+
    geom_point() +
    geom_line()  + 
    labs(x = "True effect", y  = "Bias", title = titles[i]) +
    geom_hline(yintercept = 0, linetype = "dashed", 
               linewidth  = 1,
               color = "black") +
    scale_color_manual(values = oka_col) +
    facet_wrap(~type) +
    custom_theme(20)


  p4[[i]] <- make_plot(err, "true_param", "rmse", "Root mean squared error", 
                       titles[[i]], x_label  = "True effect")
  p5[[i]] <- make_plot(err, "true_param", "sd_err", "Standard deviation of errors",
                       titles[[i]], x_label  = "True effect")
  
  #make_plot(err, "true_param", "variance", "Variance of errors",
  #          titles[[i]], x_label  = "True effect")
  
 av_bias <- err %>% 
    group_by(model)%>%
    summarise(averg_bias = mean(bias))
 p6[[i]] = ggplot(av_bias, aes(model,averg_bias)) +
    geom_point(size =3) +
    labs(x = " ", y  = "Average bias across taxa",title = titles[i]) +
    custom_theme(n)
}
#######################################################################
pp1  = (p1[[1]]|p1[[2]]) + plot_layout(guides = "collect")
pp2  = (p2[[1]]|p2[[2]]) + plot_layout(guides = "collect")
pp3  = (p3[[1]]/p3[[2]]) + plot_layout(guides = "collect")
pp4  = (p4[[1]]|p4[[2]]) + plot_layout(guides = "collect")
pp5  = (p5[[1]]|p5[[2]]) + plot_layout(guides = "collect")
pp6  = (p6[[1]]|p6[[2]]) + plot_layout(guides = "collect")
#######################################################################
width =  15; height =  6; dpi = 300
ggsave("fig/para_conf_wdth.png", plot = pp1, 
       width = width, 
       height = height, 
       dpi = dpi)



ggsave("fig/trend_long.png", plot = pp2, 
       width = width, 
       height = height, 
       dpi = dpi)

ggsave("fig/bias_long2.png", plot = pp3, 
       width = 15, 
       height = 15, 
       dpi = dpi)

ggsave("fig/rmse_long.png", plot = pp4, 
       width = width, 
       height = height, 
       dpi = dpi)

ggsave("fig/sdr_long.png", plot = pp5, 
       width = width, 
       height = height, 
       dpi = dpi)

ggsave("fig/ave_bias_long.png", plot = pp6, 
       width = width, 
       height = height, 
       dpi = dpi)
#######################################################################
















for(i in 1:length(strg)){
########################################################
  path  =   paste0(path0,strg[[i]],"results/") 
  dd    =   load_data(path) 

  conf =  dd$confint
  ddd  =  do.call(rbind,conf)

ggplot(ddd, aes(x = param_name, y = CI_width,  
                  color = model, group = model)) +
    geom_point(size = 2, alpha = 0.5) +  
    geom_line(linewidth =0.7)   +
    custom_theme(14) +
    theme(axis.text.x = element_blank()) +
    labs(x = "taxa", y  = "confidence width")
  
  #######################################################
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
  p4[[i]] = ggplot(var, aes(model, average_value)) +
    geom_point(size = 3*sz) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), size=sz, width = wd) +
    geom_hline(yintercept = (var["RRzi",1]), linewidth = lw, linetype = "dashed", color = "red") +
    custom_theme(n) +
    labs(title = titles[[i]],
         x = " ",
         y = "Average variance of error across taxa"
         #"Comparison of average bias across taxa",x = " ",y = "Average Bias"
    )
  
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

var_p <- (p4[[1]]|p4[[2]]|p4[[3]]) + plot_layout(guides = "collect")  
ggsave("fig/var_plt.png", plot = var_p, 
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


















p1[[i]] <-  ggplot(ddd, aes(x = param_name, y = CI_width,  
                            color = model, group = model)) +
  geom_point(size = 2, alpha = 0.5) +  
  geom_line(linewidth =0.7)   +
  custom_theme(nn) +
  scale_color_manual(values = oka_col) +  
  theme(axis.text.x = element_blank()) +
  labs(title  = titles[i] , x = "taxa", 
       y  = "confidence width") 

p2[[i]] <- ggplot(ddd, aes(x = param_name, y = average_estimate, 
                           color = model, group = model)) +
  geom_point(size = 1, alpha = 0.5) +  
  geom_line(linewidth = 0.7)    +
  geom_point(aes(x = param_name, y = true_param, color = "true effect"), 
             size = 3, shape = 17) + # ylim(-5, 5)  +
  scale_color_manual(values = oka_col) +  
  labs(title = titles[[i]],
       x = "Taxa",   
       y = "Average  effect estimate",
       color = "Model") +  
  theme_bw(base_size = n) +
  theme(
    plot.title = element_text(hjust = 0.5, size = n),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),  
    # axis.ticks.x = element_blank(),
    text = element_text(size = n, family = "Roboto")
  )  
