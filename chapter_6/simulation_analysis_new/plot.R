library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(RhpcBLASctl)
library(Matrix)
library(foreach)
library(huge)
library(here)
source("func.R")
#####################################################
fig_path=   "fig/"
pp0 = pp1  = pp2  =  pp3 = list()
scale  =   1; nn =  18
conf_level <-  0.95; in_sd  = scale*6.7;   alpha  <-  0.05
z_score    <-  qnorm(conf_level + (1 - conf_level)/2)
okabe_ito_palette <- c("orange", "black","violet","blue")
titles  =   c("30 subjects, 100 taxa and 3 time points",
              "50 subjects, 200 taxa and 4 time points")

dodge_width <- 0.5
dodge <- position_dodge(width = dodge_width)

#####################################################
for(param_index   in   1:2){
  source("simulation_analysis_new/initial_param0.R")
  path = paste0("simulation_analysis_new/sim_data/",nsubj,"_",ntaxa,"_",ntime,"/results/")
  path
  #####################################################
  filenames   =    c("coverage_us.rds","coverage_rr.rds",
                     "coverage_nb.rds","coverage_zinb.rds")
  
  dd          =    load_models(path,filenames)
  names(dd)   =    c("US", "RR","NB", "ZINB")
  #####################################################
  if(FALSE){
  test <- bind_rows(
    US = dd$US[[1]],
    RR = dd$RR[[1]],
    NB = dd$NB[[1]]  %>%  select(-sd_err),
    ZINB = dd$ZINB[[1]]  %>%  select(-sd_err),
    .id = "model"
  )
  
  ggplot(test, aes(x = param_name)) +
    geom_pointrange(
      aes(y = est_param, ymin = lwr, ymax = upr, color = model),
      position = dodge,
      linewidth = 0.2,
      size = 0.2
    ) +
    geom_point(
      aes(y = true_param, shape = "True Effect"),
      color = "black", size = 1,
      position = dodge
    ) +
    scale_shape_manual(name = "", values = c("True Effect" = 18)) +
    scale_color_manual(values = okabe_ito_palette) +   
    facet_wrap(~model, scales = "free_y") +
    coord_flip() +
    labs(
      y = "Estimate ± CI",
      x = "Taxa",
      color = "Model",
      title =  titles[param_index]
    ) +
    custom_theme(14) +
    theme(legend.position = "right",
          axis.text.y = element_blank())
  }
  #####################################################
  dd[c("US", "RR")] <- lapply(dd[c("US", "RR")], function(lst) {
    lapply(lst, function(x) {
      x$lwr_old     =   x$lwr
      x$upr_old     =   x$lwr
      x$sd_err      =   (x$est_param - x$lwr) / z_score
      x$lwr         =   x$est_param - z_score*in_sd*x$sd_err
      x$upr         =   x$est_param + z_score*in_sd*x$sd_err
      x
    })
  })
  ########################################################
  sub <- intersect(
    Reduce(intersect, lapply(dd$NB, `[[`, "param_name")),
    Reduce(intersect, lapply(dd$ZINB, `[[`, "param_name"))
  )

  dd <- lapply(dd, function(lst) {
    lapply(lst, function(x) x[x$param_name %in% sub, ])
  })
  ###########################################################################
  d1 = lapply(dd, function(x){
    lapply(x, function(y){
      y$padjust        =   p.adjust(y$pvalue, "BH")
      y$coverage       =   ifelse(y$lwr  < y$true_param & y$true_param < y$upr,1, 0)
      y$pval_reject    =   ifelse(y$padjust  < alpha,1, 0)
      y
    })
  })
 

    
  dd_all <- lapply(names(d1), function(model) {
    data.frame(
      true_param =  d1$US[[1]]$true_param,  
      avg_est    =  rowMeans(sapply(d1[[model]], function(x) x$est_param)),
      width      =  rowMeans(sapply(d1[[model]], function(x) x$width)),
      coverage   =  rowMeans(sapply(d1[[model]], function(x) x$coverage)),
      power      =  rowMeans(sapply(d1[[model]], function(x) x$pval_reject)),
      model      =  model
    )
  })
  ###########################################################################
  coverage_list   =  lapply(names(d1), function(model){
    coverage_matrix = sapply(d1[[model]], function(x) x$coverage)
    p_hat <- mean(coverage_matrix)
    se    <- sqrt(p_hat * (1 - p_hat) / length(coverage_matrix))
    data.frame(lwr   =  p_hat - z_score * se,
               upr   =  p_hat + z_score * se,
              sd_err =  se) })
  
  coverage_dd       =   do.call(rbind,coverage_list)
  coverage_dd$est   =   unlist(lapply(dd_all, function(x){mean(x$coverage)}))
  coverage_dd$model =   names(d1)
  ###########################################################################
  power_list     =  lapply(names(d1), function(model){
    power_matrix =  sapply(d1[[model]], function(x) x$pval_reject)
                           colMeans(power_matrix)[1:200]#[1:min(unlist(
                             #lapply(d1, length)))]
                           })
  
  names(power_list)  =  names(d1)
  power_df <- map_dfr(power_list, ~data.frame(power = .x), .id = "model")
  
  # power_df$taxon  =  paste0("taxon",1:183)
  # lapply(power_list, length)
  
  power_df$sim  = 1:length(power_list$US)
  pp0[[param_index]] <- ggplot(power_df, aes(x = sim, y = power, color = model)) +
    geom_point(size = 2, alpha = 0.5) +  
    geom_line(linewidth =0.7)   +
    labs(x = "simulation index",
         y  = "average statistical power across taxa",
         title =  titles[param_index]) +
    scale_color_manual(values = okabe_ito_palette) +  
    custom_theme(nn)
  
  pow    =  power_df %>% group_by(model) %>% summarise(ave_mean = mean(power))
  print(pow)
  #power_dd$est   =   unlist(lapply(dd_all, function(x){mean(x$power)}))
  #power_dd$model =   names(d1)
  ###########################################################################
  pp1[[param_index]] <- ggplot(coverage_dd, aes(x = model)) +
                       geom_pointrange(aes(y = est, ymin = lwr, ymax = upr),
                              position = dodge, linewidth = 1,size = 0.2, 
                              color = "black")  +
                      geom_hline(yintercept = 0.95, linetype = "dashed",
                                 color = "blue", linewidth =1) +
                      labs(x = " ", y = "average coverage across taxa",
                           title =  titles[param_index]) +
                           custom_theme(nn) 
   ###########################################################################
  ddf <- do.call(rbind, dd_all)
  ddf$type =  if_else(ddf$model %in% c("US", "RR"), "RR and US", "NB and ZINB")
  cover    =  ddf %>% group_by(model) %>% summarise(ave_mean = mean(coverage))
  print(cover)
  ############################################################################
  pp2[[param_index]] <- ggplot(ddf, aes(x = true_param, y = width,  
                  color = model, group = model)) +
    geom_point(size = 2, alpha = 0.5) +  
    geom_line(linewidth =0.7)   +
    scale_color_manual(values = okabe_ito_palette) +   
    labs(x = "true effect", y  = "average confidence width", 
         title =  titles[param_index]) +
    custom_theme(nn) 
  ############################################################################
  pp3[[param_index]] = ggplot(ddf, aes(x =  true_param, y = coverage, color = model)) +
    geom_point(size = 2)+
    geom_line() +
    custom_theme(nn) +
    scale_color_manual(values = okabe_ito_palette) +   
    facet_wrap(~type) +
    labs(x = "true effect size",
         title =  titles[param_index]) 
    
}

p0  =  (pp0[[1]]|pp0[[2]]) + plot_layout(guides = "collect")
p1  =  (pp1[[1]]|pp1[[2]]) + plot_layout(guides = "collect")
p2  =  (pp2[[1]]|pp2[[2]]) + plot_layout(guides = "collect")
p3  =  (pp3[[1]]|pp3[[2]]) + plot_layout(guides = "collect")

width =  12; height =  5; dpi = 300
ggsave("fig/power_long.png", plot = p0, 
       width = width, 
       height = height, 
       dpi = dpi)

ggsave("fig/ave_cover_long.png", plot = p1, 
       width = width, 
       height = height, 
       dpi = dpi)


ggsave("fig/ave_width_long.png", plot = p2, 
       width = width, 
       height = height, 
       dpi = dpi)

# pp2[[param_index]] <- ggplot(cover , aes(x =  model, y =  ave_mean)) +
#   geom_point(size = 2)+
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "blue", linewidth =1) +
#   custom_theme(nn) +
#   labs(x = " ",  y = "average coverage across taxa", title =  titles[param_index]) +
#   scale_color_manual(values = okabe_ito_palette) 

# ggplot(ddf, aes(x = param_name, y = CI_width,  
#                 color = model, group = model)) +
#   geom_point(size = 2, alpha = 0.5) +  
#   geom_line(linewidth =0.7)   +
#   custom_theme(14) +
#   theme(axis.text.x = element_blank()) +
#   labs(x = "taxa", y  = "confidence width")

##################################################################
pp11  =  (pp1[[1]]/pp1[[2]]) +   plot_layout(guides = "collect")
pp22  =  (pp2[[1]]/pp2[[2]]) +   plot_layout(guides = "collect")
pp33  =  (pp3[[1]]/pp3[[2]]) +   plot_layout(guides = "collect")


width =  8; height =  9; dpi = 300
ggsave("fig/ave_confwidth.png", plot = pp11, 
       width = width, 
       height = height, 
       dpi = dpi)

ggsave("fig/ave_covv.png", plot = pp22, 
       width = width, 
       height = height, 
       dpi = dpi)

ggsave("fig/covv_per_taxa.png", plot = pp33, 
       width = 9, 
       height = height, 
       dpi = dpi)
##################################################################

coverage_dd <- lapply(names(d1), function(model) {
  data.frame(
    true_param = d1$US[[1]]$true_param,  
    coverage   = rowMeans(sapply(d1[[model]], function(x) x$coverage)),
    model      = model
  )
})

########################################################
path0  =  "simulation_analysis_new/sim_data/"
####Bias and RMSE calculation
strg    =   c("30_100_3/","50_200_4/")
########################################################
n = 24; p1  =   p2   =  p3 =  p4  = p5  = p6 = p7 = list()
BIAS = VARIAN  =  list()
sz =2; wd  = 0.8; lw = 2
titles  =   c("30 subjects, 100 taxa, 3 time points",
              "50 subjects, 200 taxa, 4 time points")
plt1 = list()
nn = 18
plt  = list()
for(i in 1:length(strg)){
  ########################################################
  path  =   paste0(path0,strg[[i]],"results/") 
  dd    =   load_data(path) 
  
  conf =  dd$confint
  ddd  =  do.call(rbind,conf)
  ddd$type <- if_else(ddd$model %in% c("RR", "US"), "RR and US", "NB and ZINB")
  
  dodge_width <- 0.5
  dodge <- position_dodge(width = dodge_width)
  
  # names(d1$RR[[1]])
  # ddd  =  cbind(d1$US[[1]], d1$RR[[1]], d1$NB[[1]], d1$ZINB[[1]])
  # View(test)
plt1[[i]] = ggplot(ddd, aes(x = param_name)) +
    geom_pointrange(
      aes(y = average_estimate, ymin = lwr, ymax = upr, color = model),
      position = dodge,
      linewidth = 0.2,
      size = 0.2
    ) +
    geom_point(
      aes(y = true_param, shape = "True Effect"),
      color = "black", size = 1,
      position = dodge
    ) +
    scale_shape_manual(name = "", values = c("True Effect" = 18)) +
    scale_color_manual(values = c("RR" = "#1b9e77", "US" = "#d95f02", "NB" = "#7570b3", "ZINB" = "#e7298a")) +
    facet_wrap(~model, scales = "free_y") +
    coord_flip() +
    xlim(levels(ddd$param_name)) +
    labs(
      #y = "Estimate ± CI",
      y = "Average effect size estimate (with confidence interval)",
      x = "Taxa",
      color = "Model",
      title =  titles[i]
    ) +
    custom_theme(14) +
    theme(legend.position = "right",
          axis.text.y = element_blank())
}

plt1[[2]] 
width =  8; height =  6; dpi = 300
ggsave("fig/confint1.png", plot = plt1[[1]], 
       width = width, 
       height = height, 
       dpi = dpi)

ggsave("fig/confint2.png", plot = plt1[[2]], 
       width = width, 
       height = height, 
       dpi = dpi)
    
  
  
plt[[i]] <-  ggplot(ddd, aes(x = param_name, y = CI_width,  
                  color = model, group = model)) +
                  geom_point(size = 2, alpha = 0.5) +  
                  geom_line(linewidth =0.7)   +
                  custom_theme(nn) +
                  theme(axis.text.x = element_blank()) +
                  labs(title  = titles[i] , x = "taxa", 
                       y  = "confidence width") 

    err  =   lapply(error$RR$full_summary_dd)
View(error)

dd_est  =  do.call(rbind,dd$est_dd)

ggplot(ddd, aes(x = param_name, y = average_estimate, 
                color = model, group = model)) +
  geom_point(size = 1, alpha = 0.5) +  
  geom_line()    +
geom_point(aes(x = param_name, y = true_param, color = "true group effect"), 
           size = 1.5, shape = 17)
# ylim(-5, 5)  +
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


}

pp = (plt[[1]]|plt[[2]]) + plot_layout(guides = "collect")
width =  8; height =  6; dpi = 300
ggsave("fig/confint1.png", plot = plt1[[1]], 
       width = width, 
       height = height, 
       dpi = dpi)

ggsave("fig/confint2.png", plot = plt1[[2]], 
       width = width, 
       height = height, 
       dpi = dpi)

ggsave("fig/para_conf_wdth.png", plot = pp, 
       width = width, 
       height = height, 
       dpi = dpi)

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

df$coverage =   ifelse(df$lwr  < df$true_param & df$true_param < df$upr,1, 0)
y$p.adj     =   p.adjust(y$pvalue, "BH")

x$padjust  =  p.adjust(x$pvalue, method = "BH")
x$pow      =  ifelse(x$padjust < sig, 1, 0)