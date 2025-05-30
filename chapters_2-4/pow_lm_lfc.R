source("reproducible/power/Load_Packages.R")
source("reproducible/power/fitting_fun.R")
source("reproducible/power/simulate_fun.R")
source("reproducible/power/utils.R")
source("reproducible/power/read_dataset.R")

fig_path = "reproducible/power/figures/"

#predict power
res <- foreach(k =  1:7,.combine = "rbind",.packages = "tidyr") %do%{
  mod20  =  gam_model_samp[[k]][["20samples"]] 
  mod50  =  gam_model_samp[[k]][["50samples"]] 
  mod100 =  gam_model_samp[[k]][["100samples"]]
  mod150 =  gam_model_samp[[k]][["150samples"]]
  mod200 =  gam_model_samp[[k]][["200samples"]]
  
  pp         =  data.frame(lmean_abund = rep(5,3), abs_lfc = 2:4)
  pp$power0  =  predict(mod20,  newdata = pp,type = "response")
  pp$power1  =  predict(mod50,  newdata = pp,type = "response")
  pp$power2  =  predict(mod100, newdata = pp,type = "response")
  pp$power3  =  predict(mod150, newdata = pp,type = "response")
  pp$power4  =  predict(mod200, newdata = pp,type = "response")
  
  sample_size = rep(c(20, seq(50,200,50)),each=nrow(pp))
  pp_long <- gather(pp, key = "sample", value = "power", -c(lmean_abund,abs_lfc))
  pp_long$sample_size =  sample_size
  pp_long$dataset     =  rep(names(gam_model_samp)[k], nrow(pp_long))
  
  pp_long
  
}

ggplot(res, aes(x = sample_size, y=power, group=abs_lfc, color=as.factor(abs_lfc))) +
  geom_point() +
  geom_line() + 
  scale_color_discrete(name = "logfoldchange") +  
  facet_wrap(~ dataset,scales = "free") +
  labs(x = "sample size", color = "logfoldchange") +
  custom_theme(10) 
  
ggsave(paste0(fig_path,"increase_ss.png"))


#predict power
res2 <- foreach(k =  1:7,.combine = "rbind",.packages = "tidyr") %do%{
  mod500   =  gam_model_otu[[k]][["500otu"]] 
  mod1000  =  gam_model_otu[[k]][["1000otu"]]
  mod1500  =  gam_model_otu[[k]][["1500otu"]]

  pp         =  data.frame(lmean_abund = rep(5,3), abs_lfc = 2:4)
  pp$power0  =  predict(mod500,  newdata = pp,type = "response")
  pp$power1  =  predict(mod1000,  newdata = pp,type = "response")
  pp$power2  =  predict(mod1500, newdata = pp,type = "response")

  
  notu = rep(c(seq(500,1500,500)),each=nrow(pp))
  pp_long <- gather(pp, key = "sample", value = "power", -c(lmean_abund,abs_lfc))
  pp_long$notu =  notu
  pp_long$dataset     =  rep(names(gam_model_otu)[k], nrow(pp_long))
  
  pp_long
  
}

ggplot(res2, aes(x = notu, y = (power), group= abs_lfc, color = as.factor(abs_lfc))) +
  geom_point() +
  geom_line() + 
  scale_color_discrete(name = "logfoldchange") +  # Change the legend name here
  facet_wrap(~ dataset, scales = "free") +
  labs(x = "sample size", color = "logfoldchange") +
  custom_theme(10) 

ggsave("reproducible/power/figures/increase_ss.png")


#' simulate from the 

pp         =  data.frame(lmean_abund  = true_lmean, 
                         abs_lfc = true_lfoldchange)

pp$power1  =  predict(mod50,  newdata = pp, type = "response")
pp$power2  =  predict(mod100, newdata = pp, type = "response")
pp$power3  =  predict(mod150, newdata = pp, type = "response")
pp$power4  =  predict(mod200, newdata = pp, type = "response")

b = (pp$power2 /  pp$power1)#*100

range(b)
plot(b)










