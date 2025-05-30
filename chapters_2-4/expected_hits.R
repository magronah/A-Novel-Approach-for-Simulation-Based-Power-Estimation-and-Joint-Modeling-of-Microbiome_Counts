source("reproducible/power/Load_Packages.R")
source("reproducible/power/simulate_fun.R")
source("reproducible/power/utils.R")
source("reproducible/power/read_dataset.R")

fig_path = "reproducible/power/figures/"
#################Expected number of hits##############
## show this with the confidence intervals, divide by sqrt(n)
sample_vec = seq(10,100,10)
 
l = length(contour_data); ll = length(contour_data[[1]])

#i =j=1 
grid_len =50
ncores  =  10
cl      =   makeCluster(ncores)
registerDoParallel(cl)

result_list <- foreach (i = 1:l, .combine = "rbind") %:%
  foreach (j = 1:ll, .combine = "rbind") %do% {
    dd = contour_data[[i]][[j]]
    nam = names(contour_data)[i]
    fit_2d = readRDS(paste0(path,"gam_model/50_samples/",nam,".rds"))
    
    pp      =   with(dd,
                     expand.grid(lmean_abund = seq(min(lmean_abund),
                                                   max(lmean_abund),
                                                   length  = grid_len),
                                 abs_lfc   =  seq(min(abs_lfc),
                                                  max(abs_lfc),
                                                  length  =  grid_len)))
    # #predict power
    power <- predict(fit_2d, newdata = pp,type = "response")
    df = data.frame(expected = nrow(pp) * mean(power))
    #print(mean(power))
    #df = data.frame(expected = nrow(dd) * mean(dd$pval_reject))
    df$sample_size = sample_vec[j]
    df$dataset   =  names(contour_data)[i]
    df
  }

stopCluster(cl)
unregister_dopar() 

#View(contour_data[[1]][[1]])
ggplot(result_list, aes(x= sample_size, y= (expected/100),
                        group=(dataset),colour=dataset)) +
  geom_point() +
  geom_line() +
  labs(x="sample size", y = "expected number of significant taxa",
       color ="data accession no.")+
  custom_theme(10)

ggsave(paste0(fig_path,"expected_hit.png"))
#######################################################################
nsim = 1; notu=100; disp_scale = 0.3
nsamp_vec1 = sample_vec

i=4
  nam       =   names(dispersion_param_list)[i] 
  dispersion_param      =   dispersion_param_list[[i]]
  logmean_param         =   logmean_param_list[[i]]
  logfoldchange_param   =   logfoldchange_param_list[[i]]

  dd = foreach(k = 1:length(nsamp_vec1)) %do% {
    source("reproducible/power/simulate_fun.R")
    nsamp     =   nsamp_vec1[[k]]
    countdata_sim_fun(logmean_param, logfoldchange_param, 
                      dispersion_param,
                      nsamp_per_group=nsamp,
                      disp_scale = disp_scale,
                      ncont = NULL,ntreat = NULL,
                      notu, nsim, 
                      maxlfc_iter = 10000,
                      seed = 100)
    
  }
  names(dd) = paste0("sample", nsamp_vec1)  


grid_len =30
ncores  =  10
cl      =   makeCluster(ncores)
registerDoParallel(cl)

dd=data.frame(lmean_abund      =   unlist(read_data(dd, "logmean_list")),
abs_lfc          =   unlist(read_data(dd, "logfoldchange_list")) )
dd$sample_size = rep(nsamp_vec1, each=notu)

View(lmean_abund)
i=4
#sample_s = seq(5,50,5)
result_list <- foreach (i = 1:l, .combine = "rbind",.packages = "scam") %:%
kk= foreach (j = 1:ll, .combine = "rbind") %do% {
    #dd = contour_data[[i]][[j]]
    nam = names(contour_data)[i]
    fit_2d = readRDS(paste0(path,"sample_size_cal/gam_fit/",nam,".rds"))
plot(fit_2d)
    pp      =   with(dd,
                     expand.grid(lmean_abund = seq(min(lmean_abund),
                                                   max(lmean_abund),
                                                   length  = grid_len),
                                 abs_lfc   =  seq(min(abs_lfc),
                                                  max(abs_lfc),
                                                  length  =  grid_len)))
    
    pp$sample_size = rep(sample_vec[i], nrow(pp))
    # #predict power
    power <- predict(fit_2d, newdata = dd,type = "response")
    dd$power = power
    df = data.frame(expected = nrow(dd) * mean(power))
    print(list(mean(power), pp$sample_size[[j]]))
    #df = data.frame(expected = nrow(dd) * mean(dd$pval_reject))
    df$sample_size = sample_vec[j]
   # df$dataset   =  names(contour_data)[i]
    df
  }

ggplot(result_list, aes(x= sample_size, y= (expected/100),
                        group=(dataset),colour=dataset)) +
  geom_point() +
  geom_line() +
  labs(x="sample size", y = "expected number of significant taxa",
       color ="data accession no.")+
  custom_theme(10)

######## average power with squat number of hits   ##############
plt = list()
for(i in 1:7){
nam  =   names(dispersion_param_list)[i]

dispersion_param      =   dispersion_param_list[[i]]
logmean_param         =   logmean_param_list[[i]]
logfoldchange_param   =   logfoldchange_param_list[[i]]

###### Data simulation 
nsamp = 100; nsim = 100; notu=1000; disp_scale = 0.3

otu_data   =  countdata_sim_fun(logmean_param, logfoldchange_param, 
                                dispersion_param,
                                nsamp_per_group=nsamp,
                                disp_scale = disp_scale,
                                ncont = NULL,ntreat = NULL,
                                notu, nsim, 
                                seed = 100)

pow_vec20 = pow_vec50 = pow_vec100 = pow_vec150 =  pow_vec200 = list()

mod20  =  gam_model_samp[[i]][["20samples"]]
mod50  =  gam_model_samp[[i]][["50samples"]] 
mod100 =  gam_model_samp[[i]][["100samples"]]
mod150 =  gam_model_samp[[i]][["150samples"]]
mod200 =  gam_model_samp[[i]][["200samples"]]

pp <- foreach(j = 1:nsim, .combine = "rbind") %do% {
  true_lmean         =   otu_data$logmean_list[[j]]
  true_lfoldchange   =   otu_data$logfoldchange_list[[j]]
  dd = data.frame(lmean_abund = true_lmean, abs_lfc =true_lfoldchange)

  pow0  =  (predict(mod20,  newdata = dd,type = "response"))
  pow1  =  (predict(mod50,  newdata = dd,type = "response"))
  pow2  =  (predict(mod100, newdata = dd,type = "response"))
  pow3  =  (predict(mod150, newdata = dd,type = "response"))
  pow4  =  (predict(mod200, newdata = dd,type = "response"))
  
  pow_vec20[[j]]  =  pow0  
  pow_vec50[[j]]  =  pow1
  pow_vec100[[j]] =  pow2
  pow_vec150[[j]] =  pow3
  pow_vec200[[j]] =  pow4
}

l = list(unlist(pow_vec20),  unlist(pow_vec50), unlist(pow_vec100),
         unlist(pow_vec150), unlist(pow_vec200))

dd   =    data.frame(do.call(rbind,l))
dd   =    data.frame(min  =  apply(dd,1,min), max  =  apply(dd,1,max))

quant = do.call(rbind, lapply(l,function(x) quantile(x, probs = c(0.50, 0.60,0.70, 0.80))))
dd    =  cbind(dd,quant)

sample_size     =  c(20,seq(50,200, 50))  
dd$sample_size  =  factor(sample_size, levels=unique(sample_size))

mean_list <- list(sapply(pow_vec20, mean), sapply(pow_vec50, mean),
                  sapply(pow_vec100, mean), sapply(pow_vec150, mean),
                  sapply(pow_vec200, mean))

dd$sd     =   unlist(lapply(mean_list, sd))
dd$mean   =   unlist(lapply(mean_list, mean))

names(dd)[names(dd) == "50%"]   =  "quantile50"
names(dd)[names(dd) == "60%"]   =  "quantile60"
names(dd)[names(dd) == "70%"]   =  "quantile70"
names(dd)[names(dd) == "80%"]   =  "quantile80"

plt[[i]] =  ggplot(dd, aes(x = sample_size, y = mean,group=sample_size)) +
          geom_point()  +  
          geom_line(aes(y = quantile50), color = "blue",group = 1) +
          geom_line(aes(y = quantile60), color = "black",group = 1) +
          geom_line(aes(y = quantile70), color = "green",group = 1) +
          geom_line(aes(y = quantile80), color = "red",group = 1) 
  

}


p1  =  plt[[1]] +  
  geom_line(aes(y = quantile50, color = "50th Percentile", group = 1)) +
  geom_line(aes(y = quantile60, color = "60th Percentile", group = 1)) +
  geom_line(aes(y = quantile70, color = "70th Percentile", group = 1)) +
  geom_line(aes(y = quantile80, color = "80th Percentile", group = 1)) +
  scale_color_manual(name = "Quantiles", 
                     values = c("50th Percentile" = "blue", 
                                "60th Percentile" = "black", 
                                "70th Percentile" = "green", 
                                "80th Percentile" = "red")) +
  custom_theme(10) + labs(x="sample size", y="overall power") 
p2  =  plt[[2]] +  
  geom_line(aes(y = quantile50, color = "50th Percentile", group = 1)) +
  geom_line(aes(y = quantile60, color = "60th Percentile", group = 1)) +
  geom_line(aes(y = quantile70, color = "70th Percentile", group = 1)) +
  geom_line(aes(y = quantile80, color = "80th Percentile", group = 1)) +
  scale_color_manual(name = "Quantiles", 
                     values = c("50th Percentile" = "blue", 
                                "60th Percentile" = "black", 
                                "70th Percentile" = "green", 
                                "80th Percentile" = "red")) +
  custom_theme(10) + labs(x="sample size", y="overall power") 

p5  =  plt[[5]] +  
  geom_line(aes(y = quantile50, color = "50th Percentile", group = 1)) +
  geom_line(aes(y = quantile60, color = "60th Percentile", group = 1)) +
  geom_line(aes(y = quantile70, color = "70th Percentile", group = 1)) +
  geom_line(aes(y = quantile80, color = "80th Percentile", group = 1)) +
  scale_color_manual(name = "Quantiles", 
                     values = c("50th Percentile" = "blue", 
                                "60th Percentile" = "black", 
                                "70th Percentile" = "green", 
                                "80th Percentile" = "red")) +
  custom_theme(10) + labs(x="sample size", y="overall power") 

p6  =  plt[[6]] +  
  geom_line(aes(y = quantile50, color = "50th Percentile", group = 1)) +
  geom_line(aes(y = quantile60, color = "60th Percentile", group = 1)) +
  geom_line(aes(y = quantile70, color = "70th Percentile", group = 1)) +
  geom_line(aes(y = quantile80, color = "80th Percentile", group = 1)) +
  scale_color_manual(name = "Quantiles", 
                     values = c("50th Percentile" = "blue", 
                                "60th Percentile" = "black", 
                                "70th Percentile" = "green", 
                                "80th Percentile" = "red")) +
  custom_theme(10) + labs(x="sample size", y="overall power")

p
p2  =  plt[[2]] + custom_theme(10) + labs(x="sample size", y="overall power") +p
p7  =  plt[[6]] + custom_theme(10)+ labs(x="sample size", y="overall power")+p


(p1|p2)/(p5|p6) + plot_layout(guides = "collect")  
ggsave(paste0(fig_path,"overall_pow.png"))





(plt[[1]]|plt[[2]]|plt[[3]]|plt[[4]])/(plt[[5]]|plt[[6]]|plt[[7]])  +  plot_layout(guides = "collect")  

plt[[1]]|plt[[3]]
  geom_line(aes(y = max), color = "green",group = 1) 

  geom_line(aes(y  = quantile25,group = 1)) +
  geom_line(aes(y  = quantile50,group = 1)) 
  geom_line(aes(y  = min,group = 1)) +
  geom_line(aes(y  = max,group = 1)) 
  

  scale_y_continuous(limits = c(NA,0.1),oob = scales::squish) +
  geom_linerange(aes(ymin = mean - 1.96*sd, ymax = mean + 1.96*sd)) +
  geom_linerange(aes(x = sample_size+5,ymin = min, ymax = max), position = position_dodge(width = 3)) 

  
  
  ggplot(dd, aes(x = sample_size, y = mean)) +
    geom_point()  +
    geom_linerange(aes(ymin = mean - 1.96*sd, ymax = mean + 1.96*sd)) 
    
  

library(ggplot2) 
pd <- position_dodge(width=0.4) 
g0 + geom_linerange(aes(ymin=pred-cmult*SE,ymax=pred+cmult*SE), position=pd)


  geom_line(aes(y = quantile_25), color = "blue") +
  geom_line(aes(y = quantile_50), color = "green") 

  geom_line(aes(group = variable)) +  # Add lines for quantiles
  labs(x = "Sample Size", y = "Value") +  # Set axis labels
  scale_color_discrete(name = "Statistics") +  # Set legend title

  


pdd <- pp %>% group_by(sample_size) %>% 
  summarise_at(vars(pow), list(mean = mean,sd = sd))




 
ggplot(pdd, aes(x = sample_size, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - 1.96*sd, ymax = mean + 1.96*sd), width = 0.2) 

  geom_line(aes(y = quantile_25), color = "blue") +
  geom_line(aes(y = quantile_50), color = "green") +
  geom_line(aes(y = quantile_75), color = "red") 
  labs(x = "Sample Size", y = "Mean") +
  ggtitle("Mean with Quantiles for Each Sample Size")


pd <- pp %>%
  group_by(sample_size) %>%
  summarise_at(vars(pow), list(range = ~paste(round(range(.),3), collapse = " - "),
                               mean = mean,
                               sd = sd,
                               quantile_25 = ~quantile(., 0.25),
                               quantile_50 = ~quantile(., 0.5),
                               quantile_75 = ~quantile(., 0.75)))



sample_vec =  c(20, seq(50,200,50))
cbind(pd,t(column_ranges))




mod20  =  gam_model_samp[[k]][["20samples"]] 
mod50  =  gam_model_samp[[k]][["50samples"]] 
mod100 =  gam_model_samp[[k]][["100samples"]]
mod150 =  gam_model_samp[[k]][["150samples"]]
mod200 =  gam_model_samp[[k]][["200samples"]]

pp1 <- foreach(j = 1:nsim, .combine = "rbind") %do% {
  true_lmean         =   otu_data$logmean_list[[j]]
  true_lfoldchange   =   otu_data$logfoldchange_list[[j]]
  data.frame(lmean_abund = true_lmean, abs_lfc =true_lfoldchange)
}

pp1$power0  =  predict(mod20,  newdata = pp1,type = "response")
pp1$power1  =  predict(mod50,  newdata = pp1,type = "response")
pp1$power2  =  predict(mod100, newdata = pp1,type = "response")
pp1$power3  =  predict(mod150, newdata = pp1,type = "response")
pp1$power4  =  predict(mod200, newdata = pp1,type = "response")


dr = data.frame(min = unlist(sapply(pp_filt,min)),
           max = unlist(sapply(pp_filt,max)))

p=t(apply(pp_filt,2,  function(x)quantile(x, probs = c(0.25, 0.5, 0.75))))

ppr =cbind(dr,p)
ppr$sample_size = 
  
  
pp_filt <- pp1 %>% select(-c(lmean_abund,abs_lfc))
pr=pp1 %>%  summarise_at(vars(power0:power4), list(range = ~paste(round(range(.),3), collapse = " - "),
                                           mean = mean,
                                           sd = sd,
                                           quantile_25 = ~quantile(., 0.25),
                                           quantile_50 = ~quantile(., 0.5),
                                           quantile_75 = ~quantile(., 0.75)))  %>%
     data.frame()

View(pr)
data.frame(summary(pp_filt))

column_ranges <- sapply(pp_filt, range)
# Calculate column means for the remaining columns
average_pow <- colMeans(pp_filt)



# The number of expected number of taxa
#' For each sample size, I predict the power using the simulated data 
#' I stack all the power together and compute the mean power and multiply by the 
#' length of otu to get the number of expected number of taxa
#' 
#' 
#' I am also computing the mean power from the data directly
#' hem using the gam
#' 
# Power decreases with increasing otu

#the plot to say that the average power is below the range of true power