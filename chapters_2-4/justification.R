source("reproducible/power/Load_Packages.R")
source("reproducible/power/utils.R")
source("reproducible/power/read_dataset.R")
source("reproducible/power/simulate_fun.R")
path  =  "reproducible/power/datasets2/"
############################################################
metadata_list_obs   =  readRDS(file = paste0(path,"metadata.rds"))

nsim = 2;   notu   =  1000; nsamp  =  100
disp_vec = c(0.1,0.3,0.9)
############################################################
dd_sig  =  dd_data = list()

for(i in 1:7){
                nam       =   names(metadata_list_obs)[i]
  dispersion_param      =   dispersion_param_list[[i]]
  logmean_param         =   logmean_param_list[[i]]
  logfoldchange_param   =   logfoldchange_param_list[[i]]
  countdata_filt        =   filtered_otu_list[[i]]
  
  for(j in 1:length(disp_vec)){
    disp_scale = disp_vec[j]
    
    ###### Data simulation 
    otu_data   =  countdata_sim_fun(logmean_param, logfoldchange_param, 
                                    dispersion_param,
                                    nsamp_per_group=nsamp,
                                    disp_scale = disp_scale,
                                    ncont = NULL,ntreat = NULL,
                                    notu, nsim, 
                                    seed = 100)
    
    
    countdata_list        =   otu_data$countdata_list
    metadata_list         =   otu_data$metadata_list
  
    pp1 = lapply(countdata_list, function(x){apply(x,1,var)})
    pp2 = lapply(countdata_list, function(x){apply(x,1,mean)})
    
   dd_sig[[j]]    =  data.frame(var  =  log(unlist(pp1)), 
                        mean  =  log(unlist(pp2)),
                       scale  =  rep(disp_scale, length(unlist(pp1))))
  }
   DD    =   do.call(rbind,dd_sig)
   DD$names  =  rep(nam, nrow(DD))
   dd_data[[i]] =  DD
}


df = do.call(rbind,dd_data)
df$scale  = as.factor(df$scale)


dd_list <- Map(function(x, name) {
  dd <- data.frame(mean = log(apply(x, 1, mean)), var = log(apply(x, 1, var)))
  dd$names <- rep(name, nrow(dd))
  dd
}, filtered_otu_list, names(filtered_otu_list))

dd_lis =  do.call(rbind, dd_list)

ddf = df[(df$scale) == "0.3",]

n  = 10
p1 <- ggplot(ddf, aes(x = var, colour = scale)) + 
  geom_density(lwd = 1) +
  geom_density(data = dd_lis, aes(x = var), 
               colour = "black", linetype = "dashed", lwd = 2) +
  theme_bw(base_size = n) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = n, family = "Roboto"),
    axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
    axis.text.y = element_text(family = "Roboto", size = n, color = "black")
  ) +
  facet_wrap(~names, scales = "free")
p1

p2 <- ggplot(ddf, aes(x = mean, colour = scale)) + 
  geom_density(lwd = 0.8) +
  geom_density(data = dd_lis, aes(x = mean), 
               colour = "black", linetype = "dashed", lwd = 2) +
  theme_bw(base_size = n) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = n, family = "Roboto"),
    axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
    axis.text.y = element_text(family = "Roboto", size = n, color = "black")
  ) +
  facet_wrap(~names, scales = "free")

p2
####################################################

### compare variance of simulation and observed
plt1[[i]]  = compare_dataset(countdata_sim_list =  countdata_list[1:5],
                             countdata_filt, method = "var",n=11)

plt2[[i]] = compare_dataset(countdata_sim_list =  countdata_list[1:5],
                            countdata_filt, method = "mean",n=11)



plt1 =  plt2 =  dd_sig = list()
disp_scale = 0.3; nsim = 3
notu   =  1000; nsamp  =  100
for(i in 1:7){
  nam  =   names(metadata_list_obs)[i]
  
  dispersion_param      =   dispersion_param_list[[i]]
  logmean_param         =   logmean_param_list[[i]]
  logfoldchange_param   =   logfoldchange_param_list[[i]]
  countdata_filt        =   filtered_otu_list[[i]]
  
    ###### Data simulation 
    otu_data   =  countdata_sim_fun(logmean_param, logfoldchange_param, 
                                    dispersion_param,
                                    nsamp_per_group=nsamp,
                                    disp_scale = disp_scale,
                                    ncont = NULL,ntreat = NULL,
                                    notu, nsim, 
                                    seed = 100)
    
    
    countdata_list      =   otu_data$countdata_list
    metadata_list         =   otu_data$metadata_list
    true_lmean_list       =   otu_data$logmean_list
    true_lfoldchange_list    =   otu_data$logfoldchange_list
    treat_countdata_list    =   otu_data$treat_countdata_list
    control_countdata_list   =   otu_data$control_countdata_list
    
    pp1 = lapply(countdata_list, function(x){apply(x,1,var)})
    pp2 = lapply(countdata_list, function(x){apply(x,1,mean)})
   
    ### compare variance of simulation and observed
    plt1[[i]]  = compare_dataset(countdata_sim_list =  countdata_list[1:nsim],
                                 countdata_filt, method = "var",n=11)
    
    plt2[[i]] = compare_dataset(countdata_sim_list =  countdata_list[1:nsim],
                                countdata_filt, method = "mean",n=11)
}

(plt1[[1]] | plt1[[2]] | plt1[[3]]) /
  (plt1[[4]] | plt1[[5]] | plt1[[6]]) /
  (plt1[[7]] | plot_spacer() | plot_spacer()) + 
  plot_layout(guides = "collect") 


(plt2[[1]] | plt2[[2]] | plt2[[3]]) /
  (plt2[[4]] | plt2[[5]] | plt2[[6]]) /
  (plt2[[7]] | plot_spacer() | plot_spacer()) + 
  plot_layout(guides = "collect") 


CV = disp = list()
res = foreach(i = 1:7, .combine = "rbind") %do%{
  nam         =   names(logmean_list)[i]
  mean        =   2^logmean_list[[i]]
  dispersion  =   dispersion_list[[i]]
  CV[[i]]     =   data.frame(cv =  coeff_var(mean, dispersion), 
                             data_accession_no  =   rep(nam, length(dispersion)))
  disp[[i]]   =   data.frame(dispersion = dispersion, 
                             data_accession_no  =   rep(nam, length(dispersion)))
}


CV = do.call(rbind,CV)
disp = do.call(rbind,disp)
ggplot(CV, aes(x = data_accession_no, y= cv)) +
  geom_violin(scale = "width") +
  geom_point(alpha= 0.1) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size =10, family = "Roboto"),
    axis.text.x = element_text(family = "Roboto",angle = 90, size =10, color = "black"),
    axis.text.y = element_text(family = "Roboto", size =10, color = "black")
  )+ 
  labs(x="data accession number", y ="coefficient of variation")



ggplot(CV, aes(x = data_accession_no, y= cv)) +
  geom_violin(scale = "width") +
  geom_point(alpha= 0.1) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size =10, family = "Roboto"),
    axis.text.x = element_text(family = "Roboto",angle = 90, size =10, color = "black"),
    axis.text.y = element_text(family = "Roboto", size =10, color = "black")
  )+ 
  labs(x="data accession number", y ="coefficient of variation")


####################################################################################
#Fit lm
combined_list <- mapply(
  function(lmc, lfc) list(lmc = lmc, lfc = lfc),
  lmc = logmean_list,
  lfc = logfoldchange_list,
  SIMPLIFY = FALSE
)

mod=  lapply(combined_list, function(x){lm(x$lfc ~ x$lmc)})

summ_mod = lapply(mod, function(x){summary(x)})
names(summ_mod)

summ_mod[["PRJNA355023"]]
summ_mod[["PRJNA453621"]]
summ_mod[["PRJNA644763"]]


mod = lm(logfoldchange_list[["PRJNA355023"]] ~ logmean_list[["PRJNA355023"]])
summary(mod)

mod = lm(logfoldchange_list[["PRJNA453621"]] ~ logmean_list[["PRJNA453621"]])
summary(mod)

mod = lm(logfoldchange_list[["PRJNA644763"]] ~ logmean_list[["PRJNA644763"]])
summary(mod)

