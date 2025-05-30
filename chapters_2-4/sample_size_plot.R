source("reproducible/power/Load_Packages.R")
source("reproducible/power/utils.R")
source("reproducible/power/read_dataset.R")
source("reproducible/power/simulate_fun.R")
path = "reproducible/power/datasets2/sample_size_cal/"
fig  = paste0(path,"fig/")
library(tidytext)
library(metR)
################################################################################
ss  = c(1,2,7)
dd_lfc  = foreach(i =ss, .combine = "rbind")%do%{
  nam   =    names(logfoldchange_list)[i]
  dd    =    data.frame(logfoldchange = logfoldchange_list[[i]])
  dd$names = rep(nam, nrow(dd))
  dd
}

nm = 20
p1 = ggplot(dd_lfc, aes(names,logfoldchange,group=names)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +  # Add violin plot with some transparency
  geom_boxplot(width = 0.1, color = "darkblue") + # Add box plot
  # geom_hline(yintercept = -1,       # Add horizontal line at y = 0
  #            color = "red",        # Set line color to red
  #            linetype = "dashed",  # Set line type to dashed
  #            linewidth = 1) +
  # geom_hline(yintercept = 1,       # Add horizontal line at y = 0
  #            color = "red",        # Set line color to red
  #            linetype = "dashed",  # Set line type to dashed
  #            linewidth = 1) +
ylab(TeX("$\\log_2$(fold change)")) +
  labs(x = " ") +
  custom_theme(nm)

ggsave(filename = paste0(fig,"lfc_range.png"), plot = p1, 
       width = 15, height = 5, dpi = 300)
################################################################################
dd_lmb  = foreach(i =ss, .combine = "rbind")%do%{
  nam   =    names(logmean_list)[i]
  dd    =    data.frame(logmean = logmean_list[[i]])
  dd$names = rep(nam, nrow(dd))
  dd
}

p2 = ggplot(dd_lmb, aes(names,logmean,group=names)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +   
  geom_boxplot(width = 0.1, color = "darkblue") +  
  ylab(TeX("$\\log_2$(mean count)")) +
  labs(x = " ") +
  custom_theme(nm) 


ggsave(filename = paste0(fig,"lmb_range.png"), 
       plot = p2, width = 15,height = 5, dpi = 300)

py = p1|p2
ggsave(filename = paste0(fig,"lmb_lfcrange.png"), 
       plot = py, width = 15,height = 5, dpi = 300)
#################################################################################
gam_mod_list  =  foreach(i = 1:7) %do% {
  nam      =   names(logmean_list)[i]
  mod      =   readRDS(paste0(path,"gam_fit/",nam, ".rds"))
  mod
}
names(gam_mod_list) = names(filtered_otu_list)
######################################################################
### estimate sample size for each power with changing effect size
lfc    =  c(0,1,2,3);  lmab = 5
target = seq(0.1, 0.99, 0.05); len = length(target)

dd = foreach(k = 1:7, .combine = "rbind") %:%
   foreach(i = 1:length(lfc) ,.combine = "rbind") %do%{
    
    nam         =   names(filtered_otu_list)[k]
    gam_fit_ss  =   gam_mod_list[[k]]
    
    abs_lfc = rep(lfc[i],len)
    lmb     = rep(lmab,len)
    
    pp   =    inverse_fun(target=target,lmb=lmb,abs_lfc=abs_lfc, 
                          model=gam_fit_ss,xmin=10, xmax=100)
    
    d    =    data.frame(sample_size = pp, lmb=lmb, 
                         lfc = abs_lfc, power=target)
    
  d$lfc  =    factor(d$lfc, levels = unique(d$lfc))
    d$dataset = rep(nam, nrow(d))
    d
  }

oka_col = c(
  "#556B2F", 
  "#E23D28", 
  "#0000FF",
  "#E69F00" # Orange
)

sub  =  c("PRJNA168470", "PRJNA355023","PRJNA687773")# "PRJNA453621")
dd_sub = dd[dd$dataset %in% sub,]

pp1 = ggplot(dd_sub, aes(sample_size,power, color=lfc)) +
  geom_point(size=3) +
  geom_line(lwd=1.5) + 
  ylim(0,1) +  
  scale_color_manual(name = TeX("|$\\log_2$(fold change)|"),
                     values = oka_col) +
  facet_wrap(~dataset, scale="free") + 
  guides(color = guide_legend(reverse = TRUE), 
         linetype = guide_legend(reverse = TRUE)) +
  labs(x = "sample size per group") +
  custom_theme(18) 

ggsave(paste0(fig,"lfc.png"), plot = pp1, width = 15, height = 5, dpi = 300)
################################################################################
### sample size estimation for different lmb
lfc    =  2;  lmb_vec = c(-1,2,7)
target = seq(0.1, 0.99, 0.05); len = length(target)

dd = foreach(k = 1:7, .combine = "rbind") %:%
  foreach(i = 1:length(lmb_vec) ,.combine = "rbind") %do%{
    
    nam         =   names(filtered_otu_list)[k]
    gam_fit_ss  =   gam_mod_list[[k]]
    
    abs_lfc = rep(lfc,len)
    lmb     = rep(lmb_vec[i],len)
    
    pp   =    inverse_fun(target=target,lmb=lmb,abs_lfc=abs_lfc, 
                          model=gam_fit_ss,xmin=10, xmax=100)
    
    d    =    data.frame(sample_size = pp, lmb=lmb, 
                         lfc = abs_lfc, power=target)
    
    d$lmb  =    factor(d$lmb, levels = unique(d$lmb))
    d$dataset = rep(nam, nrow(d))
    d
  }

#sub=c("PRJNA589343", "PRJNA644763","PRJNA687773")
dd_sub = dd[dd$dataset %in% sub,]

pp2 = ggplot(dd_sub, aes(sample_size,power, color=lmb)) +
  geom_point(size=3) +
  geom_line(lwd=1.5) + 
  ylim(0,1) +  
  scale_color_manual(name = TeX("$\\log_2$(mean count)"),values = oka_col) +
  facet_wrap(~dataset, scale="free") + 
  guides(color = guide_legend(reverse = TRUE), 
         linetype = guide_legend(reverse = TRUE)) +
  labs(x = "sample size per group") +
  custom_theme(18) 

pp2

ggsave(paste0(fig,"lmb.png"), plot = pp2, width = 15, height = 5, dpi = 300)
##########
pp3 = ggplot(dd_sub, aes(sample_size,power, color=dataset)) +
  geom_point(size=3) +
  geom_line(lwd=1.5) + 
  ylim(0,1) +  
  scale_color_manual(name = TeX("$\\log_2$(mean count)"),values = oka_col) +
  facet_wrap(~lmb, scale="free") + 
  guides(color = guide_legend(reverse = TRUE), 
         linetype = guide_legend(reverse = TRUE)) +
  labs(x = "sample size per group") +
  custom_theme(18) 
################################################################################
mean=lapply(filtered_otu_list,function(x){log2(range(rowMeans(x)))})
log2(as.numeric(mean))
#' After normalising and filtering, these are the range 
#5 to 10
sample_vec = seq(50,100,10)#c(50,70,100)
target  =   seq(0.1, 0.99, 0.05); len = length(target)
lmb_    =   10

dd = foreach(k = 1:7,.combine = "rbind") %:%
   foreach(i = 1:length(sample_vec),.combine = "rbind") %do%{
    
    nam         =   names(filtered_otu_list)[k]
    gam_fit_ss  =   gam_mod_list[[k]]
    
            ss  =    rep(sample_vec[i],len)
          lmb   =    rep(lmb_, len)
      pp        =    inverse_lfc(target=target,samp_size=ss,
                                  lmb   =  lmb, model = gam_fit_ss,
                                  xmin  =  0, xmax = 10)
          d     =    data.frame(lfc = pp, lmb=lmb, sample_size = ss, 
                                power=target)
 d$sample_size  =    factor(d$sample_size, levels = unique(d$sample_size))
          
    d$dataset   =    rep(nam, nrow(d))
    d
   }

################################################################################
oka_col = c(
  "#556B2F", 
  "#E69F00", # Orange
  "#0000FF",
  "#E23D28"
)

dd_sub  =  dd[dd$power >= 0.8,]
dd_sub$power  = factor(dd_sub$power, levels = unique(dd_sub$power)) 
(pp3 = ggplot(dd_sub, aes(sample_size, lfc, group=power,color=power)) +
  geom_point(size=3) +
  geom_line(lwd=1.5) +
  scale_color_manual(name = "power",values = oka_col) +
  facet_wrap(~dataset, scale="free") + 
  guides(color = guide_legend(reverse = TRUE), 
         linetype = guide_legend(reverse = TRUE)) +
  labs(x = "sample size per group", y = TeX("$\\log_2$(fold change)") ) +
  custom_theme(18) )


ggsave(paste0(fig,"effect_est.png"), plot = pp3, width = 15, height = 10, dpi = 300)


max_values <- dd_sub %>%
  group_by(sample_size,dataset,power) %>%
  filter(lfc == max(lfc))  
View(max_values)


for(i in sample_vec){
  p=max_values[max_values$sample_size == i,]
}

range_values <- max_values %>%
  group_by(sample_size,dataset,power) %>%
  filter(lfc == range(lfc)[1])  
View(range_values)

#' for sample size of 30, we need effect size of at least 5 to be able to detect
#' an 80% power, which is every unlikely. 
#' 
ggsave(paste0(fig,"lfc_power.png"), plot = pp3, width = 15, height = 5, dpi = 300)
################################################################################ 
###predicting log mean count is not the best
###
sample_vec  =   c(50,70,100)
target      =   seq(0.1, 0.99, 0.05); len = length(target)
abs_lfc_    =   2

dd = foreach(k = 1:7,.combine = "rbind") %:%
  foreach(i = 1:length(sample_vec),.combine = "rbind") %do%{
    
    nam         =   names(filtered_otu_list)[k]
    gam_fit_ss  =   gam_mod_list[[k]]

    ss        =    rep(sample_vec[i],len)
    abs_lfc   =    rep(abs_lfc_, len)
    pp        =    inverse_lmb(target=target,samp_size=ss,
                               abs_lfc   =  abs_lfc, model = gam_fit_ss,
                               xmin  =  -5, xmax = 5)
    d     =    data.frame(lmb = pp, abs_lfc=abs_lfc, sample_size = ss, 
                          power=target)
    d$sample_size  =    factor(d$sample_size, levels = unique(d$sample_size))
    
    d$dataset   =    rep(nam, nrow(d))
    d
  }


oka_col = c(
  "#556B2F", 
  "#E69F00", # Orange
  "#0000FF"
)

View(dd)
(pp7 = ggplot(dd, aes(lmb,power, color=sample_size)) +
    geom_point(size=3) +
    geom_line(lwd=1.5) +
    ylim(0,1) +  
    scale_color_manual(name = "sample size per group",values = oka_col) +
    facet_wrap(~dataset, scale="free") + 
    guides(color = guide_legend(reverse = TRUE), 
           linetype = guide_legend(reverse = TRUE)) +
    labs(x = TeX("$\\log_2$(mean count)")) +
    custom_theme(18) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", lwd = 1.7)) 

ggsave(paste0(fig,"lmb_power.png"), plot = pp3, width = 15, height = 5, dpi = 300)


###############################################################################
lfc  =  c(0.5,1,2); lmb0 = list(c(-1,2,7),c(-1,2,7),c(-1,2,7)) 
target  =   c(seq(0.1, 0.99, 0.05),seq(0.1, 0.99, 0.05),seq(0.1, 0.99, 0.05))
len = length(target)

dd = foreach(k = 1:7, .combine = "rbind") %:%
  foreach(i = 1:length(lfc) ,.combine = "rbind") %do%{
    
    nam         =   names(filtered_otu_list)[k]
    gam_fit_ss  =   gam_mod_list[[k]]
    
    abs_lfc = rep(lfc[i],len)
    lmb     = rep(lmb0[[i]],each=len/length(lfc))
    
    pp   =    inverse_fun(target=target,lmb=lmb,abs_lfc=abs_lfc, 
                          model=gam_fit_ss,xmin=10, xmax=100)
    
    d    =    data.frame(sample_size = pp, lmb=lmb, 
                         lfc = abs_lfc, power=target)

    d$lfc  =    factor(d$lfc, levels = unique(d$lfc))
    d$lmb  =    factor(d$lmb, levels = unique(d$lmb))
    
    d$dataset = rep(nam, nrow(d))
    d
  }

# Subset and order the data
dd_sub = dd %>% 
  filter(power == 0.8) %>% 
  arrange(sample_size)

min_values <- dd_sub %>%
  group_by(lfc,lmb) %>%
  filter(sample_size == min(sample_size)) %>%
  select(lfc, lmb, sample_size = sample_size, dataset)

custom_labels <- as_labeller(c("0.5" = "|log fold change| = 0.5", 
                               "1" = "|log fold change| = 1", 
                               "2" = "|log fold change| = 2"))
# Plot with reordered dataset within each lfc facet

oka_col = c(
  "#556B2F", 
  "#0000FF",
  "#E69F00" # Orange
)

pp4 = ggplot(dd_sub, aes(x = reorder_within(dataset, sample_size, lfc), y = sample_size,
                         group=lmb,color= lmb)) +
  geom_point(size=3) +
  geom_line(lwd=1.5) +
  scale_color_manual(name = TeX("$\\log_2$(mean count)"), values = oka_col) +
  facet_wrap(~lfc, scales = "free",labeller = custom_labels) +
  scale_x_reordered() +
  custom_theme(18)  +
 #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.x = element_blank()) +

  #geom_text(data = min_values, aes(x = -Inf, y = sample_size, 
   #                                label = paste(sprintf("%.2f",sample_size))), 
    #        hjust = -0.7, vjust = 0.5, color = "red", size=6) +
  #labs(x = "data accession number", y ="sample size per group") 
  labs(x = "dataset", y ="sample size per group") 

ggsave(paste0(fig,"ss_plus.png"), plot = pp4, width = 15, height = 5, dpi = 300)
#################################################################################
#what other plot? Expected number of significan oty and sample sizes. . expected number of  
nsim = 1; notu=1000; disp_scale = 0.3
otu_data =  list()
nsamp_vec = seq(30,190,20)

for(i in 1:7){
  nam       =   names(dispersion_param_list)[i] 
  dispersion_param      =   dispersion_param_list[[i]]
  logmean_param         =   logmean_param_list[[i]]
  logfoldchange_param   =   logfoldchange_param_list[[i]]
  
  dd = foreach(k = 1:length(nsamp_vec), .combine="rbind") %do% {
    source("reproducible/power/simulate_fun.R")
    nsamp     =   nsamp_vec[[k]]
    logmean        =   logmean_sim_fun(logmean_param,notu)
    logfoldchange  =   logfoldchange_sim_fun(logmean,logfoldchange_param)
    d  =  data.frame(lmean_abund=logmean, abs_lfc = abs(logfoldchange))
    d$sample_size = rep(nsamp, nrow(d))
    d
  }
  otu_data[[i]] = dd
}
names(otu_data) = names(dispersion_param_list)  

dd = foreach(k = 1:7, .combine = "rbind") %do% {
    nam         =   names(filtered_otu_list)[k]
    gam_fit_ss  =   gam_mod_list[[k]]
    
    newdata     =   otu_data[[i]]
 newdata$power  =    predict(gam_fit_ss, newdata= newdata, type="response")
 newdata$dataset = rep(nam, nrow(newdata))
 newdata
}

min_values <- dd %>%
  group_by(sample_size,dataset) %>%
  dplyr::summarise(average_power=mean(power))

min_values$expected  =  notu*min_values$average_power

(pp6 = ggplot(min_values, aes(x = sample_size, y= expected, color = dataset)) +
  geom_point(size=3) +
  geom_line(lwd=1.5) +
  guides(color = guide_legend(reverse = TRUE), 
         linetype = guide_legend(reverse = TRUE)) +
   #scale_color_manual(name = TeX("$\\log_2$(mean count)"), values = oka_col) +
  custom_theme(18)  +
  labs(x = "sample size per group", y= "expected number of significant taxa") )



ggsave(paste0(fig,"expected.png"), plot = pp6, width = 7, height = 5, dpi = 300)

dd$pval_reject  = ifelse(dd$power < 0.05, "Yes",  "No")
dd <- dd %>%
      mutate(pval_reject = ifelse(pval_reject == 1, "Yes", "No"))

library(metR)
dd$pvalue_reject <- factor(dd$pval_reject)
dd_sub   =    dd[dd$sample_size == 110,]

cont_breaks   = seq(0,1,0.01)
gg_2dimc <- (ggplot(dd_sub)
               + aes(lmean_abund, abs_lfc)
               + rasterise(geom_point(aes(color = pval_reject), alpha = 0.5))
               + xlab(TeX("$\\log_2$(mean counts)")) 
               + ylab(TeX("|$\\log_2$(fold change)|")) 
               + scale_colour_manual(values = c("black", "red"))
               + geom_contour(data = dd_sub,
                              aes(z=power),lwd=1,
                              breaks = cont_breaks)
               + geom_label_contour(data = dd_sub, 
                                    aes(z= power,label = sprintf("%.3f", after_stat(level))),
                                    breaks = cont_breaks
               ) +
               facet_wrap(~dataset,scales = "free")
               
  )
  gg_2dimc

  
   (ggplot(dd_sub, aes(lmean_abund, power)) +
                geom_line() 
               # + rasterise(geom_point(aes(color = pval_reject), alpha = 0.5))
               + xlab(TeX("$\\log_2$(mean counts)")) 
               + ylab(TeX("|$\\log_2$(fold change)|")) +
                 facet_wrap(~dataset,scales = "free")
               
  )
  
  
##' You need at least ... samples per group to aattain 80% or more power. 
##' In scientific research a power of at least 80% is considered enough power, 
##' so we explore how many samples are needed to get a at least 80% power. 
##' To do that we need to decide on the an effect size to use. so you can find 
##' out effect sizes needed to predict at least 80% power for each dataset 
##' so we have a plot for 3 mean abundance, where mean abundance of 5, 10 and 15 
##' and power of 80-100 and then we desice of different sample sizes or show we 
##' have  
##' You can use the mean abundance from the actual datasets
##' simulate lmb and then lf
##' predict lfc, 
##' predict lmb
##' simulate lfc and use it to predict lmb and then use the lmb to predict lfc

combined_data  =  foreach(i = 1:7) %do% {
  nam      =   names(logmean_list)[i]
  dd      =   readRDS(paste0(path,"combined_data/",nam, ".rds"))
  dd
}


for(r in 1:length(target)){
  targ  =   target[r]
  
  pp=uniroot_lmb(target_power=targ,sample_size=50,abs_lfc=5,
                 model=gam_fit_ss,xmin=-5,xmax=15)
  print(pp)
}

 



abs_lfc = rep(lfc,len)
lmb     = rep(lmb_vec[i],len)

lfc    =  2;  lmb_vec = c(-1,0,2,5)
target = seq(0.1, 0.99, 0.05); len = length(target)

dd  = foreach(i = 1:7 ,.combine = "rbind") %do%{
    nam         =   names(filtered_otu_list)[i]
    comb        =   combined_data[[i]]
    lmb         =   comb$lmean_abund
    abs_lfc     =   comb$abs_lfc
    gam_fit_ss  =   gam_mod_list[[i]]
    
    sub_comb  <- comb %>%
      group_by(sample_size) %>%
      slice_sample(n = 50)
    
    comb$power=predict(gam_fit_ss, newdata= comb, type="response")
    
    
    #ggplot(sub_comb, )
    pp   =    inverse_fun(target=target,lmb=lmb,abs_lfc=abs_lfc, 
                          model=gam_fit_ss,xmin=10, xmax=100)
    
    d    =    data.frame(sample_size = pp, lmb=lmb, 
                         lfc = abs_lfc, power=target)
    
    d$lmb  =    factor(d$lmb, levels = unique(d$lmb))
    d$dataset = rep(nam, nrow(d))
    d
  }

#' predict effect and predict lfc with it
sub=c("PRJNA589343", "PRJNA644763","PRJNA687773")
dd_sub = dd[dd$dataset %in% sub,]

ggplot(sub_comb, aes(abs_lfc,power)) +
  geom_point() +
  #geom_line()  + 
  geom_point(aes(lmean_abund,power),color= "blue") 
  #geom_line(aes(lmean_abund,power),color= "blue") 
#   ylim(0,1) #+  
  #scale_color_manual(name = TeX("$\\log_2$(mean count)"),values = oka_col) +
  #facet_wrap(~dataset, scale="free") + 
  guides(color = guide_legend(reverse = TRUE), 
         linetype = guide_legend(reverse = TRUE)) +
  labs(x = "sample size per group") +
  custom_theme(18) 

ggsave(paste0(fig,"lmb.png"), plot = pp2, width = 15, height = 5, dpi = 300)
################################################################################


## logMean vrs LFC plot
logmean <- logmean_list %>%
  imap_dfr(~ data.frame(logmean = .x, Name = paste0(.y, "...")))  


logfoldchange <- logfoldchange_list %>%
  imap_dfr(~ data.frame(logfoldchange = .x, Name = paste0(.y, "...")))  

dd = left_join(logmean,logfoldchange, by= "Name") 

p = ggplot(dd, aes(logmean, logfoldchange, group = Name)) +
  geom_point(alpha=0.1)  + 
  geom_smooth() + 
  facet_wrap(~Name)

p
  # geom_text_repel(aes(label = label),
  #                 max.overlaps = nrow(dat)) + 
  #xlab(TeX("$\\alpha$")) 
  xlab(TeX("$\\log_2$(mean abundance)")) + ylab(TeX("$\\log_2$(foldchange)")) +
  ggtitle(Names) +
  theme_bw(base_size=14) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),  
        axis.text.y = element_text(color = "black") 
  )
p


plt <- foreach(i = 1:7) %do% {
  
  dat   =   data.frame(logmean = logmean_list[[i]], 
                      logfoldchange = logfoldchange_list[[i]])
  
  p = ggplot(dat, aes(logmean, logfoldchange)) +
    geom_point(alpha=0.1)  + 
    geom_smooth() + 
    # geom_text_repel(aes(label = label),
    #                 max.overlaps = nrow(dat)) + 
    #xlab(TeX("$\\alpha$")) 
    xlab(TeX("$\\log_2$(mean abundance)")) + ylab(TeX("$\\log_2$(foldchange)")) +
    #ggtitle(data_name) +
    theme_bw(base_size=14) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(color = "black"),  
          axis.text.y = element_text(color = "black") 
    )
  p
  

}


plt[[2]]|plt[[5]]|plt[[3]] +  plot_layout(guides = "collect")  



lmean_lfoldchange_plt <- function(logmean, logfoldchange,data_name="XYZ"){
  
  dat <- data.frame(logmean, logfoldchange, label = names(logmean))
 
}



######################################################################
### Contour plot
metadata_list   =  readRDS(file = paste0(path,"metadata.rds"))
combined_data = power_data= list()
for(i in 1:7){
  #i =1
  nam =  names(metadata_list)[i]
  comb =  readRDS(paste0(path, "contour_data/100_samples/combined_data_", nam,"_.rds"))
  pow_data  =  readRDS(paste0(path, "contour_data/100_samples/power_estimate_", nam,"_.rds"))
  
  comb$dataset = rep(nam, nrow(comb))
  pow_data$dataset     =    rep(nam, nrow(pow_data))
  combined_data[[i]]  = comb
  power_data[[i]]     = pow_data
}

names(power_data) = names(combined_data) = names(metadata_list)

comb_data=do.call(rbind, combined_data)
pow_data =do.call(rbind, power_data)
comb_data$pval_reject <- factor(comb_data$pval_reject, 
                                levels= unique(comb_data$pval_reject))
comb_data <- comb_data %>%
  mutate(pval_reject = ifelse(pval_reject == 1, "Yes", "No"))

dim(comb_data)
#comb_data$pval_reject <- factor(comb_data$pval_reject, levels = c("Yes", "No"))
sampled_df <- comb_data %>% sample_n(10000, replace = FALSE)

View(sampled_df)
cont_breaks =  c(seq(0,0.1,0.04),seq(0.4,1,0.2))
gg_2dimc <- (ggplot(sampled_df)
             + aes(lmean_abund, abs_lfc)
             + rasterise(geom_point(aes(color = pval_reject), alpha = 0.5))
             + xlab(TeX("$\\log$(mean count)")) 
             + ylab(TeX("|$\\log$(fold change)|")) 
             + scale_colour_manual(values = c("black","red"))
            + geom_contour(data = pow_data,
                           aes(z=power),lwd=1,
                           breaks = cont_breaks)
            + geom_label_contour(data = pow_data,
                                 aes(z= power,label = sprintf("%.2f", after_stat(level))),
                                 breaks = cont_breaks)
            #+ facet_wrap(~dataset, scales = "free_y")
           + labs(colour = "significant taxa")  +
            custom_theme(15)  

)
gg_2dimc

ggsave(paste0(fig_path,"contour_plots.png"))

######################################################################
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

ggsave(paste0(fig,"cv.png"))

ggplot(disp, aes(x = data_accession_no, y= dispersion)) +
  geom_violin(scale = "width") +
  geom_point(alpha= 0.1) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size =10, family = "Roboto"),
    axis.text.x = element_text(family = "Roboto",angle = 90, size =10, color = "black"),
    axis.text.y = element_text(family = "Roboto", size =10, color = "black")
  )+ 
  labs(x="data accession number", y="dispersion estimate")

ggsave(paste0(fig_path,"dispersion.png"))
######################################################################
disp = foreach(i = 1:7,.combine = "rbind") %do%{
  nam =  names(logmean_list)[[i]]
  x = 2^logmean_list[[i]]
  #dispersion  =   dispersion_list[[i]]
  
  data        =   normalised_otu_list[[i]]
  var_cal     =   apply(data,1,function(x){var(x)})
  
  dispersion    =   (var_cal - x)/x^2
    
  asymptDisp  =  dispersion_param_list[[i]]$asymptDisp
  extraPois   =  dispersion_param_list[[i]]$extraPois
  disp_sim    =  dispersion_list[[i]]#dispersion_fun(x, asymptDisp,extraPois)
  
  
  cv_true  =  coeff_var(mean=x, dispersion=dispersion)
  cv_sim   =  coeff_var(mean=x, dispersion=0.3*disp_sim)
  cv = c(cv_true,cv_sim)

  dd =  data.frame(cv  =  cv, 
                   cv_type = c(rep("cv_true", length(cv_true)),
                              rep("cv_sim", length(cv_sim))), 
                   dataset = rep(nam,length(cv)))
  
  dd
}

ggplot(disp, aes(x = dataset, y= cv)) +
  geom_violin(scale = "width") +
  geom_point(alpha= 0.1) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size =10, family = "Roboto"),
    axis.text.x = element_text(family = "Roboto",angle = 90, size =10, color = "black"),
    axis.text.y = element_text(family = "Roboto", size =10, color = "black")
  )+ 
  labs(x="data accession number", y="dispersion estimate") +
  facet_wrap(~ cv_type)







######################################################################
gg_2dimc = list();   cont_breaks = c(0.4)

a_dataframes <-  lapply(contour_data, function(x) x$"100samples")
combined_a   <-  do.call(rbind, a_dataframes)
 
for(n  in   1:7){
  combined_dd         =   contour_data[[n]][["100samples"]]
  power_estimate_dd   =   power_estimate[[n]][["100samples"]]

  contour_plot_fun(combined_dd, 
                   power_estimate_dd, 
                               cont_breaks=seq(0,1,0.1))
  
  
  
  combined_dd         =    do.call(rbind,combined_dd)
  power_estimate_dd   =    do.call(rbind,power_estimate_dd)
  
  
  combined_dd$pvalue_reject      <- factor(combined_dd$pval_reject)
  combined_dd$sample_size        <- factor(combined_dd$sample_size)
  power_estimate_dd$sample_size  <-  factor(power_estimate_dd$sample_size)
  
  gg_2dimc[[n]] <- (ggplot(combined_dd)
                    + aes(lmean_abund, abs_lfc)
                    #+ geom_point(pch=".",alpha = 0.5)
                    + rasterise(geom_point(alpha = 0.5))
                    + xlab(TeX("$\\log_2$(mean counts)")) 
                    + ylab(TeX("|$\\log_2$(fold change)|")) 
                    + geom_contour(data = power_estimate_dd,
                                   aes(z=power,color =sample_size,group=sample_size),lwd=1.5,
                                   breaks = cont_breaks)
                    + geom_label_contour(data = power_estimate_dd, 
                                         aes(z= power,color =sample_size,group=sample_size,
                                             label = sprintf("%.3f", after_stat(level))),
                                         breaks = cont_breaks
                    ) 
                    + scale_colour_manual(values = c("blue", "red",  "green", "orange"))
                    
  )
  
}
names(gg_2dimc)  =  names(contour_data)
######################################################################
gg_2dimc = list();   cont_breaks = c(0.4)

for(n  in   1:7){
  combined_dd         =   contour_data[n][[1]]
  power_estimate_dd   =   power_estimate[n][[1]]
  combined_dd         =    do.call(rbind,combined_dd)
  power_estimate_dd   =    do.call(rbind,power_estimate_dd)
  
  
  combined_dd$pvalue_reject      <- factor(combined_dd$pval_reject)
  combined_dd$sample_size        <- factor(combined_dd$sample_size)
  power_estimate_dd$sample_size  <-  factor(power_estimate_dd$sample_size)
  
  gg_2dimc[[n]] <- (ggplot(combined_dd)
               + aes(lmean_abund, abs_lfc)
               #+ geom_point(pch=".",alpha = 0.5)
               + rasterise(geom_point(alpha = 0.5))
               + xlab(TeX("$\\log_2$(mean counts)")) 
               + ylab(TeX("|$\\log_2$(fold change)|")) 
               + geom_contour(data = power_estimate_dd,
                              aes(z=power,color =sample_size,group=sample_size),lwd=1.5,
                              breaks = cont_breaks)
               + geom_label_contour(data = power_estimate_dd, 
                                    aes(z= power,color =sample_size,group=sample_size,
                                        label = sprintf("%.3f", after_stat(level))),
                                    breaks = cont_breaks
               ) 
               + scale_colour_manual(values = c("blue", "red",  "green", "orange"))
               
  )
  
}

names(gg_2dimc)  =  names(contour_data)
  
(gg_2dimc[[1]]|gg_2dimc[[2]]|gg_2dimc[[3]]) +  plot_layout(guides = "collect")  

#contour plot 



path = "reproducible/power/datasets/"
countdata_list_obs  =  readRDS(file = paste0(path,"data.rds"))
metadata_list_obs   =  readRDS(file = paste0(path,"metadata.rds"))
######################################################
combined_data_list = power_estimate_list = list()
nsamp  =  100
for(i in 1:7){
  nam  =   names(metadata_list_obs)[i]
  combined_data     =   readRDS(paste0(path,"combined_data_",nam,"_", nsamp,".rds"))
  power_estimate    =   readRDS(paste0(path,"power_estimate_",nam,"_", nsamp,".rds"))
  
  combined_data$dataset =   rep(nam, nrow(combined_data))
  power_estimate$dataset =   rep(nam, nrow(power_estimate))
  
  combined_data_list[[i]]       =  combined_data
  power_estimate_list[[i]] =  power_estimate

}
names(combined_data_list)   =   names(metadata_list_obs)
names(power_estimate_list)  =   names(metadata_list_obs)

combined_df       <- do.call(rbind, combined_data_list)
power_estimate_df <- do.call(rbind, power_estimate_list)

combined_df$pvalue_reject <- factor(combined_df$pval_reject)



#Contour plot
 

