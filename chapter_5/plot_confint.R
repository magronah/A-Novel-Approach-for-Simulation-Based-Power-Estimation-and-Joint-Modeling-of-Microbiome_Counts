library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(here)
source("func2.R")
fig_path=   "fig/"
########################################################
path    =   "200_600/confint/"
titles  =   "100 subjects per group and 600 taxa"
########################################################
# Define filenames to load
filenames <- c(
  "rr.rds",
  "rrzi.rds",
  "us.rds",
  "uszi.rds",
  "nbmm.rds",
  "zinbmm.rds",
  "deseq.rds"  
)

# Load models for autism data
dd        <-  load_models(path, filenames)
names(dd)  =    sub("\\.rds$", "", filenames)
names(dd)  <- c("RR","RRzi","US","USzi","NB","ZNB","DE" )
#####################################################################
nrw   =   600;  ncl = 6
dd_nbmm     <-  discard(dd[["NB"]],   ~ any(dim(.x) != c(nrw, ncl)))
dd_zinbmm   <-  discard(dd[["ZNB"]], ~ any(dim(.x) != c(nrw, ncl)))
######################################################################
common_names1  <-   intersect(names(dd_nbmm), names(dd_zinbmm))
common_names   <-   intersect(common_names1, names(dd$RRzi))

dd_nbmm       <-   dd_nbmm[common_names]
dd_zinbmm     <-   dd_zinbmm[common_names]
dd_deseq      <-   dd$DE[common_names]
dd_rr         <-   dd$RR[common_names]
sub_dd        <-   dd[grep("US|RR", names(dd))]
us_rr         <-   lapply(sub_dd, function(x){x[common_names]})
dd_list       <-   list(NB  =  dd_nbmm, ZNB = dd_zinbmm)
######################################################################
## compute the average proportion of significant taxa
sig1   <-   lapply(us_rr, function(y){x <- map_dfc(y, ~ .x$pvalue)
                                    apply(x, 2, p.adjust, method = "BH")})

sig2   <-   lapply(dd_list, function(y){x <- map_dfc(y, ~ .x$pvalue)
                           apply(x, 2, p.adjust, method = "BH")})

sig3   <-   list(DE = map_dfc(dd_deseq, ~ pull(.x, padj)))

sig   <-    c(sig1, sig2, sig3)
######################################################################
alpha  = 0.05
ddf1   = data.frame(lapply(sig, function(x){mean(colMeans(x < alpha))}))
ddf11  = data.frame(lapply(sig, function(x){
  sd(colMeans(x < alpha))/length(colMeans(x < alpha))}))
dddf  =  data.frame(power = t(ddf1))
dddf$se  =  t(ddf11)

ddf =   dddf %>%
       data.frame() %>%
       rownames_to_column("model") %>%
       mutate(upr =  power + 1.96*se,
              lwr =  power - 1.96*se)
##############################################################
size = 3; nn= 13; width =  7; height =  5; dpi = 300
plt0 = ggplot(ddf, aes(model, power)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +  
  custom_theme(18) +
  labs(x =" ", y = "Average number of p-values less than 0.05")

ggsave("fig/pvalhits_rr.png", plot = plt0, 
       width = width, 
       height = height, 
       dpi = dpi)
##############################################################
##confidence width and coverage plot
true_param    =   readRDS("200_600/true_param.rds")
##################################################################
cov_list  <-   c(dd_list, us_rr, DE = list(dd_deseq))
##################################################################
ddd <- lapply(cov_list, function(x){
  ddf <- lapply(x, function(y){
    df         =   left_join(true_param,y, by = "param_name")
    df$cov     =   ifelse(df$lwr  < df$true_param & df$true_param < df$upr,1, 0)
    df
  })
  names(ddf)  = names(x)
  ddf
})
##############################################################
covv   <-   lapply(ddd, function(y){map_dfc(y, ~ .x$cov)})
dd_cov   <-   data.frame(lapply(covv, function(x){mean(rowMeans(x))}))
dd_cov2   <-   data.frame(lapply(covv, function(x){(rowMeans(x))}))

dd_covt  =  data.frame(coverage=t(dd_cov))
dd_covt$se <-  unlist(apply(dd_cov2, 2, function(x){sd(x) / sqrt(length(x))}))

ddf1 =  dd_covt %>%
  as.data.frame() %>%
  rownames_to_column("model") %>%
  mutate(lwr  = coverage - 1.96*se,
         upr = coverage  + 1.96*se)

# ddf1 =  t(dd_cov) %>%
#    as.data.frame() %>%
#   setNames("coverage") %>%
#   rownames_to_column("model") %>%
#   mutate(lwr  = coverage - 1.96*se)


plt1 =  ggplot(ddf1, aes(model, coverage)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1) +  
  custom_theme(18) +
  labs(x =" ", y = "coverage") +
  ggtitle("Average coverage for all taxa")

plt1 = ggplot(ddf1, aes(model, coverage,
                        ymin = lwr, ymax = upr)) +
  geom_pointrange(position = position_dodge(width = 1),
                  size = 0.5) +
  geom_point() +
  custom_theme(18) +
  labs(x =" ", y = "coverage") +
  ggtitle("Average coverage for all taxa")


ggsave("fig/coverage_rr.png", plot = plt1, 
       width = width, 
       height = height, 
       dpi = dpi)

#########################################################################
width_dfs <- lapply(cov_list, function(model) {
  width_df <- as.data.frame(lapply(model, function(sim) pull(sim, width)))
  data.frame(wdmean =  mean(rowMeans(width_df)),
             se     = sd(rowMeans(width_df))/length(rowMeans(width_df)))
})

width_dd  =  do.call(rbind, width_dfs) %>%  
             data.frame()  %>%
             rownames_to_column("model")  %>%
             mutate(lwr  = wdmean - 1.96*se,
                     upr = wdmean  + 1.96*se) 
            
plt2 = ggplot(width_dd, aes(model, wdmean,
                        ymin = lwr, ymax = upr)) +
  geom_pointrange(position = position_dodge(width = 0.1), 
                  size = 0.5) +
  geom_point() +
  custom_theme(18) +
  ggtitle("Average confidence width across simulations for all taxa") +
  labs(x =" ", y = "Confidence width")

ggsave("fig/cov_width_rr.png", plot = plt2, 
       width = width, 
       height = height, 
       dpi = dpi)
##################################################################
pl=plt2|plt1
ggsave("fig/cov_nwidth_rr.png", plot = pl, 
       width = 13, 
       height = height, 
       dpi = dpi)

