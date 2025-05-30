library(tidyr)
library(tidyverse)
library(tibble)
library(here)
library(DESeq2)
library(RhpcBLASctl)
library(glmmTMB)
##############################################################
path   =   paste0(getwd(),"/real_data/")
source(paste0(path,"/fun.R"))
##############################################################
data        =   readRDS(paste0(path,"aut_data.rds"))
metadata    =   readRDS(paste0(path,"aut_metadata.rds"))

nam        =  "PRJNA644763"
count_dd   =   data[[nam]]
meta_dd    =   metadata[[nam]] %>% 
               setNames(c("subject", "group"))

dd_filt      =   filter_fun(count_dd,meta_dd,abund_thresh=10, 
                            sample_thresh=5)
##############################################################
dim(dd_filt)
dim(count_dd)
##############################################################
dd_long     =   df_long(dd_filt)
dd_long     =   left_join(dd_long, meta_dd, by ="subject")  
###########################################################
pp   =  deseqfun(dd_filt,meta_dd,alpha_level=0.1,ref_name="ASD")

normalizer =   data.frame(normalizer = sizeFactors(pp$object)) %>% 
  rownames_to_column("subject")

df      =    left_join(dd_long, normalizer, by ="subject")  
###########################################################
form      =   count ~ 1 + us(1 + group|taxon) +  
                     rr(0 + taxon | subject,2) + 
                     offset(normalizer)

##Fit the model
par_ctrl <- glmmTMBControl(
  parallel = list(n = 10, autopar = TRUE)
)

#par_ctrl <- glmmTMBControl(
#  parallel = list(n = 6, autopar = TRUE),
#  optArgs  = list(maxit = 10000)
#)

gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")

options(glmmTMB_openmp_debug = TRUE)
blas_set_num_threads(1)

tt1 = system.time(
    fit1  <-  glmmTMB(form, data = df,
                   family  = nbinom2, 
                   prior   = gprior,
                   REML    = FALSE,
                   control = par_ctrl
  )
)

saveRDS(tt1, file=paste0(path,"autism/data/rr_runtime",nam,".rds"))
saveRDS(fit1, file=paste0(path,"autism/data/rr_mod",nam,".rds"))


tt2 = system.time(
  fit2 <- glmmTMB(formula = form, 
                 data = df, 
                 family = nbinom2, 
                 ziformula = ~1, # + (1 | taxon),
                 prior = gprior, 
                 REML = FALSE, 
                 control = par_ctrl)
)
saveRDS(tt2, file=paste0(path,"autism/data/rrzi_runtime",nam,".rds"))
saveRDS(fit2, file=paste0(path,"autism/data/rrzi_mod",nam,".rds"))
###############################################################
form      =   count ~ 1 + us(1 + group|taxon) +    offset(normalizer)

tt3 = system.time(
  fit3  <-  glmmTMB(form, data = df,
                    family  = nbinom2, 
                    prior   = gprior,
                    REML    = FALSE,
                    control = par_ctrl
  )
)

saveRDS(tt3, file=paste0(path,"autism/data/us_runtime",nam,".rds"))
saveRDS(fit3, file=paste0(path,"autism/data/us_mod",nam,".rds"))


tt4 = system.time(
  fit4 <- glmmTMB(formula = form, 
                  data = df, 
                  family = nbinom2, 
                  ziformula = ~1,# + (1 | taxon),
                  #prior = gprior, 
                  REML = FALSE, 
                  control = par_ctrl)
)


saveRDS(tt4, file=paste0(path,"autism/data/uszi_runtime",nam,".rds"))
saveRDS(fit4, file=paste0(path,"autism/data/uszi_mod",nam,".rds"))
##############################################################
tt5 = system.time(
  fit5 <- glmmTMB(formula = form, 
                  data = df, 
                  family = nbinom2, 
                  ziformula = ~1 + (1 | taxon),
                  #prior = gprior, 
                  REML = FALSE, 
                  control = par_ctrl)
)

saveRDS(tt5, file=paste0(path,"autism/data/uszi_each_runtime",nam,".rds"))
saveRDS(fit5, file=paste0(path,"autism/data/uszi_each_mod",nam,".rds"))
#####################################################################
# Transform data
otu_clr <- microbiome::transform(dd_filt, "clr")

# Compute correlation
cor_matrix <- cor(t(otu_clr), method = "spearman")

# Visualize
library(reshape2)
library(ggplot2)
cor_df <- melt(cor_matrix)

ggplot(cor_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = "black", mid = "yellow", high = "blue") +
  labs(x = "ASV", y = "ASV") +
  theme_bw()


betazi    =  fit$fit$par[["betazi"]]
betadisp  =  
theta     =  fit$fit$par[names(fit$fit$par) == "theta"]

unique(names(fit$fit$par))
#################################
# Generate data
ntaxa    =  max(dim(dd_filt))
grp_pars =  list(name_vec =  as.vector(unique(met_dd$group)), 
                size_vec  =  c(10,10))

param    =   list(ntaxa = 100, 
                  beta_pars =  list(beta= fixef(fit),
                                    betadisp = fit$fit$par[["betadisp"]]),
                  theta     =  theta[names(fit$fit$par) == "theta"],
                  grp_pars  =  list(name_vec = c("control", "treat"), 
                                  size_vec = c(10,10)))

                  

meta_dd <- met_data(ntaxa, grp_pars)

# Simulate new data using the specified parameters
sim_count = simulate_new(
  RHSForm(form, as.form = TRUE),
  newdata =   meta_dd,
  newparams = pars,
  family = nbinom2,
  ziformula = ~ taxon,
  nsim = nsim,
  seed = seed
)


names(fit$fit$par)
#next try the other datasets as well 
#coef(fit)
# library(glmmTMB)
# setting OpenMP: threads = 6L, autopar = TRUE
# setting OpenMP: threads = 1L, autopar = FALSE
# user   system  elapsed 
# 2070.020    6.624  391.032

# setting OpenMP: threads = 6L, autopar = TRUE
# setting OpenMP: threads = 1L, autopar = FALSE
# user   system  elapsed 
# 1743.617    6.990  331.135 

#############################################
# Assuming `otu_table` is your OTU table where rows are taxa and columns are samples

# Calculate the proportion of zeros for each taxon
proportion_zeros <- rowSums(dd_filt == 0) / ncol(dd_filt)

# Convert to a data frame for ggplot2
proportion_zeros_df <- data.frame(Taxon = rownames(dd_filt), ProportionZeros = proportion_zeros)

# Create the histogram
library(ggplot2)
n = 13
ggplot(proportion_zeros_df, aes(x = ProportionZeros)) +
  geom_histogram(binwidth = 0.005, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    #title = "Histogram of Proportion of Zeros in Each Taxon",
    x = "Proportion of Zero Counts",
    y = "Number of Taxa"
  ) +
  custom_theme(n)

ggsave(paste0("zero_prob.png"))  # Save the plot


 

