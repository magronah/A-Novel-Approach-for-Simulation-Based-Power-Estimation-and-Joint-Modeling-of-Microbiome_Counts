setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")
library(Matrix)
library(huge)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
source("func2.R")
source("initial_param0.R")
####################################################################
ntaxa  = 600 ; nsubj = 200
path   =   paste0("cov2/",nsubj,"_",ntaxa, "/")
path
####################################################################
## combined p-value calculation
nbmm_pval     =   readRDS(paste0(path,"nbmm_pval.rds"))
zinbmm_pval   =   readRDS(paste0(path,"zinbmm_pval.rds"))

nbmm_pval <- as.data.frame(lapply(nbmm_pval, function(col) as.numeric(unlist(col))))
zinbmm_pval <- as.data.frame(lapply(zinbmm_pval, function(col) as.numeric(unlist(col))))


combine_nbmm1    =  apply(nbmm_pval, 1, fishers_combined)
combine_nbmm2    =  apply(nbmm_pval, 1, stouffer_combine)

plot(log(combine_nbmm1 +  1e-6),log(combine_nbmm2 +  1e-6))
abline(0,1)

combine_zinbmm1    =   apply(zinbmm_pval, 1, fishers_combined)

summary(combine_zinbmm1[-103])



combine_zinbmm2    =  apply(zinbmm_pval, 1, stouffer_combine)

which(is.na(combine_zinbmm2), arr.ind = TRUE)

plot(log(combine_zinbmm1 +  1e-6),log(combine_zinbmm2 +  1e-6))
abline(0,1)


range(combine_zinbmm2)
#############################################################################
ggplot(dd, aes(x = taxa)) +
  geom_density() +
  labs(title = "Density Plots for First 50 Columns", 
       x = "Value", y = "Density") +
  theme_minimal()

deseq_dd   =   readRDS(paste0(path,"nbmm.rds"))
dd = as.data.frame(t(deseq_dd))

# Load the necessary library
library(boot)
i=100
set.seed(123)
bootstrap_samples <- as.numeric(deseq_dd[i,])
param_function <- function(data, indices) {
  return(mean(data[indices]))
}

bootstrap_results <- boot(data=bootstrap_samples, statistic=param_function, R=1000)

bca_ci <- boot.ci(bootstrap_results, type="bca")
point_estimate <- bootstrap_results$t0  # Your combined estimate

print(point_estimate)
print(paste("BCa Confidence Intervals: ", bca_ci$bca[4], " to ", bca_ci$bca[5]))


set.seed(123)

# Calculate the mean and variance
mean_estimate <- mean(bootstrap_samples)
variance_estimate <- var(bootstrap_samples)
n <- length(bootstrap_samples)

# James-Stein shrinkage factor
shrinkage_factor <- 1 - (variance_estimate / sum((bootstrap_samples - mean_estimate)^2))
shrinkage_factor <- max(shrinkage_factor, 0)

# Compute the James-Stein estimator
james_stein_estimate <- shrinkage_factor * bootstrap_samples + (1 - shrinkage_factor) * mean_estimate
combined_estimate_james_stein <- mean(james_stein_estimate)

# Display results
print(paste("James-Stein Combined Estimate: ", combined_estimate_james_stein))



dd_long <- pivot_longer(dd, cols = everything(), 
                        names_to = "taxa", 
                        values_to = "estimate")

dd_long_filtered <- dd_long %>% 
  filter(taxa %in% unique(taxa)[1:50])

ggplot(dd_long_filtered, aes(x = estimate)) +
  geom_density() +
  facet_wrap(~ taxa, scales = "free") +
  labs(title = "Density Plots for First 50 Columns", 
       x = "Value", y = "Density") +
  theme_minimal()


# View the result
head(dd_long,10)


colnames(dd)
colnames(t(deseq_dd))
length(deseq_dd[1,])

hist(as.numeric(deseq_dd[100,]))


deseq_dd   =   readRDS(paste0(path,"deseq_otu_meta_list_withzi_taxa.rds"))
deseq_noShrink_dd   =   readRDS(paste0(path,"deseq_noShrink_otu_meta_list_withzi_taxa.rds"))
#################################
#GAM fit

r =  readRDS(path)
gam_fit(pvalue,effect_size,mean_count,grid_len = 100, alpha_level = 0.05)

library(boot)
set.seed(123)
data <- rnorm(100, mean=50, sd=10)

param_function <- function(data, indices) {
  return(mean(data[indices]))
}

bootstrap_results <- boot(data=data, statistic=param_function, R=1000)

# Get BCa confidence intervals and point estimate
bca_ci <- boot.ci(bootstrap_results, type="bca")
point_estimate <- bootstrap_results$t0  # Your combined estimate

# Display results
print(point_estimate)
print(bca_ci)
