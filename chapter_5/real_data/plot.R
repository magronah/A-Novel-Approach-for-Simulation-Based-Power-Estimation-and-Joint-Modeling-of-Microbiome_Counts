library(glmmTMB)
library(ggplot2)
library(dplyr)
library(patchwork)
library(AICcmodavg)
library(here)
library(tidyverse)
library(purrr)
##############################################################
path1   =   paste0("real_data/")
source(paste0(path1,"fun.R"))
#########################################################
# Define filenames to load
filenames <- c(
  "rr_mod.rds",
  "rrzi_each_mod.rds",
  "us_mod.rds",
  "uszi_each_mod.rds",
  "nbmm_aicc.rds",
  "zinbmm_aicc.rds",
  "deseq_aicc.rds"  
)

# Define filenames to load
filenames <- c(
  "rr_runtime.rds",
  "rrzi_each_runtime.rds",
  "us_runtime.rds",
  "uszi_each_runtime.rds",
  "nbmm_runtime.rds",
  "zinbmm_runtime.rds",
  "deseq_runtime.rds"  
)
######################################################################
# Load models for autism data
autism_path <- paste0("real_data/autism_data/results/")
atlass_path <- paste0("real_data/atlass_data/results/")
crohn_path  <- paste0("real_data/CrohnD_data/results/")
soil_path   <- paste0("real_data/soil_data/results/")

# Load models for autism data
autism_models <- load_models(autism_path, filenames)
atlass_models <- load_models(atlass_path, filenames)
crohn_models <- load_models(crohn_path, filenames)
soil_models   <- load_models(soil_path, filenames)

# Assigning names 
mod_names <- c("RR","RRzi","US","USzi","NB","ZNB","DE")
names(autism_models)    =   names(atlass_models)  =    mod_names
names(crohn_models)    =   names(soil_models)    =    mod_names
mod_list    =   lst(autism_models,atlass_models,crohn_models,soil_models)
######################################################################
time_atlass <-  data.frame(time = unlist(lapply(atlass_models, 
                                                function(x){x[["elapsed"]]/60})), 
                           model =  names(atlass_models)) %>%
  mutate(dataset = "Human Intestine Data")

time_autism <-  data.frame(time = unlist(lapply(autism_models, 
                                                   function(x){x[["elapsed"]]/60})), 
                              model =  names(autism_models)) %>%
  mutate(dataset = "Autism Data")

time_soil <-  data.frame(time = unlist(lapply(soil_models, 
                                                function(x){x[["elapsed"]]/60})), 
                           model =  names(soil_models)) %>%
  mutate(dataset = "Soil Data")

time_crohn <-  data.frame(time = unlist(lapply(crohn_models, 
                                                function(x){x[["elapsed"]]/60})), 
                           model =  names(crohn_models)) %>%
  mutate(dataset = "Crohn Disease Data")


time_dd  <-  rbind(time_atlass, time_crohn, time_soil, time_autism)

plt_time <- ggplot(time_dd, aes(x =  model, y =  log(time))) + 
  geom_point(size=2) +
  custom_theme(16) +
  labs(x = " ", y = "Runtime (mins) (log-scale)") +
  facet_wrap(~dataset, scales = "free") 

size = 3; nn= 13; width =  13; height =  7; dpi = 300

ggsave("fig/time_rr.png", plot = plt_time, 
       width = width, 
       height = height, 
       dpi = dpi)
######################################################################
##AIC comparison
caic_files <- c(US   = "caic_us.rds",
                USzi = "caic_uszi.rds",
                RR   = "caic_rr.rds",
                RRzi = "caic_rrzi.rds") 

caic_files2 <- c(US   = "caic_us.rds",
                USzi = "caic_uszi.rds") 

other_aic <- c(NB   = "nbmm_aicc.rds",
               ZNB  = "zinbmm_aicc.rds",
               DE   = "deseq_aicc.rds") 
#######################################################
# Load models for autism data
atlass_caic <- load_models(atlass_path, caic_files)
crohn_caic  <- load_models(crohn_path, caic_files)

autism_caic <- load_models(autism_path, caic_files2)
soil_caic   <- load_models(soil_path, caic_files2)

atlass_other <- load_models(atlass_path, other_aic)
crohn_other  <- load_models(crohn_path, other_aic)
autism_other <- load_models(autism_path, other_aic)
soil_other  <- load_models(soil_path, other_aic)
#######################################################
names(autism_caic) = names(soil_caic)= names(caic_files2)
names(crohn_caic) = names(atlass_caic)  = names(caic_files)

names(autism_other) = names(atlass_other) = names(other_aic)
names(crohn_other) = names(soil_other) = names(other_aic)

caic_list     =  lst(autism_caic, atlass_caic, crohn_caic, soil_caic)
#caic_list     =  lst(atlass_caic, crohn_caic)

otheraic_list =  list(autism_caic =  autism_other, 
                     atlass_caic  =  atlass_other, 
                     crohn_caic   =  crohn_other,
                     soil_caic    =  soil_other )
#######################################################
caic_df <- caic_list %>%
  map_dfr(~map_df(.x, ~tibble(Name = names(.x), Estimate = .x), 
                  .id = "Model"), .id = "Dataset") %>%
  mutate(Dataset = case_when(
    Dataset == "autism_caic"  ~ "Autism data",
    Dataset == "atlass_caic"  ~ "Human Intestine data",
    Dataset == "crohn_caic"   ~ "Crohn’s Disease data",
    Dataset == "soil_caic"    ~ "Soil data"
    #TRUE ~ Dataset  # default case, in case other datasets exist
  ))

caic_dd <- caic_df  %>% filter(Name == "caic")
other_aicc_dd <- otheraic_list %>%
              map_dfr(~map_df(.x, ~tibble(Name = names(.x), Estimate = .x), 
                              .id = "Model"), .id = "Dataset") %>%
  mutate(Dataset = case_when(
    Dataset == "autism_caic"  ~ "Autism data",
    Dataset == "atlass_caic"  ~ "Human Intestine data",
    Dataset == "crohn_caic"   ~ "Crohn’s Disease data",
    Dataset == "soil_caic"    ~ "Soil data"
    #TRUE ~ Dataset  # default case, in case other datasets exist
  ))
other_aicc_dd$Name  = "aicc"

aic_df  = rbind(caic_dd,other_aicc_dd) %>% arrange(Dataset)

aic =  aic_df   %>%
       group_by(Dataset) %>%
       mutate(Delta_AIC = Estimate - min(Estimate)) %>%
       ungroup()
#######################################################
pp = ggplot(aic, aes(Model,Delta_AIC)) +
  geom_point() +
  custom_theme(18) +
  facet_wrap(~Dataset, scale="free") +
  labs(x = " ", y = "Change in AIC")

size = 3; nn= 13; width =  13; height =  7; dpi = 300

ggsave("fig/aiccc_rr.png", plot = pp, 
       width = width, 
       height = height, 
       dpi = dpi)

 #######################################################
res     =  lapply(mod_list, function(x){
                           pp  =  list(RR = x$RR, RRzi = x$RRzi,
                                       US = x$US, USzi =x$USzi)
                           
                           dd   <-  data.frame(aictab(pp, modnames = names(pp)))
                           df1  <-  dd %>% select(Modnames,AICc)
                           na_indices <- which(is.na(df1$AICc))
                           
                           if (length(na_indices) > 0) {
                             df1$AICc[na_indices] <- sapply(na_indices, function(i) my_aicc_fun(pp[[df1$Modnames[i]]]))
                           }
                           df2  <- data.frame(Modnames =  c("NB", "ZNB", "DE"),
                                              AICc     =  c(NB   =  x$NB,
                                                            ZNB =  x$ZNB,
                                                              DE   =   x$DE))
                           df   =   rbind(df1, df2)
                           df = bind_rows(df1, df2) %>%
                             mutate(Delta_AICc = AICc - min(AICc)) %>%
                             arrange(Delta_AICc)
                           rownames(df) =  df$Modnames
                           df
                           })

res
#############################################################

aic_filtered <- aic %>%
  #filter(Model %in% c("US", "USzi")) %>%
  select(Dataset, Model, Estimate) %>%
  setNames(c("Dataset", "Model", "AICc")) %>%
  mutate(Source = "aic_conditional")

# Convert res list to a structured format
res_list <- list(
  autism = res$autism_models,
  atlass = res$atlass_models,
  crohn  = res$crohn_models,
  soil   = res$soil_models
)

res_filtered <- bind_rows(lapply(names(res_list), function(dataset) {
  df <- as.data.frame(res_list[[dataset]])  
  df <- tibble::rownames_to_column(df, "Model")  # Extract row names as a column
  df$Dataset <- paste0(dataset, "_caic")  # Match dataset naming
  df$Source <- "res"
  df
}), .id = "ID") 

rrr = res_filtered %>%
  #filter(Model %in% c("US", "USzi")) %>%
  select(Model, AICc,Dataset)  %>%
  mutate(Source = "aicc_marginal")

dddd = rbind(rrr,aic_filtered)

ggplot(dddd, aes(x = Model, y = AICc, group = Model, color = Source)) +
  #geom_line(aes(linetype = Model), size = 1) +
  geom_point(size = 3) +
  facet_wrap(~Dataset, scales = "free_y") +
  theme_bw() +
  labs(title = "Comparison of AIC for US and USzi Models",
       x = "Source",
       y = "AIC Estimate",
       color = "type")

#' the 4th data set fitting zi for each taxa leads to worse AICc
#' Nothing changes for the 3rd dataset in terms of the order but each zi is better
#' 
#############################################################
plt   =   list()
nam   =    c("Autism data",
             "Atlass data", 
             "Crohn disease data", 
             "Soil data")

for(i in 1:length(res)){
  ddf      =      res[[i]]
  ddf$Delta_AICc[ddf$Delta_AICc == 0] <- 1e-6
  plt[[i]] = ggplot(ddf, aes(Modnames,Delta_AICc)) + 
            geom_point() +
             scale_y_log10() +
              custom_theme(12) +
             labs(x ="", y = "Change in AICc", title = nam[i])
}

(plt[[1]]|plt[[2]])/(plt[[3]]|plt[[4]])
#############################################################
############Confidence Intervals
CI_filenames <- c(
  "CI_rr.rds",
  "CI_rrzi.rds",
  "CI_us.rds",
  "CI_uszi.rds",
  "CI_nbmm.rds",
  "CI_zinbmm.rds",
  "CI_deseq.rds")


CI_filenames1 <- c(
  "CI_rr.rds",
  "CI_rrzi_each.rds",
  "CI_us.rds",
  "CI_uszi.rds",
  "CI_nbmm.rds",
  "CI_zinbmm.rds",
  "CI_deseq.rds"
)

CI_filenames11 <- c(
  "CI_rr.rds",
  "CI_rrzi_wald3.rds",
  "CI_us.rds",
  "CI_uszi.rds",
  "CI_nbmm.rds",
  "CI_zinbmm.rds",
  "CI_deseq.rds"
)
##############################################################
autism_confint <-  load_models(autism_path, CI_filenames)
atlass_confint <-  load_models(atlass_path, CI_filenames11)
crohn_confint  <-  load_models(crohn_path, CI_filenames)
soil_confint   <-  load_models(soil_path, CI_filenames11)

CI_filenames11 <- c("CI_nbmm.rds",
                    "CI_zinbmm.rds",
                    "CI_deseq.rds",
                    "CI_us_wald3.rds",
                    "CI_uszi_wald3.rds",
                    "CI_rr_wald3.rds",
                    "CI_rrzi_wald3.rds")


autism_confint <-  load_models(autism_path, CI_filenames11)
atlass_confint <-  load_models(atlass_path, CI_filenames11)
crohn_confint  <-  load_models(crohn_path, CI_filenames11)
soil_confint   <-  load_models(soil_path, CI_filenames11)

View(atlass_confint)
# Assigning names 
CI_names <- c("RR","RRzi","US","USzi","NB","ZINB","DE")
names(autism_confint)    =   names(atlass_confint)  =    CI_names
names(crohn_confint)    =   names(soil_confint)    =    CI_names
CI_list    =   lst(autism_confint,atlass_confint,crohn_confint,soil_confint)


kk=lapply(CI_list, function(x){ 
  lapply(x, function(y){
   y$padjust    =  p.adjust(y$pvalue, method = "BH")
   y$p_reject   =  if_else(y$padjust < 0.05, 1, 0)
   y$n       =  length(y$p_reject)
   y$x       =   sum(y$p_reject)
   y$power   =   mean(y$p_reject)
   y
   pk =  apply(y, 1, function(x){
     binom.test(x = as.numeric(x[["x"]]), n= as.numeric(x[["n"]]), 
                p =  as.numeric(x[["power"]]))$conf.int
   })
   pq  = cbind(y,t(pk))
   colnames(pq) =  c(colnames(y),c("pow_lwr","pow_upr"))
   pq
  })
})

library(dplyr)
library(purrr)

pp_df <- imap_dfr(kk, function(model_list, dataset_name) {
  imap_dfr(model_list, function(y, model_name) {
    y %>%
      summarise(
        power = unique(power),
        lwr = unique(pow_lwr),
        upr = unique(pow_upr)
      ) %>%
      mutate(model = model_name, dataset = dataset_name) %>%
      mutate(dataset = recode(dataset,
                              "autism_confint" = "Autism data",
                              "atlass_confint" = "Human intestine data",
                              "crohn_confint" = "Crohn disease data",
                              "soil_confint" = "Soil data"))
    
  })
})


library(ggplot2)

gp_s = ggplot(pp_df, aes(x = model, y = power, ymin = lwr, ymax = upr)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
  labs(
    #title = "Power Estimates with 95% Confidence Intervals",
    x = " ",
    y = "Estimated average power across taxa",
    color = "Dataset"
  ) +
  custom_theme(18) +
  facet_wrap(~dataset, scale = "free")

ggsave(paste0("fig/power_realdata.png"), plot = gp_s, width = 12, height = 6, dpi = 300)


# ppq <- lapply(kk, function(z){
#  kk= lapply(z, function(y){
#     pk =  apply(y, 1, function(x){
#          binom.test(x = as.numeric(x[["x"]]), n= as.numeric(x[["n"]]), 
#                  p =  as.numeric(x[["power"]]))$conf.int
#     })
#     t(pk)
#   })
# })
names(kk)
class(kk$autism_confint)
#View(ppq$autism_confint$RR)
#################################################################
combine <- lst(
  autism_combined =  bind_rows(autism_confint, .id = "type"),
  atlass_combined =  bind_rows(atlass_confint, .id = "type"),
  crohn_combined  =  bind_rows(crohn_confint, .id = "type"),
  soil_combined   =  bind_rows(soil_confint, .id = "type")
)

combine_filt  =  combine
combine_filt   <- lapply(combine, filter_complete_separation)
names(combine_filt) <- c("Autism data",
                         "Human intestine data",
                         "Crohn disease data",
                         "Soil data") #names(combine)

View(combine_filt$`Soil data`)
################################################################
pplt  =  list()
for(i in 1:length(combine_filt)){
  ddff      =   combine_filt[[i]]
  plt_tile  =   names(combine_filt)[i]
  # dd1 <- ddff %>%
  #   group_by(type) %>%
  #   summarise(mean_width = range(pvalue, na.rm = TRUE))

  dd <- ddff %>%
    group_by(type) %>%
    summarise(mean_width = mean(width, na.rm = TRUE),
              sd_err     =    sd(width)/length(width),
              lwr        =    mean_width -  1.96*sd_err,
              upr        =    mean_width + 1.96*sd_err)
  
  pplt[[i]] <-     ggplot(dd, aes(x=type, y = mean_width,ymin = lwr, ymax = upr)) +
    geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
              geom_point(size =1) +
              #geom_line() +
              custom_theme(12) +
              ggtitle(plt_tile) +
              #theme(axis.ticks.x = element_blank(),   
             #       axis.text.x = element_blank()) +
              labs(x = " ", y = "Average confidence width", color = "model") +
    custom_theme(18)
}

library(patchwork)
pr = (pplt[[1]]|pplt[[3]])/(pplt[[2]]|pplt[[4]]) +   plot_layout(guides = "collect") 
ggsave(paste0("fig/confwidth_realdata.png"), plot = pr, width = 12, 
       height = 9, dpi = 300)

#################################################################
meant_count_file  <-  "mean_count.rds"
autism_mean_count <-  load_models(autism_path, meant_count_file)
atlass_mean_count <-  load_models(atlass_path, meant_count_file)
crohn_mean_count  <-  load_models(crohn_path,  meant_count_file)
soil_mean_count   <-  load_models(soil_path, meant_count_file)
#################################################################

#mod =autism_models$US
pppr = lapply(atlass_confint[names(atlass_confint) != "DE"], function(x){ 
  pval  =   as.numeric(x$pvalue)
  print(range(pval))
  x$p_adj <- p.adjust(pval, method = "bonferroni")
  x} )       


pqp=autism_confint[[1]]$pvalue
p.adjust(pqp, method = "BH")


ntaxa_vec    <-    c(autism_mod  =  length(autism_nbmm_mod[[1]]$responses),
                   atlass_mod  =   length(atlass_nbmm_mod[[1]]$responses),
                   crohn_mod  =  length(crohn_nbmm_mod[[1]]$responses),
                   soil_mod    =  length(soil_nbmm_mod[[1]]$responses))




mean_count_list <- list(
                  autism_mod  =  readRDS(paste0(autism_path,"mean_count.rds")),
                  atlass_mod  =   readRDS(paste0(atlass_path,"mean_count.rds")),
                  crohn_mod  =   readRDS(paste0(crohn_path,"mean_count.rds")),
                  soil_mod    =   readRDS(paste0(soil_path,"mean_count.rds"))
                  )

View(dd_full)
dd_full   =   wald_confint(mod = crohn_models$RR)

grp_dd    <-  extract_grp_effect(dd_full,ntaxa)

ggplot(grp_dd,aes(param_name, width, group=1)) +
  geom_point() +
  geom_line() +
  custom_theme(12)


# Add model type
dd$model =   rep("wald",nrow(dd))


library(zCompositions)
library(compositions)
library(corrplot)
library(reshape2)
library(ggplot2)
library(dplyr)
#############################################################
# Step 1: Create Synthetic Microbiome Dataset
set.seed(123)
# Parameters
n_sites <- 5
n_groups <- 2
n_taxa <- 10
n_samples_per_site <- 20

# Create site, group, and sample IDs
site <- rep(paste0("Site", 1:n_sites), each = n_samples_per_site)
group <- rep(rep(c("Control", "Treatment"), each = n_samples_per_site / 2), times = n_sites)
sample_id <- paste0("Sample", 1:(n_sites * n_samples_per_site))

# Generate taxon counts with site-specific patterns
taxa_counts <- matrix(0, nrow = length(site), ncol = n_taxa)
colnames(taxa_counts) <- paste0("Taxon", 1:n_taxa)

for (taxon in 1:n_taxa) {
  for (s in 1:n_sites) {
    site_samples <- which(site == paste0("Site", s))
    mean_count <- 5 + s * 2 + rnorm(1, 0, 2)
    taxa_counts[site_samples, taxon] <- rpois(length(site_samples), lambda = mean_count)
  }
}

# Introduce group effects (e.g., treatment increases certain taxa counts)
group_effect_taxa <- c(1, 2)  # Group effects on Taxon1 and Taxon2
for (taxon in group_effect_taxa) {
  taxa_counts[group == "Treatment", taxon] <- taxa_counts[group == "Treatment", taxon] +
    rpois(sum(group == "Treatment"), lambda = 5)
}

# Create a data frame
microbiome_data <- as.data.frame(taxa_counts)
microbiome_data$site <- site
microbiome_data$group <- group
microbiome_data$sample_id <- sample_id

# Step 2: Apply CLR Transformation to Address Compositionality
# Replace zeros using Bayesian multiplicative replacement
taxa_counts_no_zeros <- cmultRepl(taxa_counts, label = 0, method = "CZM")

# CLR transformation
clr_transformed <- apply(taxa_counts_no_zeros, 1, function(row) clr(as.numeric(row)))
clr_transformed <- t(clr_transformed)  # Transpose to match original structure

# Convert CLR-transformed data to a data frame
clr_data <- as.data.frame(clr_transformed)
colnames(clr_data) <- paste0("Taxon", 1:n_taxa)
clr_data$site <- microbiome_data$site
clr_data$group <- microbiome_data$group
clr_data$sample_id <- microbiome_data$sample_id

# Step 3: Explore Taxon Correlations Across Sites (CLR-Transformed Data)
correlation_by_site_clr <- lapply(unique(clr_data$site), function(s) {
  site_data <- clr_data[clr_data$site == s, 1:n_taxa]
  cor(site_data)
})
View(clr_data)
#############################################################
crohn_path <- paste0(getwd(), "/real_data/CrohnD_data/results/")

crohn_dir <- dirname(crohn_path)
source(paste0(crohn_dir,"/prep_data.R"))
#############################################################
# Load required packages
library(zCompositions)
library(compositions)
library(corrplot)
library(reshape2)
library(ggplot2)
library(dplyr)
library(glmmTMB)

# Ensure your data contains columns: count, age, group, taxon, subject, normalizer
microbiome_data        <-  as.matrix(t(pp$data$countdata))
microbiome_data        <-   cbind(microbiome_data, pp$data$meta_data)
stopifnot(rownames(microbiome_data) == pp$data$meta_data$subject)
rownames(microbiome_data)  <- NULL

# Step 2: CLR Transformation to Address Compositionality
# Replace zeros using Bayesian multiplicative replacement
taxa_counts_no_zeros <- cmultRepl((microbiome_data[, grep("taxon", colnames(microbiome_data))]), 
                                  label = 0, method = "CZM")

# CLR transformation
clr_transformed <- apply(taxa_counts_no_zeros, 1, function(row) clr(as.numeric(row)))
clr_transformed <- t(clr_transformed)  # Transpose to match original structure

# Convert CLR-transformed data to a data frame
clr_data <- as.data.frame(clr_transformed)
colnames(clr_data) <- grep("taxon", colnames(microbiome_data), value = TRUE)
clr_data$age <- microbiome_data$age
clr_data$group <- microbiome_data$group
clr_data$subject <- microbiome_data$subject
clr_data$normalizer <- microbiome_data$normalizer

subject_data = cor(clr_data[,grep("taxon", colnames(clr_data))])

# Step 5: Visualize Correlations Across Subjects (CLR-Transformed Data)
correlation_by_subject <- lapply(unique(clr_data$subject), function(s) {
  subject_data <- clr_data[clr_data$subject == s, grep("taxon", colnames(clr_data))]
  cor(subject_data)
})



# Visualize correlation matrices
par(mfrow = c(2, 3))
for (i in seq_along(correlation_by_subject)) {
  corrplot(correlation_by_subject[[i]], method = "color", 
           title = paste0("CLR-Subject ", i), mar = c(0, 0, 2, 0))
}

# Step 6: Check for Group-Level Effects (CLR-Transformed Data)
# Melt CLR-transformed data for group-level comparison
clr_melted_data <- melt(clr_data, id.vars = c("age", "group", "subject"), 
                        variable.name = "taxon", value.name = "clr_count")

# Calculate mean CLR-transformed abundance by group
group_means_clr <- clr_melted_data %>%
  group_by(group, taxon) %>%
  summarise(mean_clr = mean(clr_count), .groups = "drop")

# Visualize group differences for each taxon
ggplot(group_means_clr, aes(x = taxon, y = mean_clr, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Mean CLR-Transformed Taxon Abundance by Group", 
       x = "Taxon", y = "Mean CLR Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


















# Step 2: Apply CLR Transformation to Address Compositionality
# Replace zeros using Bayesian multiplicative replacement
dd1     =     dd
rownames(dd1) <- NULL

taxa_count1 =  as.matrix(dd1)
taxa_counts_no_zeros <- cmultRepl(taxa_count1, label = 0, method = "CZM")
#############################################################
# CLR transformation
clr_transformed <- apply(taxa_counts_no_zeros, 1, function(row) clr(as.numeric(row)))
clr_transformed <- t(clr_transformed)  # Transpose to match original structure

# Convert CLR-transformed data to a data frame
clr_data <- as.data.frame(clr_transformed)
colnames(clr_data) <- paste0("Taxon", 1:n_taxa)
clr_data$site <- microbiome_data$site
clr_data$group <- microbiome_data$group
clr_data$sample_id <- microbiome_data$sample_id
#############################################################

# Step 3: Explore Taxon Correlations Across Sites (CLR-Transformed Data)
correlation_by_site_clr <- lapply(unique(clr_data$site), function(s) {
  site_data <- clr_data[clr_data$site == s, 1:n_taxa]
  cor(site_data)
})

# Visualize correlation matrices
par(mfrow = c(2, 3))
for (i in seq_along(correlation_by_site_clr)) {
  corrplot(correlation_by_site_clr[[i]], method = "color", 
           title = paste0("CLR-Site ", i), mar = c(0, 0, 2, 0))
}

# Step 4: Check for Group-Level Effects (CLR-Transformed Data)
# Melt CLR-transformed data for group-level comparison
clr_melted_data <- melt(clr_data, id.vars = c("site", "group", "sample_id"), 
                        variable.name = "taxon", value.name = "clr_count")

# Calculate mean CLR-transformed abundance by group
group_means_clr <- clr_melted_data %>%
  group_by(group, taxon) %>%
  summarise(mean_clr = mean(clr_count), .groups = "drop")

# Visualize group differences for each taxon
ggplot(group_means_clr, aes(x = taxon, y = mean_clr, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Mean CLR-Transformed Taxon Abundance by Group", 
       x = "Taxon", y = "Mean CLR Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
