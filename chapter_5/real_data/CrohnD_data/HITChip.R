library(ANCOMBC)
library(phyloseq)
library(tidyverse)

# Set seed for reproducibility
set.seed(123)

# Create OTU table
otu_table <- as.data.frame(matrix(sample(1:1000, 50, replace = TRUE), nrow = 10))
colnames(otu_table) <- paste0("OTU", 1:5)
rownames(otu_table) <- paste0("Subject", 1:10)

# Create metadata
metadata <- data.frame(
  subject = rownames(otu_table),
  group = c(rep("A", 5), rep("B", 5)) # Two groups: A and B
)
rownames(metadata) <- rownames(otu_table)

# Convert to phyloseq object
OTU = otu_table(as.matrix(otu_table), taxa_are_rows = FALSE)
sample_data = sample_data(metadata)
physeq = phyloseq(OTU, sample_data)

# Fit the ANCOMBC2 model with random effects
set.seed(123)
out = ancombc2(data = otu_table, taxa_are_rows = FALSE,
               tax_level = NULL,
               meta_data = metadata,
               rand_formula = "(1 + group | taxon)",
               fix_formula = "group",
               p_adj_method = "holm")


# Extract and view results
res <- out$res
print(res)


library(ANCOMBC)
library(phyloseq)
library(tidyverse)
###########################################################
set.seed(123)
otu_table <- as.data.frame(matrix(sample(1:1000, 50, replace = TRUE), nrow = 10))
colnames(otu_table) <- paste0("OTU", 1:5)
rownames(otu_table) <- paste0("Subject", 1:10)

metadata <- data.frame(
  subject = rownames(otu_table),
  group = c(rep("A", 5), rep("B", 5)) # Two groups: A and B
)

rownames(metadata) <- rownames(otu_table)


OTU = otu_table(as.matrix(otu_table), taxa_are_rows = FALSE)
sample_data = sample_data(metadata)
physeq = phyloseq(OTU, sample_data)

class(otu_table)
set.seed(123)
out = ancombc2(data = otu_table, taxa_are_rows = FALSE,
               tax_level = NULL,
               meta_data = metadata,
               rand_formula = NULL,
               fix_formula = "group",
               p_adj_method = "holm")


names(out)
View(out$res)

prv_cut = 0.10,
lib_cut = 1000, main_var = "bmi", adj_formula = NULL, 
rand_formula = NULL, lme_control = NULL, struc_zero = TRUE,
neg_lb = TRUE, alpha = 0.05, n_cl = 2, verbose = TRUE)
res = out$res
################################################################
data(atlas1006, package = "microbiome")

# Subset to baseline
pseq = phyloseq::subset_samples(atlas1006, time == 0)

# Re-code the bmi group
meta_data = microbiome::meta(pseq)
meta_data$bmi = recode(meta_data$bmi_group,
                       obese = "obese",
                       severeobese = "obese",
                       morbidobese = "obese")

# Note that by default, levels of a categorical variable in R are sorted 
# alphabetically. In this case, the reference level for `bmi` will be 
# `lean`. To manually change the reference level, for instance, setting `obese`
# as the reference level, use:
meta_data$bmi = factor(meta_data$bmi, levels = c("obese", "overweight", "lean"))
# You can verify the change by checking:
# levels(meta_data$bmi)

# Create the region variable
meta_data$region = recode(as.character(meta_data$nationality),
                          Scandinavia = "NE", UKIE = "NE", SouthEurope = "SE", 
                          CentralEurope = "CE", EasternEurope = "EE",
                          .missing = "unknown")

phyloseq::sample_data(pseq) = meta_data

# Subset to lean, overweight, and obese subjects
pseq = phyloseq::subset_samples(pseq, bmi %in% c("lean", "overweight", "obese"))
# Discard "EE" as it contains only 1 subject
# Discard subjects with missing values of region
pseq = phyloseq::subset_samples(pseq, ! region %in% c("EE", "unknown"))

set.seed(123)
out = ancom(data = pseq, tax_level = "Family", meta_data = NULL,
            p_adj_method = "holm", prv_cut = 0.10,
            lib_cut = 1000, main_var = "bmi", adj_formula = "age + region", 
            rand_formula = NULL, lme_control = NULL, struc_zero = TRUE,
            neg_lb = TRUE, alpha = 0.05, n_cl = 2, verbose = TRUE)
res = out$res
#########################################

















# Run ANCOMBC
ancombc_results <- ancom(phyloseq = physeq,
                           formula = "group", # Use group as the variable of interest
                           p_adj_method = "holm", # Adjust p-values using Holm's method
                           zero_cut = 0.90, # Threshold for filtering low-abundance taxa
                           lib_cut = 1000, # Minimum library size to retain samples
                           group = "group", # Specify group variable
                           struc_zero = TRUE, # Structural zeros for zero-inflation
                           global = TRUE)

# Extract results
res <- ancombc_results$res

# Estimated fold changes
fc_table <- res$beta
fc_table

 

library(ANCOMBC)
library(tidyverse)
data(atlas1006, package = "microbiome")

x1         =   phyloseq(otu_table(atlas1006))
x2         =   phyloseq(tax_table(atlas1006))
meta_data  =   microbiome::meta(atlas1006)

meta_data$bmi = recode(meta_data$bmi_group,
                       obese = "obese",
                       severeobese = "obese",
                       morbidobese = "obese")


meta_data$bmi = factor(meta_data$bmi, levels = c("obese", "overweight", "lean"))


