source("reproducible/power/read_dataset.R")
source("reproducible/power/Load_Packages.R")
source("reproducible/power/simulate_fun.R")
source("reproducible/power/utils.R")
path = "reproducible/power/datasets/"
####################################################################
metadata_list     =    readRDS(paste0(path,"metadata.rds"))
data_list         =    readRDS(paste0(path,"data.rds"))
####################################################################
pp   =   compare_dataset(data_list,data_list[[1]],method = "mean")
pp   =   compare_dataset(data_list,data_list[[1]],method = "var")
####################################################################
otu_dataset_list     =    readRDS(paste0(path,"otu_dataset_list.rds"))
filtered_otu_list    =    otu_dataset_list[["filtered_otu"]]

count_data = t(filtered_otu_list[[1]])
lapply(filtered_otu_list, dim)

# Assume `dataset_list` is your named list of count matrices
plot_results <- lapply(names(filtered_otu_list), function(ds_name) {
  plot_diagnostics(filtered_otu_list[[ds_name]], ds_name)
})
names(plot_results) <- names(filtered_otu_list)

# Example: View plots for PRJNA168470
print(plot_results$PRJNA168470$p1)
print(plot_results$PRJNA168470$p2)
print(plot_results$PRJNA168470$p3)

# Or save to PDF
pdf("diagnostic_plots_PRJNA168470.pdf")
print(plot_results$PRJNA168470$p1)
print(plot_results$PRJNA168470$p2)
print(plot_results$PRJNA168470$p3)
dev.off()



prop_zero_taxa <- rowMeans(count_data == 0)
hist(prop_zero_taxa, main = "Proportion of Zero Counts per Taxon", xlab = "Proportion Zero", col = "skyblue")

prop_zero_samples <- colMeans(count_data == 0)
hist(prop_zero_samples, main = "Proportion of Zero Counts per Sample", xlab = "Proportion Zero", col = "salmon")


# Convert to data frame
otu_dims   =    lapply(filtered_otu_list, dim)
otu_df <- do.call(rbind, lapply(names(otu_dims), function(name) {
  data.frame(
    Dataset = name,
    Taxa = otu_dims[[name]][1],
    Samples = otu_dims[[name]][2]
  )
}))

ggplot(otu_df, aes(x = Samples, y = Taxa, label = Dataset)) +
  geom_point(size = 3, color = "steelblue") +
  geom_text(vjust = -0.5, hjust = 0.5) +
  labs(
    title = "Number of Samples vs Taxa in OTU Datasets",
    x = "Number of Samples",
    y = "Number of Taxa"
  ) +
  theme_bw()

lapply(filtered_otu_list, dim)

