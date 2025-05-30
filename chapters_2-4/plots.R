library(ggrepel)
source("reproducible/power/Load_Packages.R")
source("reproducible/power/simulate_fun.R")
source("reproducible/power/utils.R")
source("reproducible/power/fitting_fun.R")
source("reproducible/power/read_dataset.R")
#path = "reproducible/power/datasets/sample_size_cal/"
########################Test####################################################
path  =  "reproducible/power/datasets2/"
countdata_list_obs  =  readRDS(file = paste0(path,"data.rds"))
metadata_list_obs   =  readRDS(file = paste0(path,"metadata.rds"))
################################################################################
#####lfc versus lm plot
c = c(2,5,3)
dd = foreach(i = c,.combine = "rbind") %do% {
  nam             =   names(logmean_list)[i]
  logmean         =   logmean_list[[i]]
  logfoldchange   =   logfoldchange_list[[i]]
         dd       =   data.frame(logmean, logfoldchange)
  dd$datatype     =   rep(nam, nrow(dd))
  dd
}

library(ggpmisc)
View(dd[dd$datatype == "PRJNA453621",])
unique(dd$datatype)
View(dd)

ggplot(dd, aes(logmean, logfoldchange)) +
  geom_point(alpha=0.1) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE)+
  #geom_smooth() +
  xlab(TeX("$\\log_2$(mean abundance)")) +
  ylab(TeX("$\\log_2$(fold change)")) +
  facet_wrap(~datatype, scale="free") +
  custom_theme(14)

################################################################################
i          =   2
meta       =   metadata_list_obs[[i]]
count      =   filtered_otu_list[[i]]
meta       =   metadata_list_obs[[i]]
deseq      =   deseqfun(count,meta,ref_name="NT")

res             =   deseq$deseq_estimate
logmean         =   deseq$logmean
logfoldchange   =   deseq$logfoldchange

dd     =    data.frame(logmean, logfoldchange)
dd1    =   subset(dd, abs(logfoldchange) > 2)
dd2     =   subset(dd1, logmean < 7 &  logmean > 1.5)

pp_rr0=ggplot(dd, aes(logmean, logfoldchange)) +
  geom_point(alpha=0.3) 
  geom_text_repel(aes(label = rownames(dd2)),
                  max.overlaps = nrow(dd2)) +
  xlab(TeX("$\\alpha$"))

nt     =    meta[meta$Groups == "NT", ]
asd    =    meta[meta$Groups == "ASD", ]

control    =   count[,colnames(count) %in% nt$Samples]
treatment  =   count[, colnames(count) %in% asd$Samples]

stopifnot(colnames(control) == nt$Samples)
stopifnot(colnames(treatment) == asd$Samples)
###################################################################
cont_zeros   =   control[which(rowMeans(control)==0),    ]
treat_zeros  =   treatment[which(rowMeans(treatment)==0), ]

dd_sub1      =   dd[rownames(dd) %in% rownames(cont_zeros), ]
stopifnot(rownames(dd_sub1) == rownames(cont_zeros))

dd_sub2      =   dd[rownames(dd) %in% rownames(treat_zeros), ]
stopifnot(rownames(dd_sub2) == rownames(treat_zeros))

pp_rr1  = ggplot(dd_sub1, aes(logmean, logfoldchange)) +
  geom_point(alpha=0.3) 
  geom_text_repel(aes(label = rownames(dd_sub1)),
                  max.overlaps = nrow(dd_sub1)) +
  xlab(TeX("$\\alpha$"))
  
pp_rr2 = ggplot(dd_sub2, aes(logmean, logfoldchange)) +
         geom_point(alpha=0.3) 
  geom_text_repel(aes(label = rownames(dd_sub2)),
                  max.overlaps = nrow(dd_sub2)) +
    xlab(TeX("$\\alpha$"))

ras0 <- rasterize(pp_rr0, layers='point', dpi=300)
ras1 <- rasterize(pp_rr1, layers='point', dpi=300)
ras2 <- rasterize(pp_rr2, layers='point', dpi=300)

ras1|ras2|ras0


ppp1   =   res[rownames(res) %in%  rownames(dd_sub1), ]
ppp2   =   res[rownames(res) %in%  rownames(dd_sub2), ]

v      =   (ppp1[is.na(ppp1$padj),])
vv1    =    dd[rownames(dd) %in%  rownames(v), ] 


ggplot(vv1, aes(logmean, logfoldchange)) +
  geom_point(alpha=0.3) 

v2      =   (ppp2[is.na(ppp2$padj),])
vv2    =    dd[rownames(dd) %in%  rownames(v2), ] 
ggplot(vv2, aes(logmean, logfoldchange)) +
  geom_point(alpha=0.3) 

control_sub    =  control[rownames(control) %in% inspect,]
treatment_sub  =  treatment[rownames(treatment) %in% inspect,]

ggplot(dd, aes(logmean, logfoldchange)) +
  geom_point(alpha=0.3)
  #geom_text_repel(aes(label = rownames(dd)),
  #                max.overlaps = nrow(dd)) +
 # xlab(TeX("$\\alpha$"))

control       =    control_otu_list[[i]]
treatment     =    treatment_otu_list[[i]]
cont_zeros    =    control[which(rowMeans(control)==0),    ]
treat_zeros   =    treatment[which(rowMeans(treatment)==0), ]

p1    =   rownames(cont_zeros)
p2    =   rownames(treat_zeros)

dd1   =   dd[rownames(dd) %in% p1,]

dim(dd)
length(p1)
stopifnot(rownames(dd1) == p1)
#nam = names(metadata_list_obs)


###############################################################################
## Plot relationship between statistical power, sample size and log fold change (log mean abundance = 5)
nsamples =  seq(30,190,20)
folder_path  =   paste0(path,"gam_model/")

##plot quantiles and average power
#' simulate data and predict 
nsamples =  seq(30,190,20)
notu = 1000; disp_scale = 0.3; nsim= 50; nsamp = 50
folder_path  =   paste0(path,"gam_model/")

otu_data_list  = list()
for(i in 1:7){
  dispersion_param      =   dispersion_param_list[[i]]
  logmean_param         =   logmean_param_list[[i]]
  logfoldchange_param   =   logfoldchange_param_list[[i]]
  otu_data_list[[i]]   =  countdata_sim_fun(logmean_param, logfoldchange_param, 
                                            dispersion_param,
                                            nsamp_per_group=nsamp,
                                            disp_scale = disp_scale,
                                            ncont = NULL,ntreat = NULL,
                                            notu, nsim, 
                                            maxlfc_iter = 10000,
                                            seed = 100)
}

##########################################################################

pp = foreach(i =  1:7, .combine = "rbind") %:%
  foreach(j =  1:length(nsamples),.combine = "rbind") %do%{
    
    nam        =   names(metadata_list_obs)[i]
    otu_data   =   otu_data_list[[i]]
    
    nsam       =   nsamples[j]
      mod      =   readRDS(paste0(folder_path,nsam,"_samples/",
                              nam, ".rds"))
  
     newdata   =   data.frame(lmean_abund  =  unlist(otu_data$logmean_list),
                              abs_lfc      =  unlist(otu_data$logfoldchange_list)) 

    power      =   predict(mod,newdata,type="response")
  quantiles    =   quantile(power, probs = c(0.5, 0.6, 0.7, 0.8))
  
    ddf        =   data.frame(sample_size = nsam, average_power = mean(power))
    ddf        =   cbind(ddf, as.data.frame(t(quantiles))) 
    
  # Rename the quantile columns
  colnames(ddf)[3:6] <- c("quantile_50%", "quantile_60%", "quantile_70%", "quantile_80%")
  ddf$name =   nam
  ddf
}

#oka_col = c("#556B2F", "#E23D28", "#0000FF","#000000")

oka_col = c(
  "#556B2F", 
  "#E23D28", 
  "#0000FF",
  "#E69F00", 
  "#000000"
)


sub = c("PRJNA168470","PRJNA355023","PRJNA644763","PRJNA589343")
pp_sub = pp[pp$name %in% sub,]

ggplot(pp_sub, aes(sample_size, average_power)) +
  geom_point(aes(shape = "average power across all taxa"), size = 3, color = "black") +  # Use shape for points
  geom_line(aes(y = `quantile_50%`, color = "50th Quantile", group = 1), lwd = 1.5) +
  geom_line(aes(y = `quantile_60%`, color = "60th Quantile", group = 1), lwd = 1.5) +
  geom_line(aes(y = `quantile_70%`, color = "70th Quantile", group = 1), lwd = 1.5) +
  geom_line(aes(y = `quantile_80%`, color = "80th Quantile", group = 1), lwd = 1.5) +
  scale_color_manual(name = "Quantile",  # Color legend for lines
                     values = c("50th Quantile" = oka_col[1], 
                                "60th Quantile" = oka_col[2], 
                                "70th Quantile" = oka_col[3], 
                                "80th Quantile" = oka_col[4])) +
  scale_shape_manual(name = "Average Power",  # Separate shape legend for points
                     values = c("average power across all taxa" = 16)) +  # Choose a shape for points
  custom_theme(16) + 
  labs(x = "Sample Size", y = "Average Power") +
  facet_wrap(~name, scales = "free") +
  guides(color = guide_legend(order = 1, reverse = TRUE),  # Reverse color legend
         shape = guide_legend(order = 2))  

###############################################################################
##' Expected hits
oka_col_cvd = c(
  "#000000", # Black
  "#E69F00", # Orange
  "#56B4E9", # Sky Blue
  "#009E73", # Bluish Green
  "#F0E442", # Yellow
  "#0072B2", # Blue
  "#D55E00"  # Vermilion
)

pp$expected  = notu * pp$average_power

ggplot(pp, aes(sample_size,expected, group=name, color=name)) +
  geom_point(size=3) +
  geom_line(lwd=1.2)  +
  scale_color_manual(name = "data accession no.", 
                     values= oka_col_cvd)+
  custom_theme(15) + labs(x="sample size", 
                          y="expected number of significant taxa") 
################################################################################
newdata    =   data.frame(lmean_abund  = c(5,5,5),
                          abs_lfc      =  c(2,3,4)) 
pp1 = foreach(i =  1:7, .combine = "rbind") %:%
  foreach(j =  1:length(nsamples),.combine = "rbind") %do%{
    
    nam        =   names(metadata_list_obs)[i]
    nsam       =   nsamples[j]
    mod        =   readRDS(paste0(folder_path,nsam,"_samples/",
                                nam, ".rds"))
    
    power      =   predict(mod,newdata,type="response")

    ddf        =   data.frame(sample_size = rep(nsam,nrow(newdata)), 
                              power = power)
    ddf        =   cbind(newdata,ddf)
  ddf$abs_lfc  =   factor(ddf$abs_lfc, levels = unique(ddf$abs_lfc))
    ddf$name   =   rep(nam,nrow(newdata))
    ddf
  }


sub = c("PRJNA168470","PRJNA687773","PRJNA644763","PRJNA453621")
pp1_sub = pp1[pp1$name %in% sub,]

oka_col_cvd = c(
  "#E69F00", # Orange
  "#56B4E9", # Sky Blue
  "#009E73" # Bluish Green
)

xlab(TeX("$\\log_2$(mean abundance)"))

ggplot(pp1_sub, aes(sample_size,power, group=abs_lfc, color=abs_lfc)) +
  geom_point(size=3) +
  geom_line(lwd=1.2)  +
 scale_color_manual(name = TeX("$\\log_{2}$(fold change)"), 
                     values= oka_col_cvd)+
  custom_theme(18) + labs(x="sample size", 
                          y="power") +
  facet_wrap(~name) +
  guides(color = guide_legend(reverse = TRUE), 
         linetype = guide_legend(reverse = TRUE))
names(metadata_list_obs)
################################################################################
nsamples =  seq(10,190,20)
folder_path  =   paste0(path,"contour_data/")
  
pp = foreach(i =  1:7, .combine = "rbind") %:%
    foreach(j =  1:length(nsamples),.combine = "rbind") %do%{
      
      nam      =   names(metadata_list_obs)[i]; nsam = nsamples[j]

      dd2      =   readRDS(paste0(folder_path,nsam,"_samples/power_estimate_",
                                  nam, "_.rds"))
      quantiles <- quantile(dd2$power, probs = c(0.5, 0.6, 0.7, 0.8))
      
      ddf = data.frame(sample_size = nsam,
                       average_power = mean(dd2$power))
      
      ddf <- cbind(ddf, as.data.frame(t(quantiles)))

      # Rename the quantile columns for clarity
      colnames(ddf)[3:6] <- c("Quantile_50%", "Quantile_60%", "Quantile_70%", "Quantile_80%")
  
      ddf$names =   nam
      ddf
    }

ggplot(pp, aes(sample_size,average_power)) +
  geom_point()+
  geom_line(aes(y = `Quantile_50%`, color = "50th Percentile", group = 1)) +
  geom_line(aes(y = `Quantile_60%`, color = "60th Percentile", group = 1)) +
  geom_line(aes(y = `Quantile_70%`, color = "70th Percentile", group = 1)) +
  geom_line(aes(y = `Quantile_80%`, color = "80th Percentile", group = 1)) +
  scale_color_manual(name = "Quantiles", 
                     values = c("50th Percentile" = "blue", 
                                "60th Percentile" = "black", 
                                "70th Percentile" = "green", 
                                "80th Percentile" = "red")) +
  custom_theme(10) + labs(x="sample size", y="overall power") +
  facet_wrap(~names)

#####################Check power distributions#####################
nsam = 110;    folder_path  =   paste0(path,"contour_data/")
  
  
  pp = foreach(i =  1:7, .combine = "rbind") %do%{
    nam = names(metadata_list_obs)[i] 
    dd2      =   readRDS(paste0(folder_path,nsam,"_samples/power_estimate_",
                                nam, "_.rds"))
    dd2$names =   rep(nam, nrow(dd2))
    
    dd2
  }
  
  ggplot(pp, aes(power))+
    geom_density() +
    facet_wrap(~names, scales = "free")
###############################################################################
nsamples = seq(5,100,5)
folder_path <- paste0(path,"gam_model/")
##plot quantiles and average power

pp = foreach(j =  1:length(nsamples)) %do%{
    nam  = names(metadata_list_obs)[1]; nsam = nsamples[j]
    dd2      =   readRDS(paste0(folder_path,nsam,"_samples/",
                                nam, ".rds"))
    dim(dd2$model)
    #print(dd2)
  }



dd = countdata_list_obs[[1]]
sample_info = metadata_list_obs[[1]]

# Run the DESeq2 differential expression analysis
deseq          =   deseqfun(dd,sample_info,ref_name="NT")
counts(dds)

# Get the results
results <- results(dds)
raw_counts <- counts(dds)
rownames(raw_counts) = paste0("taxon", 1:nrow(raw_counts))



library(DESeq2)

# Set seed for reproducibility
set.seed(123)

# Parameters for the simulation
num_genes <- 1000      # Number of genes (OTUs)
num_samples <- 20      # Number of samples
condition <- rep(c("Control", "Treatment"), each = num_samples / 2)

# Simulate count data
# We will simulate count data such that half of the genes are differentially expressed

baseline_counts <- rbinom(n=num_genes * num_samples, size = 500, prob=0.8)
counts_matrix <- matrix(baseline_counts, nrow = num_genes, ncol = num_samples)
num_zero_entries <- round(0.5 * length(counts_matrix))  # Set 10% of the entries to zero
zero_indices <- sample(length(counts_matrix), num_zero_entries)
counts_matrix[zero_indices] <- 0

# Introduce differential expression in half of the genes
fold_change_genes <- 1:500
fold_change <- 2  # Fold change of 2 for differentially expressed genes in treatment
counts_matrix[fold_change_genes, (num_samples/2 + 1):num_samples] <- 
  counts_matrix[fold_change_genes, (num_samples/2 + 1):num_samples] * fold_change

# Create a data frame for sample information
sample_info <- data.frame(
  row.names = paste0("Sample", 1:num_samples),
  condition = factor(condition)
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = sample_info,
  design = ~ condition
)

# Run the DESeq2 differential expression analysis
dds <- DESeq(dds,sfType ="poscounts")
#counts(dds)

# Get the results
results <- results(dds)
raw_counts <- counts(dds)
rownames(raw_counts) = paste0("taxon", 1:nrow(raw_counts))

#View(raw_counts)
# Inspect the top results
head(results)

# Plot MA plot to visualize the results
plotMA(results, main="DESeq2", ylim=c(-5,5))

# Print summary of results


summary(results)

# Access the estimated fold changes and p-values
fold_changes <- results$log2FoldChange
p_values <- results$pvalue

# Display the first few fold changes and p-values
head(fold_changes)
head(p_values)
dim(raw_counts)
length(fold_changes)
dd = data.frame(mean_count = rowMeans(raw_counts[,1:6]),fold_changes=fold_changes)
View(dd)

sub =  raw_counts[1:14,1:6]
dd1 = data.frame(variance = rowVars(sub))
View(sub)

library(ggplot2)

# Data points
y_values <- c(0.2, 0.1, 0.9, 0.8, 0.7, 0.8, 0.7)
x_labels <- paste0("data", 1:7)

# Create a dataframe
data <- data.frame(x = x_labels, y = y_values)

# Plot using ggplot2
ggplot(data, aes(x = x, y = y)) +
  geom_point(alpha = 8) +
  labs(x = "Dataset", y = "Maximum Power") +
  theme_minimal() +
  custom_theme(16)
  
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16)
  )

########################Test################################
logfoldchange = (logfoldchange_list[[3]])
logmean = logmean_list[[2]]

logfoldchange_param =   logfoldchange_param_list[[2]]
logmean_param       =   logmean_param_list[[2]]

notu =1000
logmean        =   logmean_sim_fun(logmean_param,notu)
logfoldchange  =   logfoldchange_sim_fun(logmean,logfoldchange_param, 
                                         max_lfc = 15,
                                         max_iter = 100)


compare_dataset <- function(countdata_sim_list,countdata_obs,method = c("var", "mean")){
  
  # Check if method is either "var" or "mean"
  if (!(method %in% c("var", "mean"))) {
    stop("Invalid method. Please choose either 'var' or 'mean'.")
  }
  
  # Calculate variance or mean based on the chosen method
  if (method == "var") {
    calc_func <- stats::var
    xlab_text <- "$\\log_{10}$(variance of taxa)"
  } else {
    calc_func <- mean
    xlab_text <- "$\\log_{10}$(mean of taxa)"
  }
  
  # Create a list combining simulation data and observed data
  dlist <- list.append(countdata_sim_list, obs_filt = countdata_obs)
  
  # Calculate variance or mean for each dataset in the list
  vars <- dlist |> purrr::map(~ apply(., 1, calc_func)) |>
    purrr::map_dfr(~ tibble(var = .), .id = "type")
  
  # Separate observed and simulated data
  vars_obs <- vars[vars$type == "obs_filt", ]
  vars_sim <- vars[vars$type != "obs_filt", ]
  
  # Create ggplot for visualization
  p <- ggplot(vars_sim, aes(x = log(var), colour = type)) + 
    geom_density(lwd = 1.5) +
    geom_density(data = vars_obs, aes(x = log(var)), 
                 colour = "black", linetype = "dashed", lwd = 1.5) +
    xlab(TeX(xlab_text))
  
  return(p)
}
