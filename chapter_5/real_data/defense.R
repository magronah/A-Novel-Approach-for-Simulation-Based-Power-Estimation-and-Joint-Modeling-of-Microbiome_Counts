library(glmmTMB)
library(here)
library(matrixStats)
library(ggplot2)
library(patchwork)
###########################################################
path0   =   paste0("real_data/")
source(paste0(path0,"/fun.R"))
plt = "defense/fig/"
###########################################################
path1   =   paste0("real_data/CrohnD_data/results/")
path2   =   paste0("real_data/autism_data/results/")
path3   =   paste0("real_data/atlass_data/results/")
path4   =   paste0("real_data/soil_data/results/")
##############################################################
countdata  =   list(autism =  readRDS(paste0(path2,"countdata.rds")),
                    crohn  =  readRDS(paste0(path1,"countdata.rds")),
                    `human intestine`  = readRDS(paste0(path3,"countdata.rds")),
                    soil  = readRDS(paste0(path4,"countdata.rds")))
##############################################################
pp5 = lapply(countdata, dim)
combined_df5 <- do.call(rbind, lapply(seq_along(pp5), function(i) {
  data.frame(x = pp5[[i]][2], group = names(pp5)[i])
}))

p5 = ggplot(combined_df5, aes(x = group, y = (x))) +
  geom_point() +
  labs(x = "Group", y = "Number of samples") +
  custom_theme(n)

combined_df6 <- do.call(rbind, lapply(seq_along(pp5), function(i) {
  data.frame(x = pp5[[i]][1], group = names(pp5)[i])
}))

n = 20; width = 17;  height = 7; dpi = 300; size=3
p51 = ggplot(combined_df5, aes(x = group, y = (x))) +
  geom_point(size=size) +
  labs(x = "dataset", y = "Number of samples") +
  custom_theme(n)

p52 = ggplot(combined_df6, aes(x = group, y = (x))) +
  geom_point(size=size) +
  labs(x = "dataset", y = "Number of taxa") +
  custom_theme(n)

p5 =  p51|p52
ggsave(paste0(plt,"real_data_ss_ntaxa.png"), 
       plot = p5, width = width, height = height, dpi = dpi)
##############################################################
pp1 = lapply(countdata, function(x){(data.frame(x=rowMeans(x)))})
combined_df <- do.call(rbind, lapply(seq_along(pp1), function(i) {
  data.frame(x = pp1[[i]]$x, group = names(pp1)[i])
}))

p1 = ggplot(combined_df, aes(x = group, y = log(x))) +
  geom_boxplot() +
  labs(x = "dataset", y = "log(mean of taxa)") +
  custom_theme(n)
###################################################  
pp2 = lapply(countdata, function(x){(data.frame(x=sqrt(rowVars(as.matrix(x)))))})
combined_df2 <- do.call(rbind, lapply(seq_along(pp2), function(i) {
  data.frame(x = pp2[[i]]$x, group = names(pp2)[i])
}))

p2 = ggplot(combined_df2, aes(x = group, y = log(x))) +
  geom_boxplot() +
  labs(x = "dataset", y = "log(standard deviation") +
  custom_theme(n)

p12 = p1|p2
ggsave(paste0(plt,"real_data_mean_sd.png"), 
       plot = p12, width = width, height = height, dpi = dpi)
###################################################  
pp3 = lapply(countdata, function(x){(data.frame(x=colVars(as.matrix(x))))})
combined_df3 <- do.call(rbind, lapply(seq_along(pp3), function(i) {
  data.frame(x = pp3[[i]]$x, group = names(pp3)[i])
}))

ggplot(combined_df3, aes(x = group, y = log(x))) +
  geom_boxplot() +
  labs(x = "Group", y = "Row Mean") +
  custom_theme(n)
###################################################  
pp4 <- lapply(countdata, function(x) {
  mat <- as.matrix(x)
  prop_zeros <- colSums(mat == 0) / nrow(mat)
  data.frame(x = prop_zeros)
})

combined_df4 <- do.call(rbind, lapply(seq_along(pp4), function(i) {
  data.frame(x = pp4[[i]]$x, group = names(pp4)[i])
}))

ggplot(combined_df4, aes(x = group, y = (x))) +
  geom_boxplot() +
  labs(x = "Group", y = "Row Mean") +
  custom_theme(n)
###################################################  


