setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
library(NBZIMM)

if(FALSE){
path = paste0(getwd(),"/reproducible/power/datasets/")

countdata    =  readRDS(paste0(path,"otu_dataset_list.rds"))
metadata     =  readRDS(paste0(path,"metadata.rds"))

cc      =   commandArgs(trailingOnly  = TRUE)
i       =   as.integer(cc[1])

name1   =   names(countdata$filtered_otu)[i]
name2   =   names(metadata)[i]
stopifnot(name1 == name2)

count_dd  =   as.data.frame(t(countdata$filtered_otu[[i]]))
met_dd    =   metadata[[i]]


met_dd$dummy = factor(1)
mod        =   mms(y = count_dd, fixed = ~Groups,
                     random = ~ 1|dummy,
                     zi_fixed = ~1,
                     data = met_dd, method = "zinb") 

saveRDS(mod, file=paste0("~/scratch/dataset/new_sim/zero_inflation_est/mod",i, ".rds"))



l = list()
for(k in 1:7){
r = readRDS(paste0("~/scratch/dataset/new_sim/zero_inflation_est/mod",k, ".rds"))
l[[k]] = range(as.numeric(unlist(lapply(r$fit,function(x){plogis(x$zi.fit[[1]])}))))
}

l

saveRDS(l, file=paste0(getwd(),"/reproducible/new_sim/zero_inflation_est.rds"))
}


# Load necessary libraries
library(ggplot2)

zero_inflation_ranges <- readRDS(paste0(getwd(),"/reproducible/new_sim/zero_inflation_est.rds"))
#list(
#  c(0.2254422, 0.8787500),
#  c(0.09176632, 0.89722222),
#  c(0.001748252, 0.940034965),
#  c(0.04507106, 0.91902174),
#  c(0.4695353, 0.9268293),
#  c(0.003937008, 0.927559055),
#  c(0.006024096, 0.915662651)
#)

# Convert the list to a data frame with dataset number, lower, and upper bounds
zero_inflation_df <- data.frame(
  Dataset = 1:7,  # Dataset numbers
  Lower = sapply(zero_inflation_ranges, function(x) x[1]),
  Upper = sapply(zero_inflation_ranges, function(x) x[2])
)

# Add midpoint values for better visualization
zero_inflation_df$Midpoint <- (zero_inflation_df$Lower + zero_inflation_df$Upper) / 2


# Create the plot using ggplot2
ggpl1 = ggplot(zero_inflation_df, aes(x = factor(Dataset), y = Midpoint)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper), color = "blue") +
  labs(
    title = "Zero-Inflation Probability Ranges for 7 Datasets",
    x = "Dataset",
    y = "Zero-Inflation Probability"
  ) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits from 0 to 1
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) 

