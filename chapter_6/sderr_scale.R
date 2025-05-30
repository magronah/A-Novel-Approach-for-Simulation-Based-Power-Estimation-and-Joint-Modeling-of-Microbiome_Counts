library(glmmTMB)
library(ggplot2)
source("func.R")
############################################
r      =   list()
files  =  list.files("~/scratch/long/new_sim/30_100_3/us/")
for(i in 1:length(files)){
  mod =  readRDS(paste0("~/scratch/long/new_sim/30_100_3/us/mod",i,".rds"))
  ref    =  ranef(mod, condVar = FALSE)
  n1     =  length(names(ref$cond$taxon))
  n2     =  sum(grepl(".*:time", names(ref$cond$taxon)))
  dd     =  model.frame(mod)
  ntaxa  =  length(unique(dd$taxon))
  ############################################
  ss <- TMB::sdreport(mod$obj, getJointPrecision = TRUE)
  ss$jointPrecision <-  as(ss$jointPrecision, "sparseMatrix")
  inverse_mat       <-  solve(ss$jointPrecision)
  ########################################################
  se_vec <- sqrt(diag(inverse_mat))
  se_vec_b  <-  se_vec[names(se_vec)=="b"]
  ########################################################
  jp <- TMB::sdreport(mod$obj, getJointPrecision = TRUE)$jointPrecision
  ########################################################
  b_inds <- which(rownames(jp) == "b")
  jpb0   <- jp[b_inds, b_inds]
  ########################################################
  n1     <-  length(names(ref$cond$taxon))
  b_sub_inds <-  1:(n1*ntaxa)
  jpb        <-  jpb0[b_sub_inds, b_sub_inds]
  inv_mat    <-  solve(jpb)
  se_vec_b2    <- sqrt(diag(inv_mat))
  #########################################################
  grp_label <-   grep("[:].*(time|days|day|week|weeks|times)", names(ref$cond$taxon), 
                      ignore.case = TRUE, 
                      value = TRUE)
  n3   <-    length(grp_label)
  #########################################################
  p     <-  1:(n1*ntaxa)
  nselect  <- split(p, ceiling(seq_along(p) / n1))  
  subselect <- lapply(nselect, function(x) tail(x, n3))  
  grp_ind   <- unique(unlist(subselect))
  ########################################################
  se_vec1   <-   se_vec_b[grp_ind]
  se_vec2   <-   se_vec_b2[grp_ind]
  r[[i]]         <-   se_vec1/se_vec2
}

mean(unlist(lapply(r, median)))

dd  = data.frame(ratio = unlist(r), 
                 name = rep(paste0("sim", 1:length(r)), each = ntaxa))

dd$name   =  factor(dd$name, levels = unique(dd$name))
plt =  ggplot(dd, aes(x = name, y=ratio)) +
      geom_boxplot() +
  labs(
    x = "Simulation Dataset",
    y = "Ratio of Standard Errors",
    title = "Distribution of Standard Error Ratios Across Simulations"
  ) +
      custom_theme(18) 

width =  8; height =  6; dpi = 300
ggsave("fig/ratio.png", plot = plt, 
       width = width, 
       height = height, 
       dpi = dpi)

