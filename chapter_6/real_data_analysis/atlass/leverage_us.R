setwd("/home/agronahm/projects/def-bolker/agronahm/longitudinal_RR/")
library(glmmTMB)
library(ggplot2)
library(dplyr)
library(patchwork)
#library(here)
##############################################################
nnp = c("atlass", "pregnancy")
for(i in 1:2){
  path <- paste0("real_data_analysis/", nnp[i],"/results/")
  print(path)
# path <- paste0("real_data_analysis/atlass/results/")

##############################################################
fig_path = "fig/"
source("func.R")
##############################################################
mod  <-  readRDS(paste0(path, "mod_us.rds"))

par_ctrl <- glmmTMBControl(parallel = list(n = 10, autopar = TRUE))

gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")
############################################################
cc     =   commandArgs(trailingOnly  = TRUE)
i      =   as.integer(cc[1])
############################################################
mmd  <-  mod
df   <-  model.frame(mod) 
dd   <-  df |> dplyr::rename(normalizer = "offset(normalizer)")

lev  <-  leverage_brute_modified(mmd, data = dd, inds = i, eps = 0.1, 
                        fit_method = "update", scale = "response",
                        pred_method  = "predict",
                        progress = TRUE)


saveRDS(lev, )
}

file_path  =  paste0(path,"leverage/us/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}

saveRDS(lev, file= paste0(file_path,"lev",i, ".rds"))


if(FALSE){

pp[1,1]; lev; lev1

eps_vec2 <- 10^seq(-4, 0, length = 31)
lev_vals_nb <- sapply(eps_vec2, \(e) 
                      leverage_brute_modified(mmd, data = dd, 
                                              inds = 1, eps = e, 
                                              fit_method = "update", scale = "response",
                                              pred_method = "predict",
                                              progress = TRUE))
plot(eps_vec2, lev_vals_nb, log = "x")

eps_vec3 <- seq(0, 4, length = 10)
lev_vals_nb3 <- sapply(eps_vec3, \(e) 
                      leverage_brute_modified(mmd, data = dd, 
                                              inds = 1, eps = e, 
                                              fit_method = "update", scale = "response",
                                              pred_method = "predict",
                                              progress = TRUE))
pp2= c(eps_vec2,eps_vec3)

pp1= c(lev_vals_nb,lev_vals_nb3)
plot(pp2, pp1, log = "x")

model = mmd
inds = 1
evals <- 10^(-3:0)
lev_vals <- lapply(evals,
                   \(e) leverage_brute(mmd, eps = e, data = dd,inds = i,progress = TRUE, scale = "response")
)


y_hat <- function(model, data) {
  predict(model, newdata = data, type = "response", re.form = NULL)
}

data  = dd
model =  mod

library(Matrix)
leverage_brute_force1 <- function(model, data,epsilon = 5,  nn = 5) {
  n <- nrow(data)
  H <- Matrix(0, n, nn, sparse = TRUE)  
  y_pred <- y_hat(model, data)  
  yname  <- as.character(formula(model)[[2]])
  
  p0 <- with(model$obj$env, parList(last.par.best[-random]))
  p0 <- p0[lengths(p0) > 0]
  p0 <- p0[setdiff(names(p0), "b")]
  
  for (i in 1:nn) {
    #i =1; epsilon = 0.1
    data_perturb <- data
    data_perturb[[yname]][i] <-  data[[yname]][i] + epsilon
    model_up <- update(model, start = p0, data = data_perturb)
    y_pred_up <- y_hat(model_up, data_perturb)
    
    H[, i] <- (y_pred_up - y_pred) / epsilon
  }
  return(H)
}

pp = leverage_brute_force1(model, data,epsilon = 0.1,  nn = 5) 
diag(pp)
lev
############################################################
file_path  =  paste0(path,"leverage/us/")

if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
  cat("Folder created at:", file_path, "\n")
} else {
  cat("Folder already exists at:", file_path, "\n")
}

saveRDS(lev, file= paste0(file_path,"lev",i, ".rds"))

}
