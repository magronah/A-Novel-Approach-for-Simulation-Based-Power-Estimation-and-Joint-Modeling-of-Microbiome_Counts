library(reformulas)
library(glmmTMB)
library(tidyverse)
library(dplyr)

otu_meta_fun = function(sim_count, meta_dd){
  
  count_table   = list()
  met_dd     =    meta_dd
  nsim       =   length(sim_count)
  
  for(i in 1:nsim){
    meta_dd$count  =  sim_count[[i]]
    
    otu_table <- meta_dd %>%
      dplyr::select(subject, taxon, count, group)
    
    otu_table <- spread(otu_table, key = taxon, value = count) 
    colnames(otu_table) <- c("subject","group", 
                             paste0("taxon",1:(ncol(otu_table)-2)))

    countdata     =    otu_table %>%
      dplyr::select(-c(subject,group))
    
    #rownames(countdata) = (meta_dd$subject)
    
    count_table[[i]]= list(countdata = countdata)
  }
  
  names(count_table)  =   paste0("sim", 1:nsim)
  
  list(count_table = count_table, 
       meta_data  =  met_dd)
  
}


rrsim <- function(estimates = NULL, param = NULL, seed = NULL, nsim = 1){
  
  # Specify the model formula
  form <- count ~ 1 + us(1 + group | taxon) + rr(0 + taxon | subject, 2)
  
  if(!is.null(param)){
    
    ntaxa <- param$ntaxa
    nsubj <- sum(param$grp_pars$size_vec)
    
    # Precompute theta and betazi
    theta <- c(unlist(param$us_pars),
               get_theta(param$rr_pars$d, n=ntaxa, param$rr_pars$logsdvec))
    
    betazi <- rnorm(ntaxa, mean = param$zi_pars$mu, sd = param$zi_pars$sd)
    
    # Parameter list
    pars <- list(
      beta = param$beta_pars$beta, 
      theta = theta, 
      betadisp = param$beta_pars$betadisp, 
      betazi = betazi
    )
    
    # Generate data
    meta_dd <- met_data(ntaxa, param$grp_pars)
    
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
    
    res  =   otu_meta_fun(sim_count,meta_dd)
    
    res_list = list(otu_table =  res$count_table, 
                    meta_dd   =  res$meta_data,
                    count_long = sim_count)
    
    return(res_list)
  }elseif{
    
    dd_long     =   df_long(count_dd)
    dd_long     =   left_join(dd_long, meta_dd, by ="subject")  
    
    dds        =   DESeqDataSetFromMatrix(count_dd, meta_dd, ~group)
    dds        =   DESeq(dds,sfType ="poscounts",minReplicatesForReplace=Inf) 
    normalizer =   data.frame(normalizer = sizeFactors(dds)) %>% 
      rownames_to_column("subject")
    
    df      =    left_join(dd_long, normalizer, by ="subject")  
    ###########################################################
    form      =   count ~ 1 + us(1 + group|taxon) +  
      rr(0 + taxon | subject,2) + 
      offset(normalizer)
    
    ##Fit the model
    par_ctrl <- glmmTMBControl(
      parallel = list(n = 6, autopar = TRUE),
      optArgs  = list(maxit = 10000)
    )
    
    gprior  <- data.frame(prior = "gamma(2, 2.5)",
                          class = "theta_sd",
                          coef = "")
    
    options(glmmTMB_openmp_debug = TRUE)
    blas_set_num_threads(1)
    
    fit <- glmmTMB(formula = form, 
                     data = df, 
                     family = nbinom2, 
                     ziformula = ~1 + (1 | taxon),
                     prior = gprior, 
                     REML = TRUE, 
                     control = par_ctrl)
    
    pars <- list(
      beta = fit$beta, 
      theta = fit$theta, 
      betadisp = fit$betadisp, 
      betazi = fit$betazi
    )
    
    sim_count = simulate(fit)
    
    return(sim_count=sim_count,fit=fit)
    
  }
  
  
}


met_data = function(ntaxa, grp_pars){
  nsubj    =  sum(grp_pars$size_vec)
  metadata <- data.frame(subject = factor(seq(nsubj)),
                         group = rep(grp_pars$name_vec, grp_pars$size_vec))
  df <- expand.grid(taxon = factor(seq(ntaxa)), subject = factor(seq(nsubj)))
  df <- merge(metadata, df, by = "subject")
  df
}

##' @param d rank (dimension)
##' @param n full dimension (latent variables per group)
##' @param logsdvec vector of log-SDs of each factor
get_theta <- function(d, n, logsdvec) {
  
  mat <- matrix(0, nrow=n, ncol=d)
  ## replicate if length-1 ...
  if (length(logsdvec) == 1) logsdvec <- rep(logsdvec, d)
  stopifnot(length(logsdvec) == d)
  
  ## 1. pick values for each column where sum(x^2)==1
  for (i in 1:d) {
    r <- rnorm(n-i+1)
    rexp <- exp(r)
    mat[(i:n), i] <- sqrt(rexp/sum(rexp))
  }
  mat <- sweep(mat, 2, FUN = "*", exp(logsdvec))
  ## check consistency (could comment this out)
  
  stopifnot(all.equal(sqrt(colSums(mat^2)), exp(logsdvec)))

  theta <- c(
    mat[row(mat)==col(mat)],  ## diagonal elements
    mat[row(mat)>col(mat)]    ## below-diagonal elements
  )
  return(theta)
}


param    =   list(ntaxa = 100, 
                  beta_pars =  list(beta=1, betadisp=0),
                  zi_pars   =  list(mu = -2, sd = 1), 
                  us_pars   =  list(logsdvec = c(0.1, -0.2), corr = 0.5),
                  rr_pars   =  list(d = 2, logsdvec = c(0.1, -0.2)),
                  grp_pars = list(name_vec = c("control", "treat"), 
                                  size_vec = c(10,10))
)

mm=rrsim(param = param, seed = 101, nsim = 2)

View(mm)
