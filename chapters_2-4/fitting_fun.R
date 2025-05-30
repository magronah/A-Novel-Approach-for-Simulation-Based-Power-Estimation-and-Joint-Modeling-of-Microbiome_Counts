RR_fit  <- function(dflong, 
                    rank    =  2, 
                    nn      =  10,
                    gprior  = data.frame(prior = "gamma(2, 2.5)",
                                         class = "theta_sd",
                                         coef = "")){
  
  form <- as.formula(
    paste0("count ~ 1 + us(1 + group|taxon) + offset(normalizer) + rr(0 + taxon | subject,", rank, ")"))
  
  ##Fit the model
  par_ctrl <- glmmTMBControl(
    parallel = list(n = nn, autopar = TRUE)
  )
  
  options(glmmTMB_openmp_debug = TRUE)
  blas_set_num_threads(1)
  
  mod  <- tryCatch({
    glmmTMB(form,
            ziformula = ~1,
            data = dflong,
            family  = nbinom2,
            prior   = gprior,
            REML    = TRUE,
            control = par_ctrl)
  }, error = function(e) {
    message("Error in first attempt, trying again without prior...")
    glmmTMB(form, 
            ziformula = ~1,
            data = df,
            family  = nbinom2,
            REML    = TRUE,
            control = par_ctrl)
  })
  
  mod
}
##########################################################################
#' Title: compute the optimal number of components using parametric bootstrap
#'
#' @param logmean 
#' @param sig: significance level 
#' @param max.comp: maximum number of Gaussian components to check 
#' @param max.boot: maximum number of bootstraps  to test
#'
#' @return: ncomp: the optimal number of components for log mean abundance
#' @export
#'
#' @examples
optimal.comp <- function(logmean,sig=0.05,max.comp=4,max.boot=100){
  a <- boot.comp(y = logmean, max.comp = max.comp, B = max.boot, 
                 mix.type = "normalmix",epsilon = 1e-3)
  
  pvals=a$p.values; l=length(pvals)
  ncomp <- if (pvals[length(pvals)]<sig) length(pvals)+1 else length(pvals)
  return(ncomp)
}


##########################################################################

#' Title
#'
#' @param logmean: vector of log mean abundances 
#'
#' @return: components: optimal number of Gaussian components and 
#'          logmean_param: parameters for  Gaussian mixtures
#'          density_plot: plot comparing simulation from fit with observation
#' @export
#'
#' @examples
#' 
logmean_fit <- function(logmean){
  
  ncomp_opt = optimal.comp(logmean)
  if(length(ncomp_opt) == 0){stop("zero number of component")}
  
  if(ncomp_opt == 1){
    mixmdl = fitdistrplus::fitdist(logmean,"norm")
    
    param  = data.frame(sigma =  coef(mixmdl)[["sd"]],
                        mu  = coef(mixmdl)[["mean"]])
    
    sim_mean = rnorm(n = length(logmean),
                     mean = param$mu, 
                     sd   = param$sigma)
  }else{
    mixmdl = normalmixEM(logmean, k=ncomp_opt, maxrestarts=1e3)
    
    param  = data.frame(lambda =  mixmdl$lambda,
                        sigma =  mixmdl$sigma,
                        mu   =  mixmdl$mu)
    
    sim_mean=mixtools::rnormmix(n=length(logmean),lambda=mixmdl$lambda, 
                      mu =mixmdl$mu, 
                      sigma=mixmdl$sigma)
  }
  
  dat = data.frame(logmean)
  pp=ggplot(dat, aes(logmean)) +
    geom_histogram(aes(y=after_stat(density)),color=1, fill = "white") +
    geom_density(aes(colour="observation"), lwd=1) +
    geom_density(aes(x=sim_mean, colour="simulation"), lwd=1) +
    xlab(TeX("$\\log_2$(overall mean abundance)")) 
  
  list(logmean_param=param,
       density_plot=pp,
       components=ncomp_opt)
}


##################################################################
#' Title: fit non-linear model: a  + b/mean_abundance to dispersion estimates
#'
#' @param dispersion: dispersion estimates from deseq 
#' @param logmean vector of log mean abundance
#'
#' @return: param: parameters for the fit with confidence intervals
#' @export
#'
#' @examples
dispersion_fit <- function(dispersion,logmean){
  
  dat = data.frame(dispersion=dispersion,mean_abund = 2^logmean)
  
  param= nls(dispersion~ asymptDisp + extraPois/mean_abund, data = dat,
             start = list(asymptDisp = 0.1, extraPois = 0.1))
  
  dd=list(param=data.frame(asymptDisp= coef(param)[[1]],
                           extraPois = coef(param)[[2]]),
          confint=confint(param)) 
  dd
}

##################################################################
#' Title
#'
#' @param logmean: vector of log mean abundance 
#' @param logfoldchange: vector of log fold change 
#' @param ncore: number of cores to use 
#' @param max_sd_ord: the maximum order of polynomial function to fit to
#'                    standard deviation parameter. 
#'                    This must be either 1 (linear) or 2(quadratic)
#' @param max_np: maximum number of Gaussian components to check for 
#' @param minval: minimum value for DEoptim search
#' @param maxval: maximum value for DEoptim search 
#' @param NP:   
#' @param itermax:  
#' @param seed 
#'
#' @return
#' @export
#'
#' @examples
logfoldchange_fit = function(logmean,logfoldchange,ncore = 2,
                             max_sd_ord = 2, max_np=5,
                             minval = -5, maxval = 5,
                             itermax = 100,NP=800, seed  = 100){
  
  res  =  list(); cnt = 0
  
  for(sd_ord in 1:max_sd_ord){
    
    for(np in 2:max_np){
      
      l   =  (np-1)+2*np+(sd_ord+1)*np 
      #NP = 10*l 
      #per the recommendation form the authors of DEOptim
      # Set the number of parents NP to 10 times the number of parameters
      #https://cran.r-project.org/web/packages/DEoptim/DEoptim.pdf
      ################################################################
      set.seed(seed)
      
      cl <- makeCluster(ncore)
      fun_names <- c("polyfun", "dnormmix", "dnormmix0", "genmixpars")
      clusterExport(cl, c("logmean", "logfoldchange", fun_names))
      
      opt  <- DEoptim(fn = nllfun, lower = rep(minval, l), upper = rep(maxval, l),
                      vals = logfoldchange, logmean = logmean, np = np, sd_ord = sd_ord,
                      control = DEoptim.control(NP = NP, itermax = itermax,
                                                cluster = cl))
  
      pn <- gen_parnames(np = np, sd_ord = sd_ord)
      par <- setNames(opt$optim$bestmem, pn)
      aic = 2*length(par) + 2*opt$optim$bestval

      
      cnt = cnt + 1
      res[[cnt]]  =  list(par = par, np = np, sd_ord= sd_ord,aic = aic)
    }
  }
  aic_values <- sapply(res, function(x) x$aic)
  index_of_min_aic <- which.min(aic_values)
  
  res[[index_of_min_aic]]
} 

##################################################################
gam_fit <- function(deseq_est_list,
                    true_lfoldchange_list,
                    true_lmean_list,
                    grid_len = 50,
                    alpha_level=0.1){
  
  
       # p_val   =  (deseq_est_list  
       #                   %>% setNames(paste0("padjust", 
       #                                     1:length(deseq_est_list))) 
       #                   %>%  purrr::map_dfr(pull, padj))
  
   p_val   = foreach(k = 1:length(deseq_est_list),.combine = "c") %do%{
        deseq_est_list[[k]]$padj
  }
  
  pval_reject        =   (!is.na(p_val) & p_val < alpha_level)
  true_lfoldchange    =    unlist(true_lfoldchange_list)
  true_lmean_abund    =    unlist(true_lmean_list)
 
             comb     =   tibble(lmean_abund  =  true_lmean_abund,
                                 abs_lfc   =  abs(true_lfoldchange),
                                 pval_reject  =  as.numeric(pval_reject))
  
         #fit scams
         fit_2d       =    scam(pval_reject ~ s(lmean_abund, abs_lfc, bs="tedmi"),
                                 data = comb, family = binomial)
         
              pp      =   with(comb,
                               expand.grid(lmean_abund = seq(min(lmean_abund),
                                                             max(lmean_abund),
                                                            length  = grid_len),
                                           abs_lfc   =  seq(min(abs_lfc),
                                                         max(abs_lfc),
                                                         length  =  grid_len)))
  #predict power
  pp$power <- predict(fit_2d, newdata = pp,type = "response")
  
  p=list(combined_data = comb, power_estimate = pp, fit_2d=fit_2d)
  p
}



#' Title
#'
#' @param deseq_estimate_list list containing simulations for each sample size
#' @param true_lfoldchange_list 
#' @param true_lmean_list 
#' @param alpha 
#'
#' @return
#' @export
#'
#' @examples
power_fun_ss <- function(deseq_est_list,
                         true_logfoldchange,
                         true_logmean,
                      sample_vec,
                      alpha=0.1,notu){ 

  # concatenate all p-values from all the sample size
  # est_list <- do.call("c", deseq_est_list)
  # p_val <- do.call("c", lapply(est_list, function(est) est$deseq_estimate$padj))
  
  # concatenate all log foldchange, logmean from all the sample size
  #true_logfoldchange <- do.call("c", do.call("c", sim_logfoldchange_list))
  #true_logmean  <-  do.call("c", do.call("c", sim_logmean_list))
  
  #deseq_est_list = deseq_sample_size[[1]]
  # find p-values that were rejected
  p_val =  deseq_est_list$padj
  
  pval_reject   =   (!is.na(p_val) & p_val < 0.1)
  #length(p_val)
  # create a table with all the information
  comb   =   tibble(lmean_abund  =   true_logmean,
                    abs_lfc      =   abs(true_logfoldchange),  
                    pval_reject  =   as.numeric(pval_reject))
  
  comb$sample_size = deseq_est_list$sample_size #rep(sample_vec, each = nsim*notu)
  #' fit GAM with covariates as tensor product (ie,interaction between
  #' log mean abundance and absolute log fold changes 
  #' and then a spline for the sample sizes
  #' log mean abundance and log fold changes are related directly,hence the 
  #' interaction but sample size is not quite related to log mean abundance
  #' and log fold changes directly
  
  df =  length(sample_vec) -1 # degrees of freedom
  comb$lss = log2(comb$sample_size)
  
  # fit_3d <- scam(pval_reject ~ s(lmean_abund, abs_lfc,bs="tedmi") + s(lss, k = df),
  #                data = comb, family = binomial)
  # fit_3d <- scam(pval_reject ~ s(lmean_abund, abs_lfc,bs="tedmi") + s(sample_size,bs="mpi"),
  #                data = comb, family = binomial)
  
  fit_3d <- scam(pval_reject ~ s(lmean_abund, abs_lfc,bs="tedmi") +
                              s(sample_size,lmean_abund,bs="tedmi") +
                              s(sample_size,abs_lfc,bs="tedmi"),
                             data = comb, family = binomial)
  
  list(combined_data=comb, gam_mod = fit_3d)
}


