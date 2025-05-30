zinbmm_confint = function(mod, mean_count,
                          group_label = "grouptreat",
                          conf_level = .95){
  
  sd_err =  as.numeric(unlist(lapply(mod$fit,function(x) 
  {summary(x)$tTable[group_label, "Std.Error"]})))
  
  est    =   as.numeric(unlist(lapply(mod$fit,function(x) 
  {summary(x)$tTable[group_label, "Value"]})))
  
  pvalue =  as.numeric(unlist(lapply(mod$fit,function(x) 
  {summary(x)$tTable[group_label, "p-value"]})))
  
  dd         =   data.frame(est_param = est)
  z_score    =   qnorm(conf_level + (1 - conf_level)/2)
  
  dd$lwr     =   est  -  z_score*sd_err
  dd$upr     =   est  +  z_score*sd_err
  dd$width   =   dd$upr  - dd$lwr
  
  dd$pvalue      =   pvalue  
  dd$param_name  =   names(mod$fit)
  dd
}



deseq_wald_confint = function(deseq_est,conf_level = .95){
  foldchange =   log(2)*deseq_est$log2FoldChange
  dd         =   data.frame(est_param = foldchange)
  z_score    =   qnorm(conf_level + (1 - conf_level)/2)
  sd_err     =   log(2)*deseq_est$lfcSE
  dd$lwr     =   foldchange  -  z_score*sd_err
  dd$upr     =   foldchange  +  z_score*sd_err
  dd$width   =   dd$upr  - dd$lwr
  dd$param_name  =  deseq_est$param_name
  dd$padj        =   deseq_est$padj  
  dd
}


#'
#' @param mod 
#' @param ntaxa 
#' @param conf_level 
#'
#' @return
#' @export
#' 
#' @examples
wald_confint = function(mod, conf_level = .95, 
                        in_sd = 1, ntaxa,
                        mean_count, mod_name,
                        sdreport_path,
                        fullsdr_path){
  
  grp_ind   =   seq(2, 2*(ntaxa), 2)
  #pred    =   predict(mod, type = "latent", se.fit = TRUE) 
  #est     =   pred$fit[grp_ind]
  #sd_err  =   pred$se.fit[grp_ind]
  ########################################################
  ss <- TMB::sdreport(mod$obj, getJointPrecision = TRUE)
  saveRDS(ss, file = paste0(sdreport_path,mod_name,".rds"))
  #ss  <-  readRDS(paste0(path,"sdreport_",mod_name,".rds"))
  #########################################################
  ss$jointPrecision <-  as(ss$jointPrecision, "sparseMatrix")
  inverse_mat       <-  solve(ss$jointPrecision)
  
  if(any(diag(inverse_mat) < 0)){
    precision_nearPD <-  as.matrix(nearPD(ss$jointPrecision)$mat)
    precision   <-  as(precision_nearPD, "sparseMatrix")
    se_vec <- sqrt(diag(solve(precision)))
  }else{
    se_vec <- sqrt(diag(inverse_mat))
  }
  # Run nearPD
  #se_vec <- sqrt(diag(solve(nearPD(ss$jointPrecision))))
  # allFit is another option
  # https://github.com/glmmTMB/glmmTMB/blob/master/misc/allFit.R
  #se_vec <- sqrt(diag(solve(ss$jointPrecision)))
  #start_method in glmmTMBControl jitter.sd
  saveRDS(se_vec, file = paste0(fullsdr_path,mod_name,".rds"))
  #########################################################
  full_est  =   mod$fit$parfull
  est       =   full_est[names(full_est) == "b"][grp_ind]
  sd_err    =   se_vec[grp_ind]
  

  # Calculate z-score for the desired confidence level
  z_score    =    qnorm(conf_level + (1 - conf_level)/2)
  
  # Calculate confidence intervals
  dd      =    data.frame(est_param =   est,
                          lwr       =   est - z_score*in_sd*sd_err,
                          upr       =   est + z_score*in_sd*sd_err)
  # Calculate p-values
  z_stat   =  est / sd_err
  p_values =  2 * (1 - pnorm(abs(z_stat)))  # Two-tailed p-value
  
  dd$pvalue  =  p_values
  dd$width   =   dd$upr  - dd$lwr
  dd$mean_count  =  as.numeric(mean_count)
  dd$param_name  =  names(mean_count)
  dd
}
#########################################################
generate_increasing_total_variance <- function(tt_logvar, tt_cor, rr_logvar, steps = 10) {
  # Ensure input validity
  if (length(tt_logvar) != 2) {
    stop("tt_logvar must each have exactly two elements.")
  }
  
  # Base values for the fixed components
  v1 <- tt_logvar[1]
  v2 <- tt_logvar[2]
  cor <- tt_cor
  
  # Proportions for rr components
  proportions <- seq(0, 1, length.out = steps)
  
  # Generate results
  results <- lapply(proportions, function(p_rr) {
    # Calculate tt contributions (remain constant)
    tt_sd1_cont <- v1
    tt_sd2_cont <- v2 / 2
    tt_cor_cont <- cor
    
    # Adjust rr contributions (increasing)
    #rr1_sd_cont <- p_rr * rr_logvar[1]
    #rr2_sd_cont <- p_rr * rr_logvar[2]
    rr_cont  <- p_rr *rr_logvar
    
    # Calculate the total variance (increasing)
    total_variance <- tt_sd1_cont + tt_sd2_cont + tt_cor_cont + sum(rr_cont)
    # rr1_sd_cont + rr2_sd_cont
    
    # Create a data frame for this iteration
  ddf = data.frame(
      tt_sd1_cont = tt_sd1_cont,
      tt_sd2_cont = tt_sd2_cont,
      tt_cor_cont = tt_cor_cont,
      t(rr_cont),
      #rr1_sd_cont = rr1_sd_cont,
      #rr2_sd_cont = rr2_sd_cont,
      total_var = total_variance,
      prop_rr = p_rr
    )
  
  colnames(ddf)  =  c("tt_sd1_cont","tt_sd2_cont","tt_cor_cont",
                      paste0("rr",1:(length(rr_cont)),"_sd_cont"),
                           "total_var", "prop_rr")
  ddf
  
  })
  
    # Combine results into a single data frame
  results_df <- do.call(rbind, results)
  return(results_df)
}


generate_variance_proportions <- function(total_logvar, tt_logvar, tt_cor, rr_logvar, steps = 10) {
  # Ensure input validity
  if (length(tt_logvar) != 2) {
    stop("tt_logvar must each have exactly two elements.")
  }
  
  # Base values for the components
  v1 <- tt_logvar[1]
  v2 <- tt_logvar[2]
  cor <- tt_cor
  
  #rr1 <- rr_logvar[1]
  #rr2 <- rr_logvar[2]
  
  # Generate proportions for v1, v2, rr1, rr2, and cor
  proportions <- seq(0, 1, length.out = steps)
  
  # Generate results
  results <- lapply(proportions, function(p_c2) {
    # Proportion for c1 (which will be decreasing)
    p_c1 <- 1 - p_c2
    
    # Split the variance contribution between the components:
    # Adjust the values of v1, v2, rr1, rr2, and cor based on the proportions
    
    # Adjust the values for tt_logvar and rr_logvar
    tt_sd1_cont <- p_c1 * v1
    tt_sd2_cont <- p_c1 * (v2 / 2) # The second component from tt_logvar has half the weight
    tt_cor_cont <- p_c1 * cor
    
    rr_cont  <- p_c2 *rr_logvar
    #rr1_contrib <- p_c2 * rr1
    #rr2_contrib <- p_c2 * rr2
    
    # Calculate the total variance as the sum of contributions
    total_variance <- tt_sd1_cont + tt_sd2_cont + tt_cor_cont + sum(rr_cont)
    #rr1_contrib + rr2_contrib
    
    # Ensure the total variance remains constant (close to total_logvar)
    if (abs(total_variance - total_logvar) > 1e-5) {
      scaling_factor <- total_logvar / total_variance
      tt_sd1_cont <- tt_sd1_cont * scaling_factor
      tt_sd2_cont <- tt_sd2_cont * scaling_factor
      tt_cor_cont <- tt_cor_cont * scaling_factor
      rr_cont <-  rr_cont * scaling_factor
      #rr1_contrib <- rr1_contrib * scaling_factor
      #rr2_contrib <- rr2_contrib * scaling_factor
    }
    
  dd1=  data.frame(
      tt_sd1_cont = tt_sd1_cont,
      tt_sd2_cont = tt_sd2_cont,
      tt_cor_cont = tt_cor_cont,
      total_var =  tt_sd1_cont + tt_sd2_cont + tt_cor_cont + sum(rr_cont),
      prop_tt = p_c1,
      prop_rr = p_c2
    )
  
  dd2=  data.frame(
    t(rr_cont),
    total_var =  tt_sd1_cont + tt_sd2_cont + tt_cor_cont + sum(rr_cont)
  )
  
  colnames(dd2)  =  c(paste0("rr",1:(length(dd2)-1),"_sd_cont"),"total_var")
  dd =right_join(dd1,dd2, by="total_var")
  dd
  })
  
  # Combine results into a single data frame
  results_df <- do.call(rbind, results)
  return(results_df)
}

df_long = function(dd, otu_names = "sp", subject_name = "subject", ntaxa){
  if(ncol(dd) != ntaxa){
    dd  =  as.data.frame(t(dd))
  }else{
    dd  =  data.frame(dd)
  }
  df   =  dd %>%
    rownames_to_column(subject_name)
  ddd =   pivot_longer(df,
                       cols = starts_with(otu_names),
                       names_to  = otu_names,
                       values_to = "count")
  names(ddd)  = c(subject_name,"taxon","count")
  ddd
}

###########################################
## goal: pick 'theta' parameters for a reduced-rank model in a sensible way
##' @param d rank (dimension)
##' @param n full dimension (latent variables per group)
##' @param logsdvec vector of log-SDs of each factor
get_theta_rr <- function(d, n, logsdvec, seed = NULL) {
  mat <- matrix(0, nrow=n, ncol=d)
  ## replicate if length-1 ...
  if (length(logsdvec) == 1) logsdvec <- rep(logsdvec, d)
  ## otherwise fail
  stopifnot(length(logsdvec) == d)
  ## 1. pick values for each column where sum(x^2)==1
  for (i in 1:d) {
    ## we don't care about identifiablity here, so we
    ## can pick as many N(0,1) values as we need and rescale
    ## (unlike if we were trying to *estimates* these parameters,
    ##  would need to constrain one value for identifiability,
    ##  e.g. set the first element to 0 (without loss of generality?)
    set.seed(seed)
    r <- rnorm(n-i+1)
    # the other terms are coming from the rnorm, and they do not 
    # affect the standard deviations of the reduced rank term
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


gam_fit <- function(pvalue, effect_size, mean_count,
                    grid_len = 100, alpha_level = 0.05){
  
  pval_reject  =   (!is.na(pvalue) & pvalue < alpha_level)
  
  comb      =   tibble(lmean_count  =  log(mean_count),
                       abs_lfc      =  abs(effect_size),
                       pval_reject  =  as.numeric(pval_reject))


  #fit scams
  fit_2d       =    scam(pval_reject ~ s(lmean_count, abs_lfc, bs="tedmi"),
                         data = comb, family = binomial)
  
  #fit_2d       =    scam(pval_reject ~ s(lmean_count, abs_lfc),
  #                       data = comb, family = binomial)

  pp   =   with(comb,
                expand.grid(lmean_count = seq(min(lmean_count),
                                              max(lmean_count),
                                              length  = grid_len),
                            abs_lfc   =  seq(min(abs_lfc),
                                             max(abs_lfc),
                                             length  =  grid_len)))
  
  #predict power
  pp$power <- predict(fit_2d, newdata = pp,type = "response")
  
  p=list(combined_data = comb, power_estimate = pp, fit_2d=fit_2d)
  p
}

#################################################################
pvalue_fun =  function(dd){
  apply(dd, 1, function(x){
    2*min(c(mean(x < 0), mean(x > 0)))})
}


fishers_combined <- function(pvals) {
  pchisq(-2 * sum(log(pvals)), df = 2 * length(pvals), lower.tail = FALSE)
}

stouffer_combine <- function(pvals) {
  sum_qnorm <- sum(qnorm(1 - pvals))
  combined_z <- sum_qnorm / sqrt(length(pvals))
  pnorm(combined_z, lower.tail = FALSE)
}
###############################################################
power_pred =  function(object, newdata){
  
  power   =  data.frame(power = predict(object, 
                                        newdata = newdata, 
                                        type = "response"))
  cbind(newdata, power)
}
###############################################################

power_predict_dd =  function(mod_obj_list, newdata){
  
  pow  =  list()
  len  =  length(mod_obj_list)
  
  for(i in 1:len){
    pow[[i]]  =   predict(mod_obj_list[[i]], newdata = newdata, type = "response")
  }
  
  if(!is.null(names(mod_obj_list))){
    names(pow)   =   names(mod_obj_list)
  }
  
  power_dd      =   data.frame(power = unlist(pow),
                               model =  rep(names(pow), 
                                            each = nrow(newdata)))
  
  
  pow_dd        =   cbind(newdata, power_dd)
  pow_dd
}


gam_fit2 <- function(pval,lfoldchange,lmean_abund,
                     grid_len = 100){
  
  #pval_reject  =   (!is.na(pval) & pval < alpha_level)
  
  comb      =   tibble(lmean_abund  =  lmean_abund,
                       abs_lfc      =  abs(lfoldchange),
                       pval_reject  =  as.numeric(pval))
  
  #fit scams
  fit_2d       =    mgcv::gam(pval_reject ~ te(lmean_abund, abs_lfc),
                              data = comb, family = binomial)
  
  pp   =   with(comb,
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

otu_meta_lst_fun = function(res_dd1){
  
  res_dd2   = list()
  for(i in 1:length(res_dd1)){
    dd  =  res_dd1[[i]]
    
    otu_table <- dd %>%
      dplyr::select(subject, taxon, count, group)
    
    otu_table <- spread(otu_table, key = taxon, value = count) 
    colnames(otu_table) <- c("subject","group", paste0("taxon",1:ntaxa))
    
    met_data    =   otu_table %>%
      dplyr::select(subject,group)
    
    countdata     =    otu_table %>%
      dplyr::select(-c(subject,group))
    
    rownames(countdata) = (met_data$subject)
    
    res_dd2[[i]]= lst(countdata,met_data)
  }
  
  names(res_dd2)  =   paste0("sim", 1:nsim)
  res_dd2
}


custom_theme <- function(n) {
  theme_bw(base_size = n) +
    theme(
      plot.title = element_text(hjust = 0.5,size = n, family = "Roboto",color = "black"),
      text = element_text(size = n, family = "Roboto",color = "black"),
      axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
      axis.text.y = element_text(family = "Roboto", size = n, color = "black")
    )
}
#################################################################
load_models <- function(path, filenames) {
  file_paths <- paste0(path, filenames)
  mod_list  <- lapply(file_paths, readRDS)
}
#################################################################
load_data <- function(path, alpha = 0.05) {
 
  true_param =  readRDS(paste0(path,"true_param.rds"))
  ####################################################
  RR     =  readRDS(paste0(path, "rr.rds"))
  RRzi   =  readRDS(paste0(path, "rrzi.rds"))
  US     =  readRDS(paste0(path, "us.rds"))
  USzi   =  readRDS(paste0(path, "uszi.rds"))
  NB    =  readRDS(paste0(path, "nbmm.rds"))
  ZNB   =  readRDS(paste0(path, "zinbmm.rds"))
  
  # rownames(nbmm) =   rownames(zinbmm)  =  rownames(rrzi)
deseqL2  =  readRDS(paste0(path, "deseq.rds"))
deseq_noShrinkL2  =  readRDS(paste0(path, "deseq_noShrink.rds"))

     DE       =    deseqL2*log(2)
     DE_noSh  =    deseq_noShrinkL2*log(2)

  ####################################################
  rrl    =   dd_long(RR,   true_param,label="RR")
  rrzil  =   dd_long(RRzi, true_param,label="RRzi")
  usl    =   dd_long(US,  true_param,label="US")
  uszil  =   dd_long(USzi, true_param,label="USzi")
  deseql =   dd_long(DE, true_param,label="DE")
  deseq_noShrinkl =   dd_long(DE_noSh, true_param,label="DE_noSh")
  
  nbmml  =   dd_long(NB,  true_param,label="NB")
 zinbmml =   dd_long(ZNB, true_param,label="ZNB")
  ####################################################
 dd = lst(true_param,RR, RRzi, US, USzi, DE, DE_noSh, NB, ZNB)
  ##convert to long format
 long_dd = list(RR = rrl, RRzi = rrzil, US = usl, 
                USzi =  uszil, DE  =  deseql,
                DE_noSh  =  deseq_noShrinkl,
                NB  = nbmml, ZNB  = zinbmml)

  confint = list(
    RR      =   para_confint(rrl, true_param, alpha = alpha),
    RRzi    =   para_confint(rrzil, true_param, alpha = alpha),
    US      =   para_confint(usl, true_param, alpha = alpha),
    USzi    =   para_confint(uszil, true_param, alpha = alpha),
    DE      =   para_confint(deseql, true_param, alpha = alpha),
    DE_noSh =   para_confint(deseq_noShrinkl, true_param, alpha = alpha),
    NB      =   para_confint(nbmml, true_param, alpha = alpha),
    ZNB     =   para_confint(zinbmml, true_param, alpha = alpha)
  ) 
 
  
  error = list(
    RR      =   error_cal(RR, true_param, model = "RR"),
    RRzi    =   error_cal(RRzi, true_param, model = "RRzi"),
    US      =   error_cal(US, true_param, model = "US"),
    USzi    =   error_cal(USzi, true_param, model = "USzi"),
    DE      =   error_cal(DE, true_param, model = "DE"),
    DE_noSh   =   error_cal(DE_noSh, true_param, model = "DE_noSh"),
    NB     =   error_cal(NB, true_param, model = "NB"),
    ZNB  =   error_cal(ZNB, true_param, model = "ZNB")
  ) 
  
  lst(dd,long_dd, confint, error)

}


error_cal <- function(model_est_dd, true_param_dd, model,
                      lwr  =  0.025, upr = 0.975,  scale = 1.96, 
                      sdrr = TRUE) {
  
  est_dd  =   model_est_dd %>%
    rownames_to_column(var = "param_name")
  merge_ddd     =   left_join(true_param_dd, est_dd, by = "param_name")  
  merge_dd      =   merge_ddd  %>% dplyr::select(-param_name) 
  errr        =   data.frame(t(apply(merge_dd, 1, 
                                   function(x){(x["true_param"] -  x)})))
  
  error        =    errr[, !colnames(errr) %in% "true_param"]
  variance     =    apply(error, 1, var)
  
  df1       =   data.frame(param_name  =   merge_ddd$param_name,
                           true_param  =   merge_ddd$true_param,
                           bias        =   rowMeans(error),
                           mse         =   rowMeans(error^2),
                           variance    =   variance)
  
  df2     =    data.frame(average_value  =   mean(rowMeans(error)))
  df3     =    data.frame(average_value  =   mean(sqrt(rowMeans(error^2))))
  df4     =    data.frame(average_value  =   mean(variance))
  
  
  df1$model   =  rep(model,nrow(df1)) 
  df2$model   =  rep(model,nrow(df2)) 
  df3$model   =  rep(model,nrow(df3)) 
  df4$model   =  rep(model,nrow(df4)) 
  
 
  if(sdrr){
    df2$mean     =   mean(df1$bias)
    df2$lwr     =   mean(df1$bias) -  sd(df1$bias)/sqrt(length(df1$bias))
    df2$upr     =   mean(df1$bias) + sd(df1$bias)/sqrt(length(df1$bias))
    ##############
    rrm        =  sqrt(df1$mse)
    df3$mean  =   mean(rrm)
    df3$lwr   =   mean(rrm) -  scale*sd(rrm)/sqrt(length(rrm))
    df3$upr   =   mean(rrm) +  scale*sd(rrm)/sqrt(length(rrm))
    
    df4$mean     =   mean(df1$variance)
    df4$lwr   =   mean(df1$variance) -  scale*sd(df1$variance)/sqrt(length(df1$variance))
    df4$upr   =   mean(df1$variance) +  scale*sd(df1$variance)/sqrt(length(df1$variance))
    
    res     =   lst(full_summary_dd =   df1, error, 
                    avg_bias      =   df2, 
                    avg_mse     =   df3,
                    avg_var    =   df4)
    
    return(res)
  }else{
    df2$lwr     =   quantile(df1$bias, lwr)
    df2$upr     =   quantile(df1$bias, upr)
    
    ##############
    df3$lwr   =   quantile(sqrt(df1$mse) , lwr)
    df3$upr   =   quantile(sqrt(df1$mse) , upr)
    
    df4$lwr   =   quantile(df1$variance , lwr)
    df4$upr   =   quantile(df1$variance , upr)
    
    res     =   lst(full_summary_dd =   df1, error, 
                    avg_bias      =   df2, 
                    avg_mse     =   df3,
                    avg_var    =   df4)
    
    return(res)
  }
  
}


###########################################################
err_extract = function(data_list, extract_name){
  summaries <- lapply(data_list$error, function(x) x[[extract_name]])
  do.call(rbind, summaries)
}
###########################################################
reorganise_dd =  function(dd, name){
  
  rr_row <- dd[dd$model == "rr", ]
  rr_rep <- rr_row[rep(1, nrow(dd) - 1), ] 
  rr_rep$type  = paste0(name,1:nrow(rr_rep))
  
  other_mod  = dd %>%
    filter(model != "rr")  %>%  
    dplyr::arrange(average_value) %>%
  mutate(type = paste0(name, row_number()))
  
  result <- rbind(rr_rep, other_mod) 
  result
}
###########################################################
#' Title
#'
#' @param est_wide 
#'
#' @return
#' @export
#'
#' @examples
dd_long  =  function(dd_wide,dd_true_param,label="rr"){
  dd = (dd_wide
        |> rownames_to_column("param_name")
        |> pivot_longer(-param_name, names_to="sim", values_to = "estimate")
  )
  dd$model  =  rep(label,nrow(dd))
  df             =  left_join(dd,dd_true_param, by = "param_name")
  df
}



para_confint  = function(est_data, true_dd,  alpha = 0.05){
  dd = (est_data
        |> group_by(param_name)
        |> summarise(lwr = quantile(estimate, alpha/2),
                     upr = quantile(estimate, 1-alpha/2), 
                     average_estimate  =  mean(estimate))
        |> mutate(param_name = factor(param_name, levels = param_name))
  )
  
  dd$model  =   rep(unique(est_data$model), nrow(dd))
  ddd      =   left_join(dd, true_dd, by = "param_name")
  
  ddd    =  ddd %>% 
    
  mutate(
    param_name = factor(param_name, levels = unique(param_name[order(true_param)])),
    coverage = ifelse(true_param > lwr & true_param < upr, 1, 0)
    )
  ddd$CI_width =   ddd$upr - ddd$lwr
  
  ddd
}

## goal: pick 'theta' parameters for a reduced-rank model in a sensible way
##' @param d rank (dimension)
##' @param n full dimension (latent variables per group)
##' @param logsdvec vector of log-SDs of each factor
get_theta_corrRR <- function(d, n, logsdvec) {
  mat <- matrix(0, nrow=n, ncol=d)
  if (length(logsdvec) == 1) logsdvec <- rep(logsdvec, d)
  stopifnot(length(logsdvec) == d)
  for (i in 1:d) {
    r <- rnorm(n-i+1)
    rexp <- exp(r)
    mat[(i:n), i] <- sqrt(rexp/sum(rexp))
  }
  mat <- sweep(mat, 2, FUN = "*", exp(logsdvec))
  stopifnot(all.equal(sqrt(colSums(mat^2)), exp(logsdvec)))
  theta <- c(
    mat[row(mat)==col(mat)],  ## diagonal elements
    mat[row(mat)>col(mat)]    ## below-diagonal elements
  )
  return(theta)
}



normalizer_fun <- function(countdata,met_data,alpha_level=0.1,ref_name="NT",
                     minReplicatesForReplace = Inf, 
                     cooksCutoff = FALSE,
                     independentFiltering = FALSE,
                     do_shrinkage =  "yes", 
                     design   = ~group,
                     shrinkage_method="normal",
                     ntaxa){
  
  #check otu table is in otu by samples format
  if(nrow(countdata) != ntaxa){
    countdata = t(countdata)
  }
  
  #remove samples with zeros for all taxa (if any such sample exist)
  keep <- (colSums(countdata) > 0)
  countdata = countdata[,keep]
  met_data= met_data[keep, ]
  
  # call deseq
  dds <- DESeqDataSetFromMatrix(countdata,met_data, design = design)
  dds$group <- relevel(dds$group, ref = ref_name)
  
  dds <- DESeq(dds,sfType ="poscounts",
               minReplicatesForReplace = minReplicatesForReplace) 
  sizeFactors(dds) 
}




deseqfun <- function(countdata,met_data,alpha_level=0.1,ref_name="NT",
                     minReplicatesForReplace = Inf, 
                     cooksCutoff = FALSE,
                     independentFiltering = FALSE,
                     do_shrinkage =  "yes", 
                     design   = ~group,
                     shrinkage_method="normal",
                     ntaxa){
  
  #check otu table is in otu by samples format
  if(nrow(countdata) != ntaxa){
    countdata = t(countdata)
  }
  
  #remove samples with zeros for all taxa (if any such sample exist)
  keep <- (colSums(countdata) > 0)
  countdata = countdata[,keep]
  met_data= met_data[keep, ]
  
  # call deseq
  dds <- DESeqDataSetFromMatrix(countdata,met_data, design = design)
  dds$group <- relevel(dds$group, ref = ref_name)
  
  dds <- DESeq(dds,sfType ="poscounts",
               minReplicatesForReplace = minReplicatesForReplace) 
  
  res <- results(dds, cooksCutoff=cooksCutoff, 
                 independentFiltering=independentFiltering,
                 alpha = alpha_level)
  
  if(do_shrinkage == "no"){
    reslt   <-   res
   }else{
      reslt <- lfcShrink(dds, res=res, coef=2, type=shrinkage_method)
   }
  
  deseq_est = data.frame(reslt)
  deseq_est$dispersion = dispersions(dds)
  deseq_est$intercept  = coef(dds)[, "Intercept"]
  
  deseq_dd   =  deseq_est  %>% 
    rownames_to_column(var = "param_name")
  deseq_dd

}

# extract otu table and metadata
norm_fun <- function(dd){
  
  otu_tab  =   long_dd = meta_data = list()
  
  nsubj =  unique(dd$subject)
    t   =  unique(dd$time)
  
  for(i in 1:length(t)){
    wide_dd <- dd %>% 
      subset(time  == t[i]) %>% 
      dplyr::select(subject, taxon, count, group) %>%
      spread(key = taxon, value = count) %>%
      setNames(c("subject","group", paste0("taxon",1:ntaxa))) 
    
    otu_table  =  wide_dd %>%
      dplyr::select(paste0("taxon",1:ntaxa))
    
    metadata   =  wide_dd %>% 
      dplyr::select(c("subject","group"))
    
             otu_count  =   t(otu_table)
    colnames(otu_count) =   paste0("subject",nsubj)
    
    dds        =   DESeqDataSetFromMatrix(otu_count,metadata, ~group)
    dds        =   DESeq(dds,sfType ="poscounts",minReplicatesForReplace=Inf) 
    normalizer =   sizeFactors(dds) # one for each subject
    
    
    normalise_dd  =    data.frame(normalizer, subject = wide_dd$subject)
    dd_res        =    left_join(dd,normalise_dd, by ="subject")  
    
    long_dd[[i]]    =   dd_res
    otu_tab[[i]]    =   otu_count
    meta_data[[i]]  =   metadata
    }

  names(otu_tab)  =  names(long_dd)  =  paste0("time", t)
  list(otu_tab    =  otu_tab, long_dd= long_dd, metadata = meta_data)
}



otu_meta_lst_fun = function(res_dd1){
  
  res_dd2   = list()
  for(i in 1:length(res_dd1)){
    dd  =  res_dd1[[i]]
    
    otu_table <- dd %>%
      dplyr::select(subject, taxon, count, group)
    
    otu_table <- spread(otu_table, key = taxon, value = count) 
    colnames(otu_table) <- c("subject","group", paste0("taxon",1:ntaxa))
    
    met_data    =   otu_table %>%
      dplyr::select(subject,group)
    
    countdata     =    otu_table %>%
      dplyr::select(-c(subject,group))
    
    rownames(countdata) = (met_data$subject)
    
    res_dd2[[i]]= lst(countdata,met_data)
  }
  
  names(res_dd2)  =   paste0("sim", 1:nsim)
  res_dd2
}


meta_data = function(ntaxa, nsubj){
  metadata <- data.frame(subject = factor(seq(nsubj)),
                         group = rep(c("control", "treat"), each = nsubj/2))
  df <- expand.grid(taxon = factor(seq(ntaxa)), subject = factor(seq(nsubj)))
  df <- merge(metadata, df, by = "subject")
  df
}


get_theta_corr1 <- function(ntaxa,nsubject, mat= NULL, seed = NULL) {
  if(!is.null(mat)){C  <- mat}
  else{set.seed(seed); C <- get_corr1(ntaxa, nsubject,seed = seed)}
  C <- nearPD(C)$mat
  scale <- sqrt(fastmatrix::ldl(as.matrix(C))$d)
  cc2 <- chol(C) %*% diag(1/scale)
  cc2[upper.tri(cc2)]
}

get_corr1 <- function(ntaxa, nsubject, seed = NULL) {
  
  set.seed(seed)
  synthetic_data     <-   huge.generator(n = nsubject, d = ntaxa, graph = "scale-free")
  covariance_matrix  <- synthetic_data$sigma
  correlation_matrix <- cov2cor(covariance_matrix)
  
  return(correlation_matrix)
}


get_corr <- function(ntaxa, nsubject, sparse_thresh = 0.05, seed = NULL) {
  
  set.seed(seed)
  synthetic_data     <-   huge.generator(n = nsubject, d = ntaxa, graph = "scale-free")
  covariance_matrix  <- synthetic_data$sigma
  correlation_matrix <- cov2cor(covariance_matrix)
  
  #correlation_matrix[abs(correlation_matrix) < sparse_thresh] <- 0
  return(correlation_matrix)
}


get_theta_corr <- function(ntaxa,nsubject, mat= NULL, seed = NULL) {
  if(!is.null(mat)){C  <- mat}
  else{set.seed(seed); C <- get_corr(ntaxa, nsubject,seed = seed)}
  C <- nearPD(C)$mat
  scale <- sqrt(fastmatrix::ldl(as.matrix(C))$d)
  cc2 <- chol(C) %*% diag(1/scale)
  cc2[upper.tri(cc2)]
}


metadata <- function(ntaxa, nIndiv, ntime){
  d <- expand.grid(taxon = factor(1:ntaxa),
                   subject = factor(rep(1:nIndiv,ntaxa)))[1:(ntaxa*nIndiv), ]
  expdes <- data.frame(subject = factor(1:nIndiv), 
                       group=rep(c("control","treatment"), 
                                 each = nIndiv/2))
  dat0 <- left_join(d, expdes, by = "subject")
  
  dd <- do.call("rbind", replicate(n=ntime, dat0, simplify = FALSE))
  dd$time = rep(1:ntime, each=ntaxa*nIndiv)
  dd$nugget <- factor(1:nrow(dd)) 
  dd
}


##' I simulated logSD from a normal distribution 
##' 

get_theta_logSD <- function(n,  meanlog = 0, sdlog = 1,seed = NULL, rank = NULL) {
  set.seed(seed)
  val <- rlnorm(n, meanlog, sdlog)  # generate values from Log-normal distribution
  logSD <- log(sqrt(val))
  
  if (!is.null(rank)) {
    logSD <- logSD[1:rank]
    return(logSD)
  } else {
    return(logSD)
  }
}



# extract otu table and metadata
otu_meta_fun <- function(dd){
  
  wide_dd <- dd %>%
    dplyr::select(subject, taxon, count, group) %>%
    spread(key = taxon, value = count) %>% 
    setNames(c("subject","group", paste0("taxon",1:ntaxa)))    
#  colnames(wide_dd) = c("subject","group", paste0("taxon",1:ntaxa))
  
  otu_table  =  wide_dd %>%
    dplyr::select(paste0("taxon",1:ntaxa))
  
  metadata   =  wide_dd %>% 
    dplyr::select(c("subject","group"))
  
  otu_count  =   t(otu_table)
  dds        =   DESeqDataSetFromMatrix(otu_count,metadata, ~group)
  dds        =   DESeq(dds,sfType ="poscounts",minReplicatesForReplace=Inf) 
  normalizer =   sizeFactors(dds) # one for each subject
  
  
  normalise_dd  =    data.frame(normalizer, subject = wide_dd$subject)
  dd            =    left_join(dd,normalise_dd, by ="subject")  
  dd
  
  #list(metadata = metadata, otu_table = otu_table,wide_dd=wide_dd)
}

 
