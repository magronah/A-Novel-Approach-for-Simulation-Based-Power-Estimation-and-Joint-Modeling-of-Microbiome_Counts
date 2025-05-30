combine_model_component <- function(dd, component, models = c("RR", "US", "NB", "ZINB")) {
  do.call(rbind, lapply(models, function(model) {
    df <- dd$error[[model]][[component]]
    df$model <- model
    df
  }))
}


make_plot <- function(data, x_var, y_var,x_label, y_label, title, 
                      font_size = n, theme_type = "bw", 
                      show_true = FALSE, true_var = "true_param", 
                      use_custom_theme = TRUE) {
  
  # Start the base plot
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], 
                        color = model, group = model)) +
    geom_point(size = 1.5, alpha = 0.5) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = oka_col) +
    labs(title = title, x = x_label, y = y_label, color = "Model") +
    theme(axis.text.x = element_blank())
  
  # Optional true effect point overlay (for average estimate plot)
  if (show_true) {
    p <- p + geom_point(aes(y = .data[[true_var]], color = "true effect"), 
                        size = 3, shape = 17)
  }
  
  # Theme handling
  if (use_custom_theme) {
    p <- p + custom_theme(font_size)
  } else if (theme_type == "bw") {
    p <- p + theme_bw(base_size = font_size)
  } else if (theme_type == "minimal") {
    p <- p + theme_minimal(base_size = font_size)
  }
  
  # Final text and layout styling
  p <- p + theme(
    plot.title = element_text(hjust = 0.5, size = font_size),
    panel.grid = element_blank(),
    text = element_text(size = font_size, family = "Roboto")
  )
  
  return(p)
}


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

load_data <- function(path, alpha = 0.05) {
  
  true_param =  readRDS(paste0(path,"true_param.rds"))
  ####################################################
  RR     =  readRDS(paste0(path, "rr.rds"))
  US     =  readRDS(paste0(path, "us.rds"))
  nb     =  readRDS(paste0(path, "nb.rds"))
  znb    =  readRDS(paste0(path, "zinb.rds"))
  
  NB     =  nb$est
  ZINB   =  znb$est
  
  colnames(NB)    =  paste0("sim", 1:ncol(NB))
  colnames(ZINB)  =  paste0("sim", 1:ncol(ZINB))
  ####################################################
  common_sim <- Reduce(intersect, 
                            list(colnames(RR), colnames(US), 
                                colnames(NB), colnames(ZINB)))
  
  RR_s    =  RR  %>% select(all_of(common_sim))
  US_s    =  US  %>% select(all_of(common_sim))
  NB_s    =  NB  %>% select(all_of(common_sim))
  ZINB_s  =  ZINB  %>% select(all_of(common_sim))

  common_taxa <- Reduce(intersect, 
                       list(rownames(RR_s), rownames(US_s), 
                            rownames(NB_s), rownames(ZINB_s)))
  
  RR    =  RR_s[rownames(RR_s) %in% common_taxa,]
  US    =  US_s[rownames(US_s) %in% common_taxa,]
  NB    =  NB_s[rownames(NB_s) %in% common_taxa,]
  ZINB  =  ZINB_s[rownames(ZINB_s) %in% common_taxa,]
  true_param  =   true_param[true_param$param_name %in% common_taxa,]
  ####################################################
  rrl     =   dd_long(RR,  true_param,label="RR")
  usl     =   dd_long(US,  true_param,label="US")
  nbmml   =   dd_long(NB,  true_param,label="NB")
  zinbmml =   dd_long(ZINB, true_param,label="ZINB")
  ####################################################
  est_dd = lst(RR, US, NB, ZINB)
  
  long_dd = list(RR = rrl, US = usl, 
                 NB  = nbmml, ZINB  = zinbmml)
  
  error = list(
    RR     =   error_cal(RR, true_param, model = "RR"),
    US     =   error_cal(US, true_param, model = "US"),
    NB     =   error_cal(NB, true_param, model = "NB"),
    ZINB    =   error_cal(ZINB, true_param, model = "ZINB")) 

  confint = list(
    RR      =   para_confint(rrl, true_param, alpha = alpha),
    US      =   para_confint(usl, true_param, alpha = alpha),
    NB      =   para_confint(nbmml, true_param, alpha = alpha),
    ZINB     =   para_confint(zinbmml, true_param, alpha = alpha)) 
  
  lst(true_param,est_dd,long_dd,confint, error)
}

error_cal <- function(model_est_dd, true_param_dd, model,
                      lwr  =  0.025, upr = 0.975,  scale = 1.96, 
                      sdrr = TRUE) {
  
  est_dd  =   model_est_dd %>%
    rownames_to_column(var = "param_name")
  
  merge_ddd     =   left_join(est_dd, true_param_dd, by = "param_name")  
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

# zinbmm_confint <- function(mod, 
#                                     group_label = "pregnant:GA_Days",
#                                     conf_level = 0.95, 
#                                     mean_count = NULL) {
#   # Extract standard errors, estimates, and p-values
#   sd_err <- as.numeric(unlist(lapply(mod$fit, function(x) {
#     summary(x)$tTable[group_label, "Std.Error"]
#   })))
#   
#   est <- as.numeric(unlist(lapply(mod$fit, function(x) {
#     summary(x)$tTable[group_label, "Value"]
#   })))
#   
#   pvalue <- as.numeric(unlist(lapply(mod$fit, function(x) {
#     summary(x)$tTable[group_label, "p-value"]
#   })))
#   
#   z_score <- qnorm(conf_level + (1 - conf_level) / 2)
#   
#   # Build result data frame
#   dd <- data.frame(
#     est_param = est,
#     sd_err = sd_err,
#     lwr = est - z_score * sd_err,
#     upr = est + z_score * sd_err,
#     width = 2 * z_score * sd_err,
#     pvalue = pvalue
#   )
#   
#   len <- length(group_label)
#   
#   if (!is.null(mean_count)) {
#     if (len > 1) {
#       dd$group_label <- rep(group_label, length = len * length(mean_count))
#       mean_cnt <- rep(mean_count, each = len)
#       dd$mean_count <- as.numeric(mean_cnt)
#       dd$param_name <- names(mean_cnt)
#     } else {
#       dd$mean_count <- as.numeric(mean_count)
#       dd$param_name <- names(mean_count)
#     }
#   } else {
#     if (len > 1) {
#       dd$group_label <- rep(group_label, length = len * length(mod$fit))
#     }
#     dd$param_name <- mod$responses
#   }
#   
#   return(dd)
# }

zinbmm_confint = function(mod, mean_count,
                          group_label = "pregnant:GA_Days",
                          conf_level = .95){

  sd_err =  as.numeric(unlist(lapply(mod$fit,function(x)
  {summary(x)$tTable[group_label, "Std.Error"]})))

  est    =   as.numeric(unlist(lapply(mod$fit,function(x)
  {summary(x)$tTable[group_label, "Value"]})))

  pvalue =  as.numeric(unlist(lapply(mod$fit,function(x)
  {summary(x)$tTable[group_label, "p-value"]})))

  dd         =   data.frame(est_param = est)
  z_score    =   qnorm(conf_level + (1 - conf_level)/2)

  dd$sd_err  =   sd_err
  dd$lwr     =   est  -  z_score*sd_err
  dd$upr     =   est  +  z_score*sd_err
  dd$width   =   dd$upr  - dd$lwr

  dd$pvalue      =   pvalue
  len  =   length(group_label)

  if(len > 1){
    dd$group_label   =   rep(group_label, length = len*length(mean_count))
      mean_cnt       =   rep(mean_count, each = len)
      dd$mean_count  =   as.numeric(mean_cnt)
      dd$param_name  =   names(mean_cnt)
      return(dd)
  }else{
    dd$mean_count  =   as.numeric(mean_count)
    dd$param_name  =   names(mean_count)
    return(dd)
  }
}

my_aicc_fun  <- function(mod){
  
  loglik      =   as.numeric(logLik(mod))
  num_params  =   attr(logLik(mod), "df")  
  correction  =   2*num_params*(num_params+1)/(nobs(mod) - num_params -  1)
  
  if(is.na(loglik)){
    loglik      =    as.numeric(mod$obj$fn())
    aic         =    2*(loglik) +  2*(num_params)
    aicc        =    aic   +  correction
    return(aicc)
  }else{
    aic         =   -2*(loglik) +  2*(num_params)
    correction  =   2*num_params*(num_params+1)/(nobs(mod) - num_params -  1)
    aicc        =   aic   +  correction
    return(aicc)
  }
}

#filter_complete_separation = function(dd, thresh=100){
  
#  param_names_to_remove <- unique(dd$param_name[dd$width > thresh])
#  dd[!dd$param_name %in% param_names_to_remove, ]
#}

custom_theme <- function(n) {
  theme_bw(base_size = n) +
    theme(
      plot.title = element_text(size = n, family = "Roboto", hjust = 0.5),
      text = element_text(size = n, family = "Roboto"),
      axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
      axis.text.y = element_text(family = "Roboto", size = n, color = "black")
    )
}


#' @param mod 
#' @param ntaxa 
#' @param conf_level 
#'
#' @return
#' @export
#' 
#' @examples
wald_confint2 = function(mod, conf_level = .95, 
                        in_sd = 1, ntaxa,
                        mean_count, mod_name,
                        path){
  
  ref    =  ranef(mod, condVar = FALSE)
  ########################################################
  jp <- TMB::sdreport(mod$obj, getJointPrecision = TRUE)$jointPrecision
  saveRDS(jp, file = paste0(path,"sdreport_",mod_name,".rds"))
  ########################################################
  #extract precision matrix corresponding to the bs only
  b_inds <- which(rownames(jp) == "b")
  jpb0   <- jp[b_inds, b_inds]
  ########################################################
  # form <- ~1 + (group*time|taxon) +  (1|nugget) +
  #            (1|subject:taxon) + (taxon + 0 | subject:time) 
  #extract b's precision matrix corresponding to (group * time | taxon) term
      n1     <-  length(names(ref$cond$taxon))
  b_sub_inds <-  1:(n1*ntaxa)
  jpb        <-  jpb0[b_sub_inds, b_sub_inds]
  inv_mat    <-  solve(jpb)

  if(any(diag(inv_mat) < 0)){
    pre_nearPD <-  as.matrix(nearPD(jpb)$mat)
    se_vec     <-  sqrt(diag(solve(pre_nearPD)))
    print("precision matrix was converted to be positive
          definate using near precison")
  }else{
    se_vec     <- sqrt(diag(inv_mat))
  }
  #########################################################
  saveRDS(se_vec, file = paste0(path,"full_sdr_",mod_name,".rds"))
  #########################################################
     # n2     <-  sum(grepl(".*:time", names(ref$cond$taxon)))
  grp_label <-   grep("[:].*(time|days|day|week|weeks|times)", names(ref$cond$taxon), 
                      ignore.case = TRUE, 
                      value = TRUE)
  n2   <-    length(grp_label)
  ############################################
  p         <-   1:(n1*ntaxa)
  nselect   <-   split(p, ceiling(seq_along(p) / n1))  
  subselect <-   lapply(nselect, function(x) tail(x, n2))  
  grp_ind   <-   unique(unlist(subselect))
  #########################################################
  sd_err    <-   se_vec[grp_ind]
  saveRDS(sd_err, file = paste0(path,"sdr_",mod_name,".rds"))
  #########################################################
  # Calculate z-score for the desired confidence level
  z_score    =    qnorm(conf_level + (1 - conf_level)/2)
  
  full_est  <-   mod$fit$parfull
  est       <-   full_est[names(full_est) == "b"][grp_ind]
  stopifnot(est == unlist(t(ref$cond$taxon))[grp_ind])
  #########################################################
  # Calculate p-values
  z_stat   =  est / sd_err
  p_values =  2 * (1 - pnorm(abs(z_stat)))  # Two-tailed p-value
  #########################################################
  # Calculate confidence intervals
  dd   <-  data.frame(est_param =   est,
                      lwr       =   est - z_score*in_sd*sd_err,
                      upr       =   est + z_score*in_sd*sd_err)
  #########################################################
  dd$pvalue  =   p_values
  dd$width   =   dd$upr  - dd$lwr
  dd$sd_err  =   sd_err
  
  if(n2 > 1){
    group_label     =  names(ref$cond$taxon)[
                        grepl(".*:time", names(ref$cond$taxon))]
    dd$group_label  =   rep(group_label, length = n2*length(mean_count))
    mean_cnt        =   rep(mean_count, each = n2)
    dd$mean_count   =   as.numeric(mean_cnt)
    dd$param_name   =   names(mean_cnt)
    return(dd)
  }else{
    dd$mean_count  =   as.numeric(mean_count)
    dd$param_name  =   names(mean_count)
    return(dd)
  }
}


wald_confint3 = function(mod, conf_level = .95, 
                         in_sd = 1, ntaxa,
                         true_param,
                         mean_count, mod_name = NULL,
                         path = NULL){
  
  ref    =  ranef(mod, condVar = FALSE)
  ########################################################
  jp <- TMB::sdreport(mod$obj, getJointPrecision = TRUE)$jointPrecision
  #saveRDS(jp, file = paste0(path,"sdreport_",mod_name,".rds"))
  ########################################################
  #extract precision matrix corresponding to the bs only
  b_inds <- which(rownames(jp) == "b")
  jpb0   <- jp[b_inds, b_inds]
  ########################################################
  # form <- ~1 + (group*time|taxon) +  (1|nugget) +
  #            (1|subject:taxon) + (taxon + 0 | subject:time) 
  #extract b's precision matrix corresponding to (group * time | taxon) term
  n1     <-  length(names(ref$cond$taxon))
  b_sub_inds <-  1:(n1*ntaxa)
  jpb        <-  jpb0[b_sub_inds, b_sub_inds]
  inv_mat    <-  solve(jpb)
  
  if(any(diag(inv_mat) < 0)){
    pre_nearPD <-  as.matrix(nearPD(jpb)$mat)
    se_vec     <-  sqrt(diag(solve(pre_nearPD)))
    print("precision matrix was converted to be positive
          definate using near precison")
  }else{
    se_vec     <- sqrt(diag(inv_mat))
  }
  #########################################################
  #saveRDS(se_vec, file = paste0(path,"full_sdr_",mod_name,".rds"))
  #########################################################
  # n2     <-  sum(grepl(".*:time", names(ref$cond$taxon)))
  grp_label <-   grep("[:].*(time|days|day|week|weeks|times)", names(ref$cond$taxon), 
                      ignore.case = TRUE, 
                      value = TRUE)
  n2   <-    length(grp_label)
  ############################################
  p         <-   1:(n1*ntaxa)
  nselect   <-   split(p, ceiling(seq_along(p) / n1))  
  subselect <-   lapply(nselect, function(x) tail(x, n2))  
  grp_ind   <-   unique(unlist(subselect))
  #########################################################
  sd_err    <-   se_vec[grp_ind]
  #saveRDS(sd_err, file = paste0(path,"sdr_",mod_name,".rds"))
  #########################################################
  # Calculate z-score for the desired confidence level
  z_score    =    qnorm(conf_level + (1 - conf_level)/2)
  
  full_est  <-   mod$fit$parfull
  est       <-   full_est[names(full_est) == "b"][grp_ind]
  stopifnot(est == unlist(t(ref$cond$taxon))[grp_ind])
  #########################################################
  # Calculate p-values
  z_stat   =  est / sd_err
  p_values =  2 * (1 - pnorm(abs(z_stat)))  # Two-tailed p-value
  #########################################################
  # Calculate confidence intervals
  dd   <-  data.frame(est_param =   est,
                      lwr       =   est - z_score*in_sd*sd_err,
                      upr       =   est + z_score*in_sd*sd_err)
  #########################################################
  dd$pvalue  =   p_values
  dd$width   =   dd$upr  - dd$lwr
  
  if(n2 > 1){
    group_label     =  names(ref$cond$taxon)[
      grepl(".*:time", names(ref$cond$taxon))]
    dd$group_label  =   rep(group_label, length = n2*length(mean_count))
    mean_cnt        =   rep(mean_count, each = n2)
    dd$mean_count   =   as.numeric(mean_cnt)
    dd$param_name   =   names(mean_cnt)
    dd              =  left_join(dd,true_param, by = "param_name")
    return(dd)
  }else{
    dd$mean_count  =   as.numeric(mean_count)
    dd$param_name  =   names(mean_count)
    dd             =  left_join(dd,true_param, by = "param_name")
    return(dd)
  }
}

nb_znb_confint = function(mod, mean_count,
                          group_label = "bmi_groupsevereobese:time",
                          conf_level = .95,
                          thresh = 500){
  
  Sd_err =  as.numeric(unlist(lapply(mod$fit,function(x) 
  {summary(x)$tTable[group_label, "Std.Error"]})))
  
  Est    =   as.numeric(unlist(lapply(mod$fit,function(x) 
  {summary(x)$tTable[group_label, "Value"]})))
  
  Pvalue =  as.numeric(unlist(lapply(mod$fit,function(x) 
  {summary(x)$tTable[group_label, "p-value"]})))
  
  indx     =  which(Sd_err > thresh)
  
  if(length(Sd_err) > 0){
    sd_err   =  Sd_err[-indx]
    est      =  Est[-indx]
    pvalue   =  Pvalue[-indx]
    mean_count= mean_count[-indx]
  }else{
    sd_err   =  Sd_err
    est      =  Est
    pvalue   =  Pvalue
    mean_count= mean_count
  }
  
  dd         =   data.frame(est_param = est)
  z_score    =   qnorm(conf_level + (1 - conf_level)/2)
  
  dd$lwr     =   est  -  z_score*sd_err
  dd$upr     =   est  +  z_score*sd_err
  dd$width   =   dd$upr  - dd$lwr
  
  dd$pvalue      =   pvalue  
  dd$mean_count  =  mean_count 
  dd$param_name  =  names(mean_count)
  dd
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


filter_complete_separation = function(dd, thresh = 10){
  df  = dd  %>% filter(abs(est_param) < thresh)
  taxa_exclude =  df$param_name
  dd   %>%  filter(param_name %in% taxa_exclude)
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
                        path){
  
  ref    =  ranef(mod, condVar = FALSE)
  n1     =  length(names(ref$cond$taxon))
  n2     =  sum(grepl(".*:time", names(ref$cond$taxon)))
  ############################################
      p     <-  1:(n1*ntaxa)
   nselect  <- split(p, ceiling(seq_along(p) / n1))  
  subselect <- lapply(nselect, function(x) tail(x, n2))  
  grp_ind   <- unique(unlist(subselect))
  ########################################################
  ss <- TMB::sdreport(mod$obj, getJointPrecision = TRUE)
  saveRDS(ss, file = paste0(path,"sdreport_",mod_name,".rds"))
  ss  <- readRDS(paste0(path,"sdreport_",mod_name,".rds"))
  
  # browser()
  #########################################################
  ss$jointPrecision <-  as(ss$jointPrecision, "sparseMatrix")
  chol_decomp <- Cholesky(ss$jointPrecision)
  inverse_mat <- solve(chol_decomp, Diagonal(nrow(ss$jointPrecision)))
  
  #inverse_mat       <-  solve(ss$jointPrecision)
  
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
  
  saveRDS(se_vec, file = paste0(path,"full_sdr_",mod_name,".rds"))
  # ########################################################
  full_est     =   mod$fit$parfull
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

load_models <- function(path, filenames) {
  file_paths <- paste0(path, filenames)
  mod_list  <- lapply(file_paths, readRDS)
}

otu_meta_lst_fun = function(res_dd, ntime){
  res_dd2   = list()
  for(i in 1:length(res_dd)){
    dd  =  res_dd[[i]]
    
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
  
 # ntime           =   length(unique(met_data$time))
  names(res_dd2)  =   paste0("time", 1:ntime)
  res_dd2
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



get_corr <- function(ntaxa, nsubject, sparse_thresh = 0.05, seed = NULL) {
  
  set.seed(seed)
  synthetic_data   <-   huge.generator(n = nsubject, d = ntaxa, graph = "scale-free")
  covariance_matrix  <- synthetic_data$sigma
  correlation_matrix <- cov2cor(covariance_matrix)
  
  correlation_matrix[abs(correlation_matrix) < sparse_thresh] <- 0
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


get_theta_logSD <- function(n,  meanlog = 0, sdlog = 1,seed = NULL, rank = NULL) {
  set.seed(seed)
  val <- rlnorm(n, meanlog, sdlog)  # Log-normal distribution
  logSD <- log(sqrt(val))
  
  if (!is.null(rank)) {
    logSD <- logSD[1:rank]
    return(logSD)
  } else {
    return(logSD)
  }
}


otu_meta_fun = function(meta_dd){
  
  count_table=   list()
  met_dd     =    meta_dd
  tt         =    unique(met_dd$time)
  
  for(i in 1:length(tt)){
    
    dd <- met_dd %>% 
      subset(time == tt[i])
    
    otu_table <- dd %>%
      dplyr::select(subject, taxon, count, group)
    
    otu_table <- spread(otu_table, key = taxon, value = count) 
    colnames(otu_table) <- c("subject","group", 
                             paste0("taxon",1:(ncol(otu_table)-2)))
    
    countdata     =    otu_table %>%
      dplyr::select(-c(subject,group))
    
    #rownames(countdata) = (meta_dd$subject)
    count_table[[i]] =  countdata[, colSums(countdata != 0) > 0]
  }
  
  names(count_table)  =   paste0("time", tt)
  
  list(count_table = count_table, 
       meta_data  =  met_dd)
  
}

# # Example usage
# logSD <- get_theta_logSD(n = 100)
# print(logSD)

# get_theta_logSD <- function(n, seed = NULL, rank = NULL, prob = 0.5) {
#   set.seed(seed)
#   val    =   rgeom(n, prob=prob) + 0.1
#   logSD  =   log(sqrt(val))
#   
#   if(!is.null(rank)){
#     logSD     =   logSD[1:rank]
#     return(logSD)
#   }
#   else{return(logSD)}
# }

#' @param model fitted TMB model
#' @param data original data
#' @param epsilon perturbation size
#' @param inds indices of observations to perturb
#' @param fit_method "hack" = quick, without finalizing model object; "update" = slow, complete model fit + update
#' @param pred_method "hack" = based on linear model; "predict" = from fitted model object
#' @param scale data ("response") or linear predictor ("link") scale
#' @param progress progress bar?
#' @param opt_args additional arguments for optimizer
#' @param return_delta for diagnostics/debugging: return unscaled delta rather than delta/eps?
leverage_brute <- function(model, data, epsilon = 1e-3, inds = seq(nrow(data)),
                           fit_method = c("hack", "update"),
                           pred_method = c("hack", "predict"),
                           scale = c("response", "link"),
                           progress = FALSE,
                           opt_args = list(),
                           return_delta = FALSE
) {
  
  scale <- match.arg(scale)
  fit_method <- match.arg(fit_method)
  pred_method <- match.arg(pred_method)
   
  n <- length(inds)
  
  if (progress) pb <- txtProgressBar(max = n, style = 3)
  ## for now, compute all leverages on link scale
  y_pred <- predict(model, type = "link")
  leverage <- rep(NA_real_, n)
  yname  <- as.character(formula(model)[[2]])
  
  ## extract parameters so we can start from best vals
  p0 <- with(model$obj$env, parList(last.par.best[-random]))
  p0 <- p0[lengths(p0) > 0]
  p0 <- p0[setdiff(names(p0), "b")]  ## drop 'b' parameters
  
  X <- getME(model, "X")
  Z <- getME(model, "Z")
  
  for (j in seq_along(inds)) {
    if (progress) setTxtProgressBar(pb, j)
    i <- inds[j]
    data_perturb <- data
    data_perturb[[yname]][i] <-  data[[yname]][i] + epsilon
    if (fit_method == "hack") {
      ## quick/hacked new fit
      ## don't want to see warnings about non-integer counts
      suppressWarnings(
        newfit0 <- update(model, start = p0, data = data_perturb, verbose = FALSE, doFit = FALSE)
      )
      system.time(newfit1 <- fitTMB(newfit0, doOptim = FALSE)) ## 1 second
      system.time(newfit2 <- with(newfit1,
                                  do.call(nlminb,
                                          c(list(start=par, objective=fn, gradient=gr),
                                            opt_args)))
      )
    } else {
      ## full new fit
      suppressWarnings(
        newfit0 <- update(model, start = p0, data = data_perturb, verbose = FALSE)
      )
      newfit1 <- newfit0$obj
    }
    pp <- with(newfit1$env, parList(last.par.best[-random]))
    pp$b  <-  newfit1$report()$b
    if (pred_method == "hack") {
      y_pred_pert <- drop(X[i,] %*% pp[["beta"]] + Z[i,] %*% pp[["b"]])
    } else {
      if (fit_method == "hack") stop("can't do regular pred with hacked fit")
      y_pred_pert <- predict(newfit0, type = "link")[i]
    }
    if (scale == "response") {
      linkinv <- family(model)$linkinv
      y_pred_pert <- linkinv(y_pred_pert)
      y_pred[i] <- linkinv(y_pred[i])
    }
    leverage[j] <- (y_pred_pert - y_pred[i])
    if (!return_delta) leverage[j] <-  leverage[j] / epsilon
  }
  if (progress) close(pb)
  return(leverage)
}



#' @param model fitted TMB model
#' @param data original data
#' @param epsilon perturbation size
#' @param inds indices of observations to perturb
#' @param fit_method "hack" = quick, without finalizing model object; "update" = slow, complete model fit + update
#' @param pred_method "hack" = based on linear model; "predict" = from fitted model object
#' @param scale data ("response") or linear predictor ("link") scale
#' @param progress progress bar?
#' @param opt_args additional arguments for optimizer
#' @param return_delta for diagnostics/debugging: return unscaled delta rather than delta/eps?
leverage_brute_modified <- function(model, data, epsilon = 1e-3, inds = seq(nrow(data)),
                                    fit_method = c("hack", "update"),
                                    pred_method = c("hack", "predict"),
                                    scale = c("response", "link"),
                                    progress = FALSE,
                                    opt_args = list(),
                                    return_delta = FALSE
) {
  
  scale <- match.arg(scale)
  fit_method <- match.arg(fit_method)
  pred_method <- match.arg(pred_method)
  
  n <- length(inds)
  
  if (progress) pb <- txtProgressBar(max = n, style = 3)
  ## for now, compute all leverages on link scale
  y_pred <- predict(model, type = "response")
  leverage <- rep(NA_real_, n)
  yname  <- as.character(formula(model)[[2]])
  
  ## extract parameters so we can start from best vals
  p0 <- with(model$obj$env, parList(last.par.best[-random]))
  p0 <- p0[lengths(p0) > 0]
  p0 <- p0[setdiff(names(p0), "b")]  ## drop 'b' parameters
  
  X <- getME(model, "X")
  Z <- getME(model, "Z")
  
  for (j in seq_along(inds)) {
    if (progress) setTxtProgressBar(pb, j)
    i <- inds[j]
    data_perturb <- data
    data_perturb[[yname]][i] <-  data[[yname]][i] + epsilon
    if (fit_method == "hack") {
      ## quick/hacked new fit
      ## don't want to see warnings about non-integer counts
      suppressWarnings(
        newfit0 <- update(model, start = p0, data = data_perturb, verbose = FALSE, doFit = FALSE)
      )
      system.time(newfit1 <- fitTMB(newfit0, doOptim = FALSE)) ## 1 second
      system.time(newfit2 <- with(newfit1,
                                  do.call(nlminb,
                                          c(list(start=par, objective=fn, gradient=gr),
                                            opt_args)))
      )
    } else {
      ## full new fit
      suppressWarnings(
        newfit0 <- update(model, start = p0, data = data_perturb, verbose = FALSE)
      )
      newfit1 <- newfit0$obj
    }
    pp <- with(newfit1$env, parList(last.par.best[-random]))
    pp$b  <-  newfit1$report()$b
    if (pred_method == "hack") {
      y_pred_pert <- drop(X[i,] %*% pp[["beta"]] + Z[i,] %*% pp[["b"]])
    } else {
      if (fit_method == "hack") stop("can't do regular pred with hacked fit")
      y_pred_pert <- predict(newfit0, type = "link")[i]
    }
    if (scale == "response") {
      linkinv <- family(model)$linkinv
      y_pred_pert <- linkinv(y_pred_pert)
      #y_pred[i] <- linkinv(y_pred[i])
    }
    leverage[j] <- (y_pred_pert - y_pred[i])
    if (!return_delta) leverage[j] <-  leverage[j] / epsilon
  }
  if (progress) close(pb)
  return(leverage)
}



## need to modify src/Makevars in glmmTMB directory to contain this
## (no fopenmp!)
install_glmmTMB <- function(pkgdir, libdir, clean_src = TRUE) {
  ## save existing src/Makevars, overwrite with what we want
  flags <- c("PKG_CPPFLAGS = -DTMBAD_FRAMEWORK -DTMBAD_INDEX_TYPE=uint64_t -DTMB_MAX_ORDER=4",
             "## PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS)",
             "## PKG_CXXFLAGS=$(SHLIB_OPENMP_CXXFLAGS)")
  td <- tempdir()
  makevars <- file.path(pkgdir, "src", "Makevars")
  file.rename(makevars, file.path(td, "Makevars"))
  on.exit(file.rename(file.path(td, "Makevars"), makevars))
  unlink(makevars)
  writeLines(flags, makevars)
  if (clean_src) {
    unlink(list.files(file.path(pkgdir, "src"),
                      pattern="\\.(o|so)$",
                      full.names = TRUE))
  }
  if (!dir.exists(libdir)) dir.create(libdir)
  system(sprintf("R CMD INSTALL -l %s %s",
                 libdir, pkgdir))
}

## need to install this way (with location adjusted to your liking)
## R CMD INSTALL -l ~/students/agronah/reduced_rank_mm/glmmTMB_lev glmmTMB

## from https://github.com/glmmTMB/glmmTMB/blob/leverage/misc/leverage.R
## @param diag Get diagonal only?
leverage <- function(fm, diag=TRUE) {
  library(RTMB)
  library(Matrix)
  has.random <- any(fm$obj$env$lrandom())
  obj <- fm$obj
  ## We mess with these... (cleanup on exit!)
  restore.on.exit <- c("ADreport",
                       "parameters",
                       "data")
  oldvars <- sapply(restore.on.exit, get, envir=obj$env, simplify=FALSE)
  restore.oldvars <- function(){
    for(var in names(oldvars)) assign(var, oldvars[[var]], envir=obj$env)
  }
  on.exit({restore.oldvars(); obj$retape()})
  ## #################################################################
  ## Get derivatives of prediction
  ##
  ##    mu_hat( b_hat( theta_hat(yobs), yobs) , theta_hat(yobs) )
  ##
  ## wrt. yobs.
  ## Note the three partial derivative 'paths' to consider.
  ## We can get this derivative by
  ##  1. yobs -> theta_hat(yobs)
  ##  2. theta -> mu_hat( b_hat( theta, yobs) , theta ) [fixed yobs]
  ##  3. yobs -> mu_hat( b_hat( theta, yobs) , theta ) [fixed theta]
  ## #################################################################
  ##parhat <- obj$env$last.par.best
  parhat <- fm$fit$parfull
  pl <- obj$env$parList(par=parhat)
  yobs <- obj$env$data$yobs
  obj$env$parameters <- pl
  theta <- parhat[!obj$env$lrandom()]  ## ALL top-level parameters
  b <- parhat[obj$env$lrandom()]
  if (!is.null(obj$env$spHess)) {
    Hbb <- obj$env$spHess(parhat, random=TRUE) ## Needed later for RE models
  }
  ## #################################################################
  ## 1. Get partial derivatives of theta_hat wrt to yobs
  ## Note: length(yobs) much greater that length(theta)
  ##       ==> Reverse mode AD is suitable !
  ## #################################################################
  ## Move 'yobs' from data -> parameters (preserve C++ template order!)
  obj$env$parameters <- c(list(yobs = yobs), obj$env$parameters)
  obj$env$data$yobs <- NULL
  obj$retape()
  ## New objective parameters: (yobs, b, theta)
  nobs <- length(obj$env$parameters$yobs)
  nb <- length(obj$env$random)
  ntheta <- length(obj$env$par) - nobs - nb
  TMB::config(tmbad.atomic_sparse_log_determinant=0, DLL="RTMB") ## TMB FIXME
  F <- GetTape(obj)
  r <- obj$env$random ## Which are random
  p <- tail(1:(nobs+ntheta), ntheta) ## Which are parameters *after* removing random
  ThetaHat <- F$laplace(r)$newton(p)
  J <- ThetaHat$jacobian(ThetaHat$par())
  ## Extra stuff we need in (3)
  F. <- F$jacfun() ## (yobs, [b], theta) -> (yobs, [b], theta)
  F. <- MakeTape(function(y) F.( c(y, parhat) ) [r] , yobs) ## yobs -> b
  Hby <- F.$jacfun(sparse=TRUE)(yobs)
  ## #################################################################
  ## 2. Get partial derivatives of mu_hat wrt to theta for *fixed* yobs
  ## Note: length(mu) much greater that length(theta)
  ##       ==> Forward mode AD is suitable !
  ## #################################################################
  obj$env$data$yobs <- yobs
  obj$env$parameters$yobs <- NULL
  obj$retape()
  r <- obj$env$random ## Which are now random
  F <- GetTape(obj)
  Bhat <- F$newton(r) ## theta -> bhat
  obj$env$data$doPredict <- as.double(1) ## Enable prediction of 'mu'
  obj$env$data$whichPredict <- as.double(1:nobs)
  obj$env$ADreport <- TRUE ## Modify return value from Cpp
  obj$retape()
  F <- GetTape(obj) ## (b, theta) -> mu
  ## This doesn't work:
  ## MuHat <- MakeTape(function(theta)F(c(Bhat(theta), theta)), theta)
  MuHat <- MakeTape(function(theta) {
    par <- advector(nb + ntheta) ## glmmTMB mixes order of parameters and random effects...
    r <- obj$env$lrandom()
    par[r] <- Bhat(theta)
    par[!r] <- theta
    F(par)
  } , theta)
  ## 'Adjoint trick'
  T2 <- MakeTape(function(weight) {
    WT <- MakeTape(function(th) sum(MuHat(th) * weight), theta)
    WT$jacfun()(advector(theta))
  }, rep(1,nobs))
  J2 <- T2$jacobian(yobs)
  if (diag) {
    term1 <- colSums(J*J2)
  } else {
    term1 <- t(J) %*% J2
  }
  ## #################################################################
  ## 3. Get partial derivatives of mu_hat wrt yobs for fixed theta
  ## Note: Tricky!
  ##       
  ## #################################################################
  term2 <- 0
  if (has.random) {
    F2 <- MakeTape(function(b) {
      par <- advector(nb + ntheta) ## glmmTMB mixes order of parameters and random effects...
      r <- obj$env$lrandom()
      par[r] <- b
      par[!r] <- theta
      F(par)
    }, b) ## (b) -> mu
    Hmb <- F2$jacfun(sparse=TRUE)(b) ## sparse deriv mu wrt b
    ## Implicit function theorem gives final partial deriv matrix:
    ##   - Hby %*% solve(Hbb) %*% Hbm
    ## of which we need the diagonal.
    ## Because mu and yobs link to the same random effects, all required b-cliques are part of Hbb !
    ## It follows that we can replace solve(Hbb) by its subset iH !
    if (diag) {
      if (length(b) == 0) {
        iH <- solve(Hbb)
      } else {
        iH <- TMB:::solveSubset(Hbb)
      }
      term2 <- -colSums( Hby * (  iH %*% t(Hmb) ) )
    } else {
      term2 <- -t(Hby) %*% solve(Hbb, t(Hmb))
    }
  }
  term1 + term2
}

# peakRAM_testfun <- function(nsubj = 100, ntax = 100, d = 2,
#                             include_ttt = FALSE, seed = 101,
#                             vars_include = c("nsubj", "ntax", "d", "include_ttt")) {
#   
#   if (packageVersion("glmmTMB") < "1.1.11") stop("are you using the formula_env branch?")
#   ## must assume we are using the `formula_env` branch of glmmTMB!! hard to test though
#   
#   library(peakRAM)
#   
#   ##
#   set.seed(seed)
#   dd <- expand.grid(subject = factor(seq(nsubj)),
#                     taxon = factor(seq(ntax)))
#   dd$group <- factor(ifelse(as.numeric(dd$subject) < nsubj %/% 2, "a", "b"))
#   
#   if (include_ttt) {
#     form <- y  ~ 1 + us(1 + group | taxon) + rr(0 + taxon | subject, d)
#   } else {
#     form <- y ~ 1 + rr(0 + taxon | subject, d)
#   }
#   dd$y <- simulate_new( form[-2],
#                         family = nbinom2,
#                         newdata = dd,
#                         control = list(set_formula_env = FALSE),
#                         newparams = list(beta = 1,
#                                          betadisp = 1,
#                                          theta = rep(0.1, ntheta(nsubj, ntax, d, include_ttt)))
#   )[[1]]
#   ## have to run this without parallelization, since we had to turn off
#   ## OpenMP for leverage calculations (we could try to load the full
#   ## version of glmmTMB with autopar for fitting the model, then
#   ## detach and load the glmmTMB_lev
#   p1 <- peakRAM(
#     mod <- glmmTMB(form,
#                    family = nbinom2,
#                    data = dd)
#   )
#   tmpf <- function(x, task = "model_fit") {
#     names(x) <- c("task", "time_sec", "total_RAM_Mb", "peak_RAM_Mb")
#     x$task <- task
#     v <- mget(vars_include, inherits = TRUE)
#     x <- do.call(data.frame, c(v, list(x)))
#     return(x)
#   }
#   p2 <- peakRAM(leverage(mod))
#   rbind(tmpf(p1), tmpf(p2, task = "leverage"))
# }


# ntheta <- function(nsubj = 100, ntax = 100, d = 2, include_ttt = FALSE) {
#   rr_n <- ntax*d - choose(d,2)
#   rr_n + ifelse(include_ttt, 3, 0)
# }

