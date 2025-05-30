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
  
  count_table   = list()
  met_dd     =    meta_dd
  tt       =   unique(met_dd$time)
  
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


load_data <- function(subfolder) {
  list(
    true_param = readRDS(paste0(path, subfolder, "/true_param.rds")),
    rr_est     = readRDS(paste0(path, subfolder, "/parametric_rr.rds")),
    nbmm_est   = readRDS(paste0(path, subfolder, "/parametric_nbmm.rds")),
    zinbmm_est = readRDS(paste0(path, subfolder, "/parametric_zinbmm.rds")),
    deseq_est  = readRDS(paste0(path, subfolder, "/parametric_deseq.rds")),
    us_est     = readRDS(paste0(path, subfolder, "/parametric_us.rds"))
  )
}
##################################################################
custom_theme <- function(n) {
  theme_bw(base_size = n) +
    theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = n, family = "Roboto"),
      axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
      axis.text.y = element_text(family = "Roboto", size = n, color = "black")
    )
}
#################################################################
pval_cal_fun =  function(dd_est){
  data.frame(pvalue=apply(dd_est, 1, function(x){
    2*min(c(mean(x < 0), mean(x > 0)))})) 
  }
#################################################################  
  
##coverage 
coverage_cal = function(dd){
  df <- dd %>%
    mutate(
      coverage = if_else(true_param > lwr & true_param < upr, 1, 0)
    )
  df
}
##################################################################
contour_plot_fun <- function(combined_data, 
                             power_estimate, 
                             cont_breaks,cust_theme = 11,
                             title_name = NULL){
  
  combined_data$pvalue_reject <- factor(combined_data$pval_reject)
  
  gg_2dimc <- (ggplot(combined_data)
               + aes(log2(lmean_abund), abs_lfc)
               + rasterise(geom_point(aes(color = pvalue_reject), alpha = 0.5))
               + xlab(TeX("$\\log_2$(mean counts)")) 
               + ylab(TeX("$\\log_2$|(fold change)|")) 
               + scale_colour_manual(values = c("black", "red")) 
               + geom_contour(data = power_estimate,
                              aes(z=power),lwd=1,
                              breaks = cont_breaks)
               + geom_label_contour(data = power_estimate, 
                                    aes(z= power,label = sprintf("%.3f", after_stat(level))),
                                    breaks = cont_breaks
               ) + custom_theme(cust_theme) +
                 labs(color="significant")
               
  )
  if(!is.null(title_name)){
    gg_2dimc = gg_2dimc  + ggtitle(title_name) +   theme(plot.title = element_text(hjust = 0.5))
  }
  gg_2dimc
}



gam_fit <- function(pval,lfoldchange,lmean_abund,
                    grid_len = 100, alpha_level = 0.05){
  
  pval_reject  =   (!is.na(pval) & pval < alpha_level)
  
  comb      =   tibble(lmean_abund  =  lmean_abund,
                       abs_lfc      =  abs(lfoldchange),
                       pval_reject  =  as.numeric(pval_reject))
  
  #fit scams
  fit_2d       =    scam(pval_reject ~ s(lmean_abund, abs_lfc, bs="tedmi"),
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


#' Title
#'
#' @param form 
#' @param meta_dd 
#' @param true_pars 
#' @param nsim 
#'
#' @return
#' @export
#'
#' @examples
sim_effect <- function(form, meta_dd, true_pars, nsim=10){
  res  =  foreach(i =1:nsim, .combine = "cbind") %do% {
    pars <- simulate_new(RHSForm(form, as.form = TRUE), nsim = 1, seed = i,
                         newdata   = meta_dd,
                         newparams = true_pars,
                         family     =  nbinom2, 
                         ziformula = ~1,
                         return_val = "pars"
    )
    stopifnot(theta_true == pars[names(pars)=="theta"])
    
    b_true_all      =   pars[names(pars)=="b"]
    true_b_index    =   seq(2,(2*ntaxa),2)
    true_param      =   b_true_all[true_b_index]
    true_param
  }
  rownames(res)     =   paste0("taxon", 1:ntaxa)
  colnames(res)     =   paste0("sim", 1:nsim)
  res
}


#' Title
#'
#' @param effects_dd 
#' @param theta_true 
#' @param meta_dd 
#'
#' @return
#' @export
#'
#' @examples
sim_count  =  function(effects_dd,theta_true,meta_dd){
  
  nsim_ =  ncol(effects_dd) 
  
  res  =  foreach(i = 1:nsim_) %do% {
    pars <- simulate_new(RHSForm(form, as.form = TRUE), nsim = 1, seed = i,
                         newdata   = meta_dd,
                         newparams =  true_pars,
                         family     =  nbinom2, 
                         ziformula = ~1,
                         return_val = "pars"
    )
    stopifnot(theta_true == pars[names(pars)=="theta"])
    b_true_all      =   pars[names(pars)=="b"]
    true_b_index    =   seq(2,(2*ntaxa),2)
    true_param      =   b_true_all[true_b_index]
    
    stopifnot(true_param == effects_dd[,i])
    ############################################################
    meta_dd$count <- simulate_new(RHSForm(form, as.form = TRUE),
                          newdata   = meta_dd,
                          newparams = true_pars,
                          family = nbinom2,
                          ziformula = ~1,
                          seed = i)[[1]]
    meta_dd
  }
  names(res)     =   paste0("sim", 1:nsim_)
  res
  #res            =    as.data.frame(res)
  #df   =   cbind(meta_dd,res)
  #df
}


# extract otu table and metadata
otu_meta_fun <- function(dd){
  
  wide_dd <- dd %>%
    dplyr::select(subject, taxon, count, group) %>%
    spread(key = taxon, value = count) %>%
    setnames(c("subject","group", paste0("taxon",1:ntaxa))) 
  
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

met_data = function(ntaxa, nsubj){
  metadata <- data.frame(subject = factor(seq(nsubj)),
                         group = rep(c("control", "treat"), each = nsubj/2))
  df <- expand.grid(taxon = factor(seq(ntaxa)), subject = factor(seq(nsubj)))
  df <- merge(metadata, df, by = "subject")
  df
}


deseqfun <- function(countdata,met_data,alpha_level=0.1,ref_name="NT",
                     minReplicatesForReplace = Inf, 
                     cooksCutoff = FALSE,
                     independentFiltering = FALSE, 
                     shrinkage_method="normal"){
  
  #check otu table is in otu by samples format
  if(all((met_data$subject)==colnames(countdata)) == FALSE){
    countdata = t(countdata)
  }
  
  #remove samples with zeros for all taxa (if any such sample exist)
  keep <- (colSums(countdata) > 0)
  countdata = countdata[,keep]
  met_data= met_data[keep, ]
  
  # call deseq
  dds <- DESeqDataSetFromMatrix(countdata,met_data, ~group)
  dds$group <- relevel(dds$group, ref = ref_name)
  
  dds <- DESeq(dds,sfType ="poscounts",
               minReplicatesForReplace = minReplicatesForReplace) 
  
  res <- results(dds, cooksCutoff=cooksCutoff, 
                 independentFiltering=independentFiltering,
                 alpha = alpha_level)
  
  reslt <- lfcShrink(dds, res=res, coef=2, type=shrinkage_method)
  
  deseq_est = data.frame(reslt)
  deseq_est$dispersion = dispersions(dds)
  deseq_dd   =  deseq_est  %>% 
               rownames_to_column(var = "param_name")
  deseq_dd
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
    mutate(param_name = factor(param_name, 
                               levels = unique(param_name[order(true_param)])))
  
  ddd
}


deseq_wald_confint = function(deseq_est, true_param_dd,conf_level = .95,model = "deseq_wald"){
  foldchange =  2^deseq_est$log2FoldChange
  dd         =  data.frame(foldchange = foldchange)
  z_score    =    qnorm(conf_level + (1 - conf_level)/2)
  sd_err     =   deseq_est$lfcSE
  dd$lwr     =   foldchange  -  z_score*sd_err
  dd$upr       =   foldchange  +  z_score*sd_err
  dd$model  =  rep(model,nrow(dd)) 
  dd$param_name  =  deseq_est$param_name
  dd$true_param  =   true_param_dd$true_param
  dd
}




error_cal <- function(model_est_dd, true_param_dd, model = "rr") {
  
  est_dd  =   model_est_dd %>%
    rownames_to_column(var = "param_name")
  merge_ddd     =   left_join(true_param_dd, est_dd, by = "param_name")  
  merge_dd      =   merge_ddd  %>% dplyr::select(-param_name) 
  error     =   data.frame(t(apply(merge_dd, 1, 
                                   function(x){(x["true_param"] -  x)})))
  
  df1       =   data.frame(param_name  =   merge_ddd$param_name,
                          true_param  =   merge_ddd$true_param,
                          bias        =   rowMeans(error),
                          mse         =   rowMeans(error^2))
  
  df2     =    data.frame(average_value  =   mean(rowMeans(error)))
  df3     =    data.frame(average_value  =   mean(rowMeans(error^2)))
  
  df1$model   =  rep(model,nrow(df1)) 
  df2$model   =  rep(model,nrow(df2)) 
  df3$model   =  rep(model,nrow(df3)) 
  
  res     =   lst(full_summary_dd=df1, error, avg_bias = df2, avg_mse = df3)
  
  return(res)
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
wald_confint = function(mod, ntaxa, conf_level = .95, in_sd = 1){
  pred  =   predict(mod, type = "latent", se.fit = TRUE) 
  est   =   pred$fit
  sd    =   pred$se.fit
  
  z_score    =    qnorm(conf_level + (1 - conf_level)/2)
  dd_full    =    data.frame(est_param  =   est,
                       lwr       =   est - z_score*in_sd*sd,
                       upr       =   est + z_score*in_sd*sd)
  
  indx       =    seq(2, (2*ntaxa),2)
  dd         =    dd_full[indx, ]  
    
  dd$param_name =   factor(paste0("taxon",1:ntaxa),
                           levels = paste0("taxon",1:ntaxa))
  
  dd$model =   rep("wald",nrow(dd))

  dd
}
