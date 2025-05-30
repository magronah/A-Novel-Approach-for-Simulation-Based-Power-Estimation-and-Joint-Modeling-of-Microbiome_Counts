
#' @param mod 
#' @param ntaxa 
#' @param conf_level 
#'
#' @return
#' @export
#' 
#' @examples
wald_confint3 = function(mod, conf_level = .95, 
                         in_sd = 1, ntaxa,
                         mean_count, mod_name,
                         path,
                         grp_label =  "groupNT"){
  
  ref    =  ranef(mod, condVar = FALSE)
  ########################################################
  jp <- TMB::sdreport(mod$obj, getJointPrecision = TRUE)$jointPrecision
  saveRDS(jp, file = paste0(path,"sdreport_short",mod_name,".rds"))
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
  saveRDS(se_vec, file = paste0(path,"full_sdr_short",mod_name,".rds"))
  #########################################################
  # n2     <-  sum(grepl(".*:time", names(ref$cond$taxon)))
  # grep("[:].*(group)", names(ref$cond$taxon), 
  #                     ignore.case = TRUE, 
  #                     value = TRUE)
  n2   <-    length(grp_label)
  ############################################
  p         <-   1:(n1*ntaxa)
  nselect   <-   split(p, ceiling(seq_along(p) / n1))  
  subselect <-   lapply(nselect, function(x) tail(x, n2))  
  grp_ind   <-   unique(unlist(subselect))
  #########################################################
  sd_err    <-   se_vec[grp_ind]
  saveRDS(sd_err, file = paste0(path,"sdr_short",mod_name,".rds"))
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
      grepl(".*:group", names(ref$cond$taxon))]
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


## https://github.com/glmmTMB/glmmTMB/issues/1059#issuecomment-2187660510
library(nloptr)
library(glmmTMB)

nl_algs <- (nloptr.get.default.options()[1,"possible_values"]
            |> strsplit(",")
            |> unlist()
            |> trimws()
            |> grep(pattern="NLOPT_LD", value = TRUE)
            |> grep(pattern="AUGLAG", value = TRUE, invert = TRUE)
)

nloptwrap1 <- function(start, objective, gradient, control,
                       algorithm = "NLOPT_LD_LBFGS") {
  ## what are nlminb default tolerances?
  my_obj <- function(x) objective(x) ##
  my_grad <- function(x) gradient(x) ##
  fit <- nloptr(x0 = start, eval_f = my_obj, eval_grad_f = my_grad,
                opts = list(algorithm = algorithm, xtol_rel = 1e-8))
  return(fit)
}

optim_algs <- c("BFGS", "CG", "L-BFGS-B")
allFit.methods <- data.frame(
  optimizer = rep(c("nlminb", "optim", "nloptwrap1"),
                  times = c(1, length(optim_algs), length(nl_algs))),
  alg = c("default", optim_algs, nl_algs))

## FIXME: don't write over existing control ('smart_refit' branch has better
##        storage of control info)
## FIXME: don't refit with original alg?
## FIXME: allow for setting tolerances
## FIXME: factory etc. to store warnings
## FIXME: summary/tidy methods?
## FIXME: exclude terrible/sensitive optimizers
## FIXME: expand to use optimx methods?
## FIXME: we should only re-run the optimization, not update()! -> finalizeTMB()
##        if necessary do the modular steps again, ideally re-use info from fitted model
##        optional whether we should start from best fit or from default/specified starting values
##        optional whether to include original alg
allFit <- function(fit, methods = allFit.methods) {
  res <- vector("list", length = nrow(methods))
  names(res) <- paste(methods$optimizer,
                      methods$alg,
                      sep = ".")
  for (i in seq(nrow(methods))) {
    optimizer <- methods$optimizer[i]
    alg <- methods$alg[i]
    cat("** ", i, optimizer, alg, "\n")
    ctrl <- switch(optimizer,
                   nlminb = glmmTMBControl(optimizer=nlminb),
                   optim = glmmTMBControl(optimizer=optim,
                                          optArgs = list(method = alg)),
                   nloptwrap1 = glmmTMBControl(optimizer=nloptwrap1,
                                               optArgs = list(algorithm = alg)))
    time <- system.time(curfit <- try(update(fit, control = ctrl)))
    attr(curfit, "system.time") <- time
    res[[i]] <- curfit
  }
  res
}

# x1 <- update(mod, control = glmmTMBControl(optimizer = "nloptwrap1"))
# 
# modrr <- update(mod, control = glmmTMBControl(optimizer = nloptwrap1))
# x1 <- glmmTMB(formula(mod), data = getCall(mod)$data, 
#               family = getCall(mod)$family, 
#               control = glmmTMBControl(optimizer = "nloptwrap1"))

if (FALSE) {
  ## minimal example for testing allFit
  x0 <- glmmTMB(mpg ~ cyl, data = mtcars)
  x1 <- update(x0, control = glmmTMBControl(optimizer = nloptwrap1))
  x2 <- update(x0, control = glmmTMBControl(optimizer = nloptwrap2),
               optArgs = list(algorithm = nl_algs[1]))
  
  ## not identical (don't want them to be)
  fixef(x1)$cond-fixef(x2)$cond
}

load_models <- function(path, filenames) {
  file_paths <- paste0(path, filenames)
  mod_list  <- lapply(file_paths, readRDS)
}

filter_complete_separation = function(dd){
  df  = dd  %>% filter(abs(est_param) > 10)
  taxa_exclude =  df$param_name
  dd   %>%  filter(!param_name %in% taxa_exclude)
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
# filter_fun <- function(countdata, metadata,abund_thresh=10, sample_thresh=5){
#   # if( nrow(metadata) != ncol(countdata)){
#   #   countdata =  (countdata)
#   # }
#   keep <- rowSums(countdata >= abund_thresh) >= sample_thresh
#   countdata = countdata[keep,]
#   as.data.frame(t(countdata))
# }

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
  
 grp_ind   =   seq(2, 2*(ntaxa), 2)
  #pred    =   predict(mod, type = "latent", se.fit = TRUE) 
  #est     =   pred$fit[grp_ind]
  #sd_err  =   pred$se.fit[grp_ind]
  ########################################################
  ss <- TMB::sdreport(mod$obj, getJointPrecision = TRUE)
  saveRDS(ss, file = paste0(path,"sdreport_",mod_name,".rds"))
  #ss  <-  readRDS(paste0(path,"sdreport_",mod_name,".rds"))
  #########################################################
  #class(ss$jointPrecision)
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

#'
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
  
  
  dd       =   model.frame(mod)
  ntax     =   length(unique(dd$taxon))
  lev      =   length(levels(unique(dd$group)))
  entries  =   1:(lev*ntax)
  
  groups   =   split(entries, ceiling(seq_along(entries) / lev))
  grp_ind  =   as.numeric(unlist(lapply(groups, function(x) x[-1])))
  ########################################################
  ss <- TMB::sdreport(mod$obj, getJointPrecision = TRUE)
  saveRDS(ss, file = paste0(path,"sdreport_",mod_name,".rds"))
  #ss  <-  readRDS(paste0(path,"sdreport_",mod_name,".rds"))
  #ss$jointPrecision <-  as(ss$jointPrecision, "sparseMatrix")
  #inverse_mat       <-  solve(full_matrix)
  #se_vec      <-  sqrt(diag(inverse_mat))
  #se_vec0   <-    se_vec[names(se_vec) == "b"]#[entries]
  #########################################################
  full_matrix <-  (ss$jointPrecision)
  b_indices   <-  which(grepl("^b$", rownames(full_matrix), 
                              perl = TRUE))

  b_matrix    <-  full_matrix[b_indices,b_indices]
  b_subMat    <-  b_matrix[entries,entries] 
  inverse_mat <-  solve(b_subMat)
  se_vec      <-  sqrt(diag(b_subMat))  
  
  #plot(se_vec0,se_vec1)
  #abline(0,1)
  #########################################################
  saveRDS(se_vec, file = paste0(path,"full_sdr_",mod_name,".rds"))
  # ########################################################
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

##############################################
extract_grp_effect  =  function(dd_full,ntaxa){
  indx       =    seq(2, (2*ntaxa), 2)
  dd         =    dd_full[indx, ]
  dd$width    =    dd$upr  -  dd$lwr 
  # Add parameter names
  dd$param_name =   factor(paste0("taxon",1:ntaxa),
                           levels = paste0("taxon",1:ntaxa))
  dd
}

deseq_wald_confint = function(deseq_est,mean_count, conf_level = .95){
  foldchange =   log(2)*deseq_est$log2FoldChange
  dd         =   data.frame(est_param = foldchange)
  z_score    =   qnorm(conf_level + (1 - conf_level)/2)
  sd_err     =   log(2)*deseq_est$lfcSE
  dd$lwr     =   foldchange  -  z_score*sd_err
  dd$upr     =   foldchange  +  z_score*sd_err
  dd$width   =   dd$upr  - dd$lwr
  dd$param_name  =  deseq_est$param_name
  dd$pvalue      =   deseq_est$padj  
  dd$mean_count  =  as.numeric(mean_count)
  dd$param_name  =  names(mean_count)
  dd
}

#mod =  autism_models$Zinbmm

zinbmm_confint = function(mod, mean_count,
                                group_label = "groupNT",
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
  dd$mean_count  =  as.numeric(mean_count)
  dd$param_name  =  names(mean_count)
  dd
}

custom_theme <- function(n) {
  theme_bw(base_size = n) +
    theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = n, family = "Roboto"),
      axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
      axis.text.y = element_text(family = "Roboto", size = n, color = "black")
    )
}


filter_fun <- function(countdata, metadata,abund_thresh=5, sample_thresh=3){
  
  ## sanity check
  if(all((metadata$subject)==colnames(countdata)) == FALSE){
    countdata = t(countdata)
  }
  
  ##########################################
  # filter
  dds <- DESeqDataSetFromMatrix(countdata,metadata, ~group)
  keep <- rowSums(counts(dds) >= abund_thresh) >= sample_thresh
  
  dds=dds[keep,]
  data.frame(counts(dds))
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



# extract otu table and metadata
otu_meta_fun <- function(dd, ntaxa){
  wide_dd <- dd %>%
    dplyr::select(subject, taxon, count, group) %>%
    spread(key = taxon, value = count) %>%
    setNames(c("subject","group", paste0("taxon",1:ntaxa))) 
  
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


deseqfun <- function(countdata,met_data,alpha_level=0.1,ref_name="NT",
                     minReplicatesForReplace = Inf, 
                     cooksCutoff = FALSE,
                     independentFiltering = FALSE, 
                     do_shrinkage =  "yes", 
                     shrinkage_method="normal",
                     design = ~group, 
                     ntaxa){
  

  #check otu table is in otu by samples format
  if(nrow(countdata) != ntaxa){
    countdata = t(countdata)
  }
  
  #remove samples with zeros for all taxa (if any such sample exist)
  keep <- (colSums(countdata) > 0)
  countdata = countdata[,keep]
  met_data  = met_data[keep, ]
  
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
  
  list(result = deseq_dd, 
       data  =  list(countdata =  countdata,
                     meta_data = met_data),
       object = dds)
  
}
