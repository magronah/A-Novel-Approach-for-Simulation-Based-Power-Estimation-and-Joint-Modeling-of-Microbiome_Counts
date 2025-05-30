RRSim <- function(object, dflong, 
                  nsim, seed,
                  subject_label =  "subject", 
                  group_label   =  "group", 
                  taxa_label    =  "taxon",
                  form  =  as.formula(
                    paste0("count ~ 1 + us(1 + group|taxon) + rr(0 + taxon | subject,", 2, ")"))
){
  
  required_cols  = c(subject_label,group_label,taxa_label)
  stopifnot(all(required_cols %in% names(dflong)))
  
  countdata_list  =   metadata_list  =  list()
  par_list        =   object$obj$env$parList()
  
  true_pars0      =   getME(object, "theta")
  b_true_all      =   predict(object, type = "latent")
  par_list$b      =   predict(object, type = "latent")
  
  true_pars1      =  list(beta = par_list$beta, 
                          theta = true_pars0, b = b_true_all, 
                        betadisp = par_list$betadisp)
  
  sim <- simulate_new(RHSForm(form, as.form = TRUE), nsim = nsim, seed = seed,
                        newdata    =  dflong,
                        newparams  =  par_list,#true_pars1, #par_list
                        family     =  nbinom2, 
                        return_val =  "sim"
  )

  dd   =  dflong
  for(i in 1:nsim){
    dd$count   =   sim[[i]] 
    df    =   dd_wide_fun(dd, 
                          subject_label =  subject_label, 
                          group_label   =  group_label,
                          taxa_label    =  taxa_label)
    otu_count  =  df$countdata
    met_data   =  df$metadata
    colnames(otu_count) = met_data$subject
    countdata_list[[i]] =  otu_count
    metadata_list[[i]]   =  met_data
  }
  
  names(metadata_list)  =   names(countdata_list)  = paste0("sim",1:nsim)
  list(countdata_list  =  countdata_list, 
       metadata_list  =  metadata_list)
  
}

dd_wide_fun <- function(dd, 
                        subject_label =  "subject", 
                        group_label   =  "group", 
                        taxa_label    =  "taxon"){
  ntaxa    =  length(unique(dd[[taxa_label]]))
  
  wide_dd <- dd %>%
    dplyr::select(all_of(subject_label), taxon, count, all_of(group_label)) %>%
    spread(key = taxa_label, value = count) %>%
    select(subject, group, 
           mixedsort(names(.)[3:ncol(.)], decreasing = FALSE)) #%>% 
    #setNames(c(subject_label,group_label, paste0(taxa_label,1:ntaxa)))    
  
  otu_table  =  wide_dd %>%
    dplyr::select(-c(subject, group))
  
  # otu_table  =  wide_dd %>%
  #   dplyr::select(paste0("sp",1:ntaxa))
  
  metadata   =  wide_dd %>% 
    dplyr::select(c(subject_label,group_label))
  
  otu_count  =   t(otu_table)
  
  list(countdata  =  otu_count,
       metadata   =  metadata)
  
}


dd_long_fun = function(countdata, metadata, 
                       taxa_label = "sp",
                       subject_label = "subject"){
  
  # we want count data as subject by taxa
  if(ncol(countdata) == nrow(metadata)){
    countdata  =  as.data.frame(t(countdata))
  }else{
    countdata  =  data.frame(countdata)
  }
  df   =  countdata %>%
    rownames_to_column(subject_label)
  ddd =   pivot_longer(df,
                       cols = starts_with(taxa_label),
                       names_to  = taxa_label,
                       values_to = "count")
  
  names(ddd)  = c(subject_label,"taxon","count")
  ddd
}

otu_meta_fun <- function(dd, 
                         subject_label =  "subject", 
                         group_label   =  "group", 
                         taxa_label    =  "taxon"){
  
  
  ddf  =  dd_wide_fun(dd, 
                      subject_label =  subject_label, 
                      group_label   =  group_label,
                      taxa_label    =  taxa_label)
  
  otu_count  =   ddf$countdata
  metadata   =   ddf$metadata
  dds        =   DESeqDataSetFromMatrix(otu_count,metadata, design = 
                                          as.formula(paste("~", group_label)))
  dds        =   DESeq(dds,sfType ="poscounts",minReplicatesForReplace=Inf) 
  normalizer =   sizeFactors(dds) 
  
  
  normalise_dd          =    data.frame(metadata[[subject_label]])
  names(normalise_dd)  =    (subject_label)
  normalise_dd$normalizer  =  normalizer
  
  dd            =    left_join(dd,normalise_dd, by = subject_label) 
  dd
  
}


#' Title
#' @param mean 
#' @param dispersion 
#'
#' @return coefficient of variation
#' @examples
coeff_var <- function(mean, dispersion){
  var = mean*(1 + mean*dispersion)
  sqrt(var)/mean
}
###############################################################################
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

###############################################################################
custom_theme <- function(n) {
  theme_bw(base_size = n) +
    theme(
      plot.title = element_text(hjust = 0.5,size = n, family = "Roboto", color = "black"),
      text = element_text(size = n, family = "Roboto", color = "black"),
      axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
      axis.text.y = element_text(family = "Roboto", size = n, color = "black")
    )
}

###############################################################################
#' Title
#'
#' @param countdata:  otu table 
#' @param metadata:  
#' @param abund_thresh: minimum taxa abundance threshold
#' @param sample_thresh:    minimum sample threshold
#' @return countdata_filt: filtered otu table 
#'         control_count:  otu table for control group
#'         treat_count:  otu table for treatment group
#' @examples
filter_fun <- function(countdata, metadata,
                       abund_thresh=5,
                       sample_thresh=3,
                       subject_label  =  "Sample",
                       group_label = "Groups",
                       ntaxa){
  
  ## sanity check
  # we want count data as subject by taxa
  if(ncol(countdata) == ntaxa){
    countdata  =  as.data.frame(t(countdata))
  }else{
    countdata  =  data.frame(countdata)
  }
  
  ##########################################
  # filter
  dds <- DESeqDataSetFromMatrix(countdata,metadata, as.formula(paste("~", group_label)))
  keep <- rowSums(counts(dds) >= abund_thresh) >= sample_thresh
  
  dds=dds[keep,]
  countdata_filt= data.frame(counts(dds))
  ##########################################
  # extract control and treatment otu table 
  nt    =  metadata %>%  dplyr::filter(group_label == "NT")  %>% 
    dplyr::select(all_of(subject_label)) %>%  unlist()
  
  asd   =  metadata %>% dplyr::filter(group_label == "ASD")  %>% 
    dplyr::select(all_of(subject_label)) %>% unlist()
  
  control_count   =   countdata[,colnames(countdata) %in% nt]
  treat_count     =   countdata[,colnames(countdata) %in% asd]
  
  # sanity check
  stopifnot(colnames(control_count) ==  nt)
  stopifnot(colnames(treat_count)   ==  asd)
  
  list(countdata_filt =  countdata_filt, 
       control_count  =  control_count, 
       treat_count    =  treat_count)
}

########################################################
#' Title
#'
#' @param countdata 
#' @param metadata 
#' @param alpha_level 
#' @param ref_name: reference for fold change calculation 
#' @param minReplicatesForReplace: deseq's parameter to control minimum number of
#'     replicates needed for the replacement of outliers during dispersion estimation.
#' @param cooksCutoff: deseq's outlier removal or shrinkage 
#' @param independentFiltering: deseq's independent filtering
#' @param shrinkage_method:  deseq's shrinkage method
#'
#' @return
#' @export
#'
#' @examples
deseqfun <- function(countdata,metadata,alpha_level=0.1,ref_name="NT",
                     minReplicatesForReplace = Inf, 
                     cooksCutoff = TRUE,
                     independentFiltering = TRUE, 
                     shrinkage_method="normal"){
  
  #check otu table is in otu by samples format
  if(all((metadata$Samples)==colnames(countdata)) == FALSE){
    countdata = t(countdata)
  }
  
  #remove samples with zeros for all taxa (if any such sample exist)
  keep <- (colSums(countdata) > 0)
  countdata = countdata[,keep]
  metadata= metadata[keep, ]
  
  # call deseq
  dds <- DESeqDataSetFromMatrix(countdata,metadata, ~Groups)
  dds$Groups <- relevel(dds$Groups, ref = ref_name)
  
  dds <- DESeq(dds,sfType ="poscounts",
               minReplicatesForReplace = minReplicatesForReplace) 
  
  res <- results(dds, cooksCutoff=cooksCutoff, 
                 independentFiltering=independentFiltering,
                 alpha = alpha_level)
  
  reslt <- lfcShrink(dds, res=res, coef=2, type=shrinkage_method)
  
  deseq_est = data.frame(reslt)
  disp = dispersions(dds)
  
  logfoldchange = deseq_est$log2FoldChange
  names(logfoldchange) = rownames(deseq_est)
  
  normalised_count     =  counts(dds, normalized=TRUE)
  logmean = log2(rowMeans(normalised_count))
  
  list(logfoldchange = logfoldchange,
       logmean = logmean,
       dispersion = disp,
       deseq_estimate=deseq_est,
       normalised_count = normalised_count)
}

compare_dataset <- function(countdata_sim_list, countdata_obs, 
                            method = c("var", "mean", "compare_zeros"), 
                            n = 17,
                            lwd1 = 1, lwd2=2) {
  
  # Check if method is valid
  if (!(method %in% c("var", "mean", "compare_zeros"))) {
    stop("Invalid method. Please choose either 'var', 'mean', or 'compare_zeros'.")
  }
  
  # Create a list combining simulation data and observed data
  dlist <- list.append(countdata_sim_list, obs_filt = countdata_obs)
  
  # Initialize variable for the plot
  p <- NULL
  
  if (method == "var" || method == "mean") {
    # Choose function based on method
    if (method == "var") {
      calc_func <- stats::var
      xlab_text <- "$\\log_{10}$(variance of taxa)"
    } else {
      calc_func <- mean
      xlab_text <- "$\\log_{10}$(mean of taxa)"
    }
    
    # Calculate variance or mean for each dataset in the list
    vars <- dlist |> purrr::map(~ apply(., 1, calc_func)) |>
      purrr::map_dfr(~ tibble(var = .), .id = "type")
    
    # Separate observed and simulated data
    vars_obs <- vars[vars$type == "obs_filt", ]
    vars_sim <- vars[vars$type != "obs_filt", ]
    
    # Create ggplot for visualization
    p <- ggplot(vars_sim, aes(x = log10(var), colour = type)) + 
      geom_density(linewidth = lwd1) +
      geom_density(data = vars_obs, aes(x = log10(var)), 
                   colour = "black", linetype = "dashed", linewidth = lwd2) +
      custom_theme(n) +
      # theme_bw(base_size = n) +
      # theme(
      #   plot.title = element_text(hjust = 0.5, size = n, family = "Roboto"),
      #   text = element_text(size = n, family = "Roboto"),
      #   axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
      #   axis.text.y = element_text(family = "Roboto", size = n, color = "black"),
      #   legend.position = "bottom"
      # ) +
      xlab(TeX(xlab_text)) +
      labs(color = "simulation method")
  }
  
  if (method == "compare_zeros") {
    # Calculate proportion of zeros for each row in all datasets
    zero_props <- dlist |> purrr::map(~ apply(., 2, function(row) mean(row == 0))) |>
      purrr::map_dfr(~ tibble(zero_prop = .), .id = "type")
    
    # Separate observed and simulated data for zero proportions
    zero_props_obs <- zero_props[zero_props$type == "obs_filt", ]
    zero_props_sim <- zero_props[zero_props$type != "obs_filt", ]
    
    # Create ggplot for zero proportions comparison
    p <- ggplot(zero_props_sim, aes(x = zero_prop, colour = type)) + 
      geom_density(linewidth = lwd1) +
      geom_density(data = zero_props_obs, aes(x = zero_prop), 
                   colour = "black", linetype = "dashed", linewidth = lwd2) +
      custom_theme(n) +
      # theme_bw(base_size = n) +
      # theme(
      #   plot.title = element_text(hjust = 0.5, size = n, family = "Roboto"),
      #   text = element_text(size = n, family = "Roboto"),
      #   axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
      #   axis.text.y = element_text(family = "Roboto", size = n, color = "black"),
      #   legend.position = "bottom"
      # ) +
      xlab("Proportion of Zeros") +
      labs(color = "simulation method")
  }
  
  return(p)
}

################################################################################
compare_dataset <- function(countdata_sim_list, 
                            method = c("var", "mean", "compare_zeros"), 
                            n = 17,
                            lwd1 = 1, lwd2=2) {
  
  # Check if method is valid
  if (!(method %in% c("var", "mean", "compare_zeros"))) {
    stop("Invalid method. Please choose either 'var', 'mean', or 'compare_zeros'.")
  }
  
  # Create a list combining simulation data and observed data
  dlist <- countdata_sim_list
  
  # Initialize variable for the plot
  p <- NULL
  
  if (method == "var" || method == "mean") {
    # Choose function based on method
    if (method == "var") {
      calc_func <- stats::var
      xlab_text <- "$\\log_{10}$(variance of taxa)"
    } else {
      calc_func <- mean
      xlab_text <- "$\\log_{10}$(mean of taxa)"
    }
    
    # Calculate variance or mean for each dataset in the list
    vars <- dlist |> purrr::map(~ apply(., 1, calc_func)) |>
      purrr::map_dfr(~ tibble(var = .), .id = "type")
    
    # Separate observed and simulated data
    #vars_obs <- vars[vars$type == "obs_filt", ]
    vars_sim <- vars[vars$type != "obs_filt", ]
    
    # Create ggplot for visualization
    p <- ggplot(vars_sim, aes(x = log10(var), colour = type)) + 
      geom_density(linewidth = lwd1) +
      #geom_density(data = vars_obs, aes(x = log10(var)), 
       #            colour = "black", linetype = "dashed", linewidth = lwd2) +
      custom_theme(n) +
      # theme_bw(base_size = n) +
      # theme(
      #   plot.title = element_text(hjust = 0.5, size = n, family = "Roboto"),
      #   text = element_text(size = n, family = "Roboto"),
      #   axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
      #   axis.text.y = element_text(family = "Roboto", size = n, color = "black"),
      #   legend.position = "bottom"
      # ) +
      xlab(TeX(xlab_text)) +
      labs(color = "simulation method")
  }
  
  if (method == "compare_zeros") {
    # Calculate proportion of zeros for each row in all datasets
    zero_props <- dlist |> purrr::map(~ apply(., 2, function(row) mean(row == 0))) |>
      purrr::map_dfr(~ tibble(zero_prop = .), .id = "type")
    
    # Separate observed and simulated data for zero proportions
    zero_props_obs <- zero_props[zero_props$type == "obs_filt", ]
    zero_props_sim <- zero_props[zero_props$type != "obs_filt", ]
    
    # Create ggplot for zero proportions comparison
    p <- ggplot(zero_props_sim, aes(x = zero_prop, colour = type)) + 
      geom_density(linewidth = lwd1) +
      geom_density(data = zero_props_obs, aes(x = zero_prop), 
                   colour = "black", linetype = "dashed", linewidth = lwd2) +
      custom_theme(n) +
      # theme_bw(base_size = n) +
      # theme(
      #   plot.title = element_text(hjust = 0.5, size = n, family = "Roboto"),
      #   text = element_text(size = n, family = "Roboto"),
      #   axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
      #   axis.text.y = element_text(family = "Roboto", size = n, color = "black"),
      #   legend.position = "bottom"
      # ) +
      xlab("Proportion of Zeros") +
      labs(color = "simulation method")
  }
  
  return(p)
}
###############################################################################

compare_dataset1 <- function(countdata_sim_list,countdata_obs,method = c("var", "mean"),n= 17){
  
  # Check if method is either "var" or "mean"
  if (!(method %in% c("var", "mean"))) {
    stop("Invalid method. Please choose either 'var' or 'mean'.")
  }
  
  # Calculate variance or mean based on the chosen method
  if (method == "var") {
    calc_func <- stats::var
    xlab_text <- "$\\log_{10}$(variance of taxa)"
  } else {
    calc_func <- mean
    xlab_text <- "$\\log_{10}$(mean of taxa)"
  }
  
  # Create a list combining simulation data and observed data
  dlist <- list.append(countdata_sim_list, obs_filt = countdata_obs)

    # Calculate variance or mean for each dataset in the list
  vars <- dlist |> purrr::map(~ apply(., 1, calc_func)) |>
    purrr::map_dfr(~ tibble(var = .), .id = "type")
  
  # Separate observed and simulated data
  vars_obs <- vars[vars$type == "obs_filt", ]
  vars_sim <- vars[vars$type != "obs_filt", ]
  
  
  # Create ggplot for visualization
  p <- ggplot(vars_sim, aes(x = log10(var), colour = type)) + 
    geom_density(lwd = 0.8) +
    geom_density(data = vars_obs, aes(x = log10(var)), 
                 colour = "black", linetype = "dashed", lwd = 2) +
  #  scale_color_manual(values = okabe_ito_colors) + # Apply manually specified Okabe-Ito colors
    theme_bw(base_size = n) +
    theme(
      plot.title = element_text(hjust = 0.5, size = n, family = "Roboto"),
      text = element_text(size = n, family = "Roboto"),
      axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
      axis.text.y = element_text(family = "Roboto", size = n, color = "black"),
      legend.position = "bottom"  # Remove the legend
    ) +
    xlab(TeX(xlab_text)) +
    labs(color = "simulation method")
  
  return(p)
}

#zeros for later use
# hist(zeroes_prop_metaS)    
# # Proportion of zeros
# zeroes_prop_obs       =    rowMeans(countdata_obs == 0)
# zeroes_prop_HMP       =    rowMeans(HMP == 0)
# zeroes_prop_metaS     =    rowMeans(metaSPARSim == 0)
# zeroes_prop_scaled    =    rowMeans(scaled == 0) 
# zeroes_prop_unscaled  =    rowMeans(unscaled == 0) 

read_data <- function(dataset_list, extract_name){
  extract_data_list= list()
  for(n in 1:length(dataset_list)){
    dataset <- dataset_list[[n]]
    extract_data_list[[n]] <- dataset[[extract_name]]
  }
  names(extract_data_list) <- names(dataset_list)
  extract_data_list
}

####################################
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

####################################
# compare_dataset <- function(countdata_sim_list,countdata_obs){
#   dlist <- list.append(countdata_sim_list,obs_filt=countdata_obs)
#   vars <-(dlist |> purrr::map(~apply(., 1, stats::var)) |>
#             purrr::map_dfr(~ tibble(var = .),.id = "type")) 
#   
#   vars_obs <- vars[vars$type == "obs_filt",]
#   vars_sim <- vars[vars$type != "obs_filt",]
#   
#   p = ggplot(vars_sim, aes(x = log(var), colour = type)) + geom_density(lwd=1.5) +
#      geom_density(data=vars_obs, aes(x = log(var)),colour = "black", linetype = "solid",lwd=3) +
#     xlab(TeX("$\\log_{10}$(variance of taxa)"))  
# 
#   p
# }

# nc <- max(1, getOption("ncores", round(detectCores()/2)))
# cl <- makeCluster(ncore)

initial_cond_est <- function(p0,logmean, logfoldchange, 
                                       np, sd_ord, ncore, 
                                       minval = -5, maxval = 5, 
                                       NP = 800, itermax = 1500){
  
  cl <- makeCluster(ncore)
  fun_names <- c("polyfun", "dnormmix", "dnormmix0", "genmixpars")
  clusterExport(cl, c("logmean", "logfoldchange", fun_names))
  
  opt  <- DEoptim(fn = nllfun, lower = rep(minval, length(p0)), upper = rep(maxval, length(p0)),
                  vals = logfoldchange, logmean = logmean, np = np, sd_ord = sd_ord,
                  control = DEoptim.control(NP = NP, itermax = itermax,
                                            cluster = cl))
  
  opt$optim$bestmem
  
}


## general-purpose log-likelihood function
## vectorized sum(pars*x^i)
polyfun <- function(pars, x) {
  ord <- length(pars)
  xmat <- sapply(0:(ord-1), function(i) x^i)
  xmat |> sweep(MARGIN = 2, FUN = "*", pars) |> rowSums()
}

#' generate normal mixture parameters (prob vector, mean vector, sd vector
#' for a specified set of 'x' values (logmean)
#' @param x independent variable
#' @param pars parameter vector: first logit-probs (np-1), 
#'             then mean parameters (2 per component: intercepts, slopes), 
#'             then var parameters (varord + 1 per component: intercepts, 
#'             slopes, quad coeffs, etc.)
#' @param np number of components in mixture
#' @param sd_ord order of logsd model (2 = quadratic)
#' 

genmixpars <- function(x, pars, np = 2, sd_ord = 2){
  ## FIXME: complain if pars is wrong length
  expected_pars <- (np-1)+2*np+(sd_ord+1)*np
  if (length(pars) != expected_pars) {
    stop(sprintf("number of pars (%d) != expected (%d) (np = %d, sd_ord = %d)",
                 length(pars), expected_pars, np, sd_ord))
  }
  ## softmax function
  cp <- 1  ## parameter index
  probs <- exp(c(0,pars[1:(np-1)])) # first is zero
  probs <- probs/sum(probs)
  cp <- cp + (np-1)
  ## mu: always linear (2 params per component)
  ## intercepts first, then slopes
  mupars <- pars[cp:(cp + 2*np - 1)]
  mupars <- matrix(mupars, ncol = 2)
  muvals <- apply(mupars, 1, polyfun, x = x)
  cp <- cp + 2*np
  ## logsd: similar
  logsdpars <- pars[cp:(cp + (sd_ord+1)*np - 1)]
  logsdpars <- matrix(logsdpars, ncol = (sd_ord+1))
  logsdvals <- apply(logsdpars, 1, polyfun, x = x)
  list(probs = probs, muvals = muvals, sdvals = exp(logsdvals))
}

## general-purpose normal-mixture deviate generator: takes _matrices_
## of probabilities, means, sds
rnormmix0 <- function(n, probs, muvals, sdvals) {
  np <- length(probs)
  component <- sample(np, size = n, prob = probs, replace = TRUE)
  inds <- cbind(seq(n), component)
  rnorm(n, mean = muvals[inds], sd = sdvals[inds])
}

# logmean <- rnorm(10); par =p0
myrnormmix <- function(par, logmean, ...) {
  g0 <- genmixpars(logmean, par, ...)
  do.call(rnormmix0, c(list(n = length(logmean)), g0))
}

## takes pars as three vectors (prob, mean, sd), on constrained scale
## i.e. prob (0,1); mean; sdvals (0, Inf)
dnormmix0 <- function(x, probs, muvals, sdvals, log = FALSE) {
  np <- length(probs)
  pmat <- matrix(NA, nrow = length(x), ncol = np)
  for (i in 1:np) {
    pmat[,i] <- dnorm(x, mean = muvals[,i], sd = sdvals[,i])
  }
  lik <- pmat |> sweep(MARGIN = 2, FUN = "*", probs) |> rowSums()
  if (log) log(lik) else lik
}

## takes a single vector of parameters on the unconstrained scale
##  (i.e. softmax(prob), mean, log(sd))
dnormmix <- function(x, par, logmean, ..., log = FALSE) {
  g0 <- genmixpars(logmean, par, ...)
  do.call(dnormmix0, c(list(x), g0, list(log = log)))
}

nllfun <- function(par, vals, logmean, np, sd_ord) {
  -sum(dnormmix(x = vals, par, logmean, np = np, sd_ord = sd_ord, log = TRUE))
}


gen_parnames <- function(np, sd_ord) {
  ## taken from contr.poly:
  sdlabs <- c(".1", ".L", ".Q", ".C", paste0(".", 4:10))[1:(sd_ord+1)]
  c(paste0("logitprob_", 1:(np-1)),
    c(outer(1:np, c("int", "slope"), function(x,y) sprintf("mu_%s_%d", y, x))),
    c(outer(1:np, sdlabs, function(x,y) sprintf("logsd_%s_%d", y, x))))
}



#############################################################################
#' Title
#'
#' @param metadata_list : list of metadata
#' @param countdata_list : list of otu tables
#' @param num_cores : number of cores
#' @param ref_name reference for fold change calculation
#'
#' @return
#' @export
#'
#' @examples
deseq_fun_est <-function(metadata_list,  countdata_list,
                         num_cores=10, ref_name= "control"){
  
  registerDoParallel(cores = num_cores)
  l = length(countdata_list)
  
  dds  =  foreach(i= 1:l, .packages = "DESeq2", .export = "deseqfun") %dopar%{
    countdata =  countdata_list[[i]]
    metadata  =  metadata_list[[i]]
    stopifnot(!is.na(sum(countdata)))
    stopifnot(colnames(countdata) == metadata$Samples)
    deseqfun(countdata,metadata, ref_name=ref_name,
             minReplicatesForReplace = Inf, 
             cooksCutoff = FALSE,
             independentFiltering = FALSE)
  }

  stopImplicitCluster()
  unregister_dopar() 
  names(dds)  =  names(countdata_list)
  dds
}

#############################################################
contour_plot_fun <- function(combined_data, 
                        power_estimate, 
                        cont_breaks){
  
  combined_data$pvalue_reject <- factor(combined_data$pval_reject)
 
  gg_2dimc <- (ggplot(combined_data)
               + aes(lmean_abund, abs_lfc)
               + rasterise(geom_point(aes(color = pvalue_reject), alpha = 0.5))
               + xlab(TeX("$\\log_2$(mean counts)")) 
               + ylab(TeX("|$\\log_2$(fold change)|")) 
               + scale_colour_manual(values = c("black", "red"))
               + geom_contour(data = power_estimate,
                              aes(z=power),lwd=1,
                              breaks = cont_breaks)
               + geom_label_contour(data = power_estimate, 
                                    aes(z= power,label = sprintf("%.3f", after_stat(level))),
                                    breaks = cont_breaks
               )
               
  )
  gg_2dimc
  
}

#######################################################################
inverse_fun <- function(target,lmb,abs_lfc, model,xmin, xmax) {
  
  if (length(target) != length(lmb) || length(target) != length(abs_lfc)) {
    stop("Lengths of 'target', 'lmb', and 'abs_lfc' must be the same.")
  }
  
  # Initialize an empty list to store the roots
  roots <- list()
  
  # Loop through each element and solve for the root
  for (i in 1:length(target)) {
    roots[[i]] <- uniroot(function(s) {
      predict(model,
              type = "response",
              newdata = data.frame(sample_size = s,
                                   lmean_abund = lmb[i],
                                   abs_lfc = abs_lfc[i])
      ) - target[i]
    },
    interval = c(xmin, xmax),
    extendInt = "yes")$root
  }

  unlist(roots)
}


uniroot_lmb =  function(target_power,sample_size,abs_lfc,model,xmin,xmax){
 
   root <- uniroot(function(lmb) {
    predict(model,
            type = "response",
            newdata = data.frame(sample_size = sample_size,
                                 lmean_abund = lmb,
                                 abs_lfc = abs_lfc)
    ) - target_power
  },
  interval = c(xmin, xmax),
  extendInt = "yes")$root
  root
}



inverse_lfc <- function(target,lmb,samp_size, model,xmin, xmax) {
  
  if (length(target) != length(lmb) || length(target) != length(samp_size)) {
    stop("Lengths of 'target', 'lmb', and 'abs_lfc' must be the same.")
  }
  
  # Initialize an empty list to store the roots
  roots <- list()
  
  # Loop through each element and solve for the root
  for (i in 1:length(target)) {
    roots[[i]] <- uniroot(function(s) {
      predict(model,
              type = "response",
              newdata = data.frame(sample_size =  samp_size[i],
                                   lmean_abund =  lmb[i],
                                   abs_lfc     =  s)
      ) - target[i]
    },
    interval = c(xmin, xmax),
    extendInt = "yes")$root
  }
  
  unlist(roots)
}



inverse_lmb <- function(target,abs_lfc,samp_size, model,xmin, xmax) {
  
  if (length(target) != length(abs_lfc) || length(target) != length(samp_size)) {
    stop("Lengths of 'target', 'lmb', and 'abs_lfc' must be the same.")
  }
  
  # Initialize an empty list to store the roots
  roots <- list()
  
  # Loop through each element and solve for the root
  for (i in 1:length(target)) {
    roots[[i]] <- uniroot(function(lmb) {
      predict(model,
              type = "response",
              newdata = data.frame(sample_size =  samp_size[i],
                                   lmean_abund =  lmb,
                                   abs_lfc     =  abs_lfc[i])
      ) - target[i]
    },
    interval = c(xmin, xmax),
    extendInt = "yes")$root
  }
  
  unlist(roots)
}





uniroot_ss =  function(target_power,lmean_abund, abs_lfc,model,xmin,xmax){
  
  root <- uniroot(function(ss) {
    predict(model,
            type = "response",
            newdata = data.frame(sample_size = ss,
                                 lmean_abund = lmean_abund,
                                 abs_lfc = abs_lfc)
    ) - target_power
  },
  interval = c(xmin, xmax),
  extendInt = "yes")$root
  root
}

# Function to compute the missing input
power.nb <- function(power, sample_size, logfoldchange, logmean_abundance, gam_mod) {
  
  if (is.numeric(power) && is.numeric(sample_size) && is.numeric(logfoldchange) && is.numeric(logmean_abundance)) {
    if (!missing(power) && !missing(sample_size) && !missing(logfoldchange)) {
      if (missing(logmean_abundance)) {
        # Compute the missing input
        
      }
    } else if (!missing(power) && !missing(sample_size) && !missing(logmean_abundance)) {
      # Compute the missing input
      logfoldchange <- logmean_abundance - power - sample_size
    } else if (!missing(power) && !missing(logfoldchange) && !missing(logmean_abundance)) {
      # Compute the missing input
      sample_size <- logmean_abundance - power - logfoldchange
    } else if (!missing(sample_size) && !missing(logfoldchange) && !missing(logmean_abundance)) {
      # Compute the missing input
      power <- logmean_abundance - sample_size - logfoldchange
    } else {
      stop("Please specify at least 3 inputs.")
    }
    
    # Return the computed inputs
    return(list(power = power, sample_size = sample_size, logfoldchange = logfoldchange, logmean_abundance = logmean_abundance))
  } else {
    stop("Inputs must be numeric values.")
  }
}




