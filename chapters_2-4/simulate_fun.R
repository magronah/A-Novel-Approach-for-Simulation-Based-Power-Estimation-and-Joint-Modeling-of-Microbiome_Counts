#######################################################
#' Title
#'
#' @param logmean_param log mean parameter
#' @param notu number of otus
#'
#' @return
#' @export
#'
#' @examples
logmean_sim_fun = function(logmean_param,notu){
  
  l=length(logmean_param$mu)
  if(l == 1){
    logmean = rnorm(n = notu,
                    mean = logmean_param$mu, 
                    sd   = logmean_param$sigma)
  }else{
    logmean = mixtools::rnormmix(n=notu, 
                                 lambda = logmean_param$lambda, 
                                 mu = logmean_param$mu,
                                 sigma = logmean_param$sigma)
  }
  logmean
  
}

##############################################################
#' Title
#'
#' @param logmean_sim: simulated log mean abundance 
#' @param logfoldchange_param: list containing the
#'             par = optimal parameters for log fold change fit
#'             np = optimal number of components for log fold change
#'             sd_order: order polynomial for standard deviation parameter
#'                     for log fold change
#' @param max_lfc: maximum log fold change value to be simulated
#' @param max_iter: maximum iteration while simulating logfoldchange < max_lfc
#'
#' @return: lfc: vector of  log fold change values
#' @export
#'
#' @examples
logfoldchange_sim_fun <- function(logmean_sim, logfoldchange_param,
                              max_lfc = 15, max_iter = 10000) {
  
    par   =   logfoldchange_param$par
    np    =   logfoldchange_param$np
  sd_ord  =   logfoldchange_param$sd_ord
  
  lfc <- myrnormmix(par, logmean_sim, np = np, sd_ord = sd_ord)
  r <- range(lfc);   iteration_count <- 0
 
  while (max(abs(r)) > max_lfc && iteration_count < max_iter) {
    lfc <- myrnormmix(par, logmean_sim, np = np, sd_ord = sd_ord)
    r <- range(lfc)
    iteration_count <- iteration_count + 1
  }
  
  if (iteration_count == max_iter) {
    warning("Maximum number of iterations reached without convergence.")
  }
  
  return(lfc)
}

##############################################################
#' Title
#'
#' @param logmean_param :  log mean parameters 
#' @param logfoldchange_param : log fold change parameters 
#' @param dispersion_param : dispersion parameters 
#' @param nsamp_per_group : number of samples per group
#' @param ncont: number of control samples 
#' @param ntreat: number of treatment samples 
#' @param notu : number of otus
#' @param nsim :number of simulations 
#' @param disp_scale : scale parameter for dispersion
#' @param max_lfc : maximum log fold change
#' @param maxlfc_iter : maximum number of iterations
#' @param seed : seed for simulation
#'
#' @return
#' @export
#'
#' @examples
countdata_sim_fun <- function(logmean_param, logfoldchange_param, dispersion_param, 
                              nsamp_per_group = NULL, ncont = NULL,ntreat = NULL,
                              notu, nsim = 1,   disp_scale=0.3, max_lfc = 15, 
                              maxlfc_iter = 1000, seed = NULL){
  
  # Check if both nsamp_per_group and (ncont or ntreat) are provided
  if (!is.null(nsamp_per_group) && (!is.null(ncont) || !is.null(ntreat))) {
    stop("Please specify either 'nsamp_per_group' or 'ncont' and 'ntreat', but not both.")
  }
  
  if (is.null(nsamp_per_group)) {
    if (is.null(ncont) || is.null(ntreat)) {
      stop("Please specify both 'ncont' and 'ntreat' when 'nsamp_per_group' is not provided.")
    }
  }
  
  par     =   logfoldchange_param$par
  np      =   logfoldchange_param$np
  sd_ord  =   logfoldchange_param$sd_ord
  
  pb <- txtProgressBar(0, nsim, style = 3)
  set.seed(seed)

  res = foreach(j = 1:nsim, 
                .export = c("logfoldchange_sim_fun","logmean_sim_fun",
                            "myrnormmix", "rnormmix0","genmixpars",
                            "dispersion_fun"),
                .packages = c("stats","purrr", "mixtools"),.verbose = F) %do% {
                
              setTxtProgressBar(pb, j)
              
              ####Simulate log mean and log fold change    
              logmean        =   logmean_sim_fun(logmean_param,notu)
              logfoldchange  =   logfoldchange_sim_fun(logmean,logfoldchange_param, 
                                                   max_lfc = max_lfc,
                                                   max_iter = maxlfc_iter)
              
              ####Calculate control and treatment mean abundances    
              mean_abund =  2^logmean; foldchange = 2^logfoldchange
              control    =  (2*mean_abund)/(1+ foldchange) 
              
              ####predict dispersion    
              dispersion  =   mean_abund |>
                map(function(x) dispersion_fun(x, dispersion_param$asymptDisp,
                                                 dispersion_param$extraPois)) |>  unlist()
          
              dd = data.frame(control = control,
                              treatment = control*foldchange,
                              dispersion = dispersion)
                
              ####Simulate count data    
              if(!is.null(nsamp_per_group)){
                countdata = data.frame(apply(dd, 1, function(x) {
                                       rnbinom(n=2*nsamp_per_group,
                                               mu = rep(c(x["control"],x["treatment"]), 
                                                        each= nsamp_per_group),
                                               size = 1/(disp_scale*x["dispersion"]))
                }))
                
                metadata =  data.frame(Samples=paste0("sample_",1:(2*nsamp_per_group)),
                                       Groups=factor(rep(c("control", "treatment"),
                                                         each=(nsamp_per_group))))
                
                rownames(countdata) =   paste0("sample_",1:(2*nsamp_per_group))
                colnames(countdata) =   paste0("otu_",1:notu)
                
                control_count       =   countdata[1:nsamp_per_group, ]
                treat_count         =   countdata[(nsamp_per_group+1):(2*nsamp_per_group), ]
                
              }else{
                n  =  ncont + ntreat
                countdata = data.frame(apply(dd, 1, function(x) {
                                    rnbinom(n = n,
                                            mu = rep(c(x["control"],x["treatment"]), 
                                                     c(ncont,ntreat)),
                                            size = 1/(disp_scale*x["dispersion"]))
                }))
                
                metadata =  data.frame(Samples=paste0("sample_",1:n),
                                       Groups=factor(rep(c("control", "treatment"),
                                                         c(ncont,ntreat))))
                
                rownames(countdata) = paste0("sample_",1:n)
                colnames(countdata) = paste0("otu_",1:notu)
                
                control_count = countdata[1:ncont,]
                treat_count   = countdata[(ncont+1):n,]
              }
              names(logmean)  =   paste0("otu_",1:notu)
              names(logfoldchange)  = paste0("otu_",1:notu)
              
          
              list(countdata     =   t(countdata), 
                   control_count =   t(control_count), 
                   treat_count   =   t(treat_count),
                   metadata      =   metadata, 
                   logmean       =   logmean, 
                   logfoldchange =   logfoldchange)
              }

 
  names(res)      =   paste0("sim_",1:nsim)
  
  ## extract individual results into a list
  list(  countdata_list  =   read_data(res,"countdata"),
         metadata_list   =   read_data(res,"metadata"),
         logmean_list    =   read_data(res,"logmean"),
         logfoldchange_list     =  read_data(res,"logfoldchange"),
         treat_countdata_list   =  read_data(res,"treat_count"),
         control_countdata_list =  read_data(res,"control_count")
         )

}

#######################################################
dispersion_fun <- function(mean_abund,asymptDisp,extraPois){
        asymptDisp + extraPois/mean_abund
}

