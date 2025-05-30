## goal: pick 'theta' parameters for a reduced-rank model in a sensible way
##' @param d rank (dimension)
##' @param n full dimension (latent variables per group)
##' @param logsdvec vector of log-SDs of each factor
gen_rr_parms <- function(n, d,logsdvec) {
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
  stopifnot(all.equal(sqrt(colSums(mat^2)), exp(logsdvec)))
  theta <- c(
    mat[row(mat)==col(mat)],  ## diagonal elements
    mat[row(mat)>col(mat)]    ## below-diagonal elements
  )
  return(theta)
}

gen_rr_parms_old <- function(n, d, sdvec = c(2,1)) {
  if (length(sdvec) != d) stop("length(sdvec) should equal d")
  F <- matrix(0, nrow = n, ncol = d)
  for (i in 1:d) {
    len <- n - i + 1
    x <- rnorm(len-1, mean = 0, sd = 2)
    ## scale so sum^2 = 1
    x_exp <- exp(c(1,x))
    x2 <- sqrt(x_exp/sum(x_exp))
    ## sqrt() automatically makes x2[1] >0,
    ##   satisfies uniqueness constraint
    F[ i:n, i] <- x2
  }
  ## scale columns by factor SD
  F <- sweep(F, 2, STATS = sdvec, FUN = "*")
  return(F)
}

rr_parms_vec <- function(F) {
  n <- ncol(F)
  result <- c()
  for (i in 1:n) {
    result <- c(result, F[-(1:i), i])
  }
  return(result)
}

rr_fac_len <- function(ntaxa, rank) {
  ntaxa*rank - choose(rank,2) - rank
}


get_theta_corr <- function(n, mat= NULL, seed = NULL) {
  if(!is.null(mat)){C  <- mat}
    else{set.seed(seed); C  <-  rlkjcorr(n=1,K=n,eta=1)}
  scale <- sqrt(fastmatrix::ldl(as.matrix(C))$d)
  cc2 <- chol(C) %*% diag(1/scale)
  cc2[upper.tri(cc2)]
}


get_theta_logSD <- function(n, seed = NULL, rank = NULL, prob = 0.5) {
  set.seed(seed)
  val    =   rgeom(n, prob=prob) + 0.1
  logSD  =   log(sqrt(val))
  
  if(!is.null(rank)){
    logSD     =   logSD[1:rank]
    return(logSD)
  }
  else{return(logSD)}
}


get_params <- function(ntaxa, beta  = 0, betadisp = 0,
                       param_type = c("full", "rr", "none"), 
                       rank = NULL, seed=NULL) {

  param_type = match.arg(param_type)
  set.seed(seed)
  if(param_type == "full"){

    n = ntaxa*(ntaxa - 1)/2
    theta <- c(
      ## log-sd of intercept, group,
      get_theta_logSD(n=2, seed = seed),
      ## covariance of intercept, group
      get_theta_corr(n = 2, seed = seed),
      ## log-sd of taxa
      get_theta_logSD(n=ntaxa, seed = seed),
      ## covariance
      get_theta_corr(n = ntaxa, seed = seed)
      )
    theta_nm <- c("intsd", "grpsd", "grp_int_cov",
                  paste0("fmod_sd", 1:ntaxa),
                  paste0("fmod_cov", seq(n)))

    names(theta) <- theta_nm
  } else if (param_type == "rr"){

    if(is.null(rank)){stop("specify the rank. Example rank=2")}

    n = ntaxa*rank - rank*(rank-1)/2
    theta = c(
      ## log-sd of intercept, group,
      get_theta_logSD(n=2, seed = seed),
      ## covariance of intercept, group
      get_theta_corr(n = 2, seed = seed),
      ## RR factors
       get_theta_logSD(n=ntaxa, seed = seed,rank=rank),
       rnorm(rr_fac_len(ntaxa, rank))
      )
    theta_nm <- c("intsd", "grpsd", "grp_int_cov",
                  paste0("RR_sd", 1:rank),
                  paste0("RR_load", seq(n-rank)))

    names(theta) <- theta_nm

  } else {
      theta <- c(
          ## log-sd of intercept, group,
          get_theta_logSD(n=2, seed = seed),
          ## covariance of intercept, group
          get_theta_corr(n = 2, seed = seed)
          )
  }
    return(list(beta = beta, theta = theta, betadisp=betadisp))
}


 


sim_glmmTMB <- function(form, dd, beta = 0,rank=NULL, seed= NULL, nsim = 1,
                        param_type=c("full", "rr", "none"), return_val = "sim") {

  param_type = match.arg(param_type)
  ntaxa  =  length(unique(dd$taxa))
  parms   =  get_params(ntaxa, beta  =  beta, param_type = param_type, 
                       rank  =  rank, seed   =  seed)

  ## how long should theta be?
  if(param_type=="full"){
    ntheta <- 3    +     ## (2 log-sd and 1 covariance for 1+group)
           ntaxa*(ntaxa+1)/2  ## (ntaxa log-sd and ntaxa*(ntaxa-1)/2 covariance 
                               ##     for 0 + ntaxa term  t )
  }else if (param_type=="rr"){
    ntheta <- 3 +  ## (2 log-sd and 1 covariance for 1+group)
           ntaxa*rank - rank*(rank-1)/2  # rr(taxa + 0 | subject,r)
  } else { ## none
      ntheta <- 3
  }
  stopifnot(length(parms$theta) == ntheta)
  
  if (return_val == "pars") {
    result <- simulate_new(form, nsim = nsim, seed = seed,
                           newdata   = dd,
                           newparams = parms,
                           return_val = "pars",
                           family = nbinom2)
    return(result)
  } else if (return_val == "sim") {
    result <- simulate_new(form, nsim = nsim, seed = seed,
                           newdata = dd,
                           newparams = parms,
                           return_val = "sim",
                           family = nbinom2)
    return(result)
  } else if (return_val == "object") {
    result <- simulate_new(form, nsim = nsim, seed = seed,
                           newdata    = dd,
                           newparams  = parms,
                           return_val = "object",
                           family = nbinom2)
    return(result)
  } else {stop("Invalid return_val. Use 'pars',object or 'sim'.")}
}

##' Simulation of theta: procedure
##' I simulated a random matrix and then constructed a positive definite matrix
##' (because variance covariance matrix needs to be positive definite)
##' simply by taking the cross product and then compute the eigen 
##' values and eigen vectors of the positive definite matrix
##' I used the log(sqrt( of the diagonals)) as my logSD parameters and 
##' the first d eigen vectors concatenated rowwise as my loading parameter
##' Is the reduced rank results ordered?
##' Is the full rank result ordered?
# n=100; p <- 0.9
# data <- rgeom(n, p)
# plot(data[order(data,decreasing = T)])
# plot(data)
# log(sqrt(data))
# order()