#############################################################################
#rm(list=ls(all=TRUE))
set.seed(79861)
library(lattice)
library(nlme)


### Set working directory  #####
## Under this directory, should have "BayesRegMicrobiome_1.0.tar.gz"
#### Path  to the R package file to install
#### Once the package is installed, do not need to install again.

#install.packages("BayesRegMicrobiome_1.0.tar.gz", repos = NULL, type="source")

### load a library
library(BayesRegMicrobiome)



### load the simulated data for Section 3. Simulation Study
load("SIM_dat.RData")
View(SIM_dat)
View(SIM_dat$Y)
## SIM_dat
#> names(SIM_dat)
#[1] "rep_K"        "Tot_N"        "n"            "P"            "J"
#[6] "Z"            "Z_1"          "Y"            "max_cat"      "Miss_ind"
#[11] "miss_cov_ind"

# n: the number of time points, t_i \in [0, T], i=1, ..., n
# rep_K (K_i): n-dim vector, the number of replicates at time point t_i
# Tot_N = \sum_i K_i: the total # of samples
# P: the number of covariates
# J: the number of OTUs
# Z: Covariates recorded at each t_i  -- (n*(P+2)) matrix where columns 1 and 2 are i and t, respectively, and the remaining columns are Z
# Z_1: Covaraites for each sample (t_i, k) -- (Tot_N*(P+2)) matrix where columns 1 and 2 are i and t, respectively, and the remaining columns are Z
# Y: OTU counts (Tot_N*J) matrix
# max_cat: # of max categories for discrete covariates, -1 for continuous.  The method can impute missing *discrete* covaraites.
# Miss_ind: indicator of where a covariate is missing at any of n time points (any missing: 1, none missing: 0)   -- (n*(P+2)) matrix where columns 1 and 2 are i and t, respectively, and the remaining columns are Z
View(SIM_dat$Z_1)
SIM_dat$Z_1[,2]  = 1
### specify fixed hyperparmeters
hyper <- fn.hyper(SIM_dat$Y, SIM_dat$Tot_N, SIM_dat$Z_1[,2], max(SIM_dat$Z_1[,2]), SIM_dat$P, SIM_dat$J)

hyper <- fn.hyper(SIM_dat$Y, SIM_dat$Tot_N, SIM_dat$Z_1[,2], max(SIM_dat$Z_1[,2]), SIM_dat$P, SIM_dat$J)

View(SIM_dat)
#> names(hyper)
#[1] "a_s"    "b_s"    "L"      "u2"     "c_r"    "a_w"    "b_w"    "v2_eta"
#[9] "d"      "L_th0"  "u2_th0" "c_th0"  "aw_th0" "bw_th0" "v2_th0" "d_th0"
#[17] "a_sig"  "b_sig"  "a_lam"  "b_lam"  "M"      "u"      "gam2"   "K"
#[25] "a_tau"  "b_tau"

# s_j \iid Ga(a_s, b_s)
# (L, u2, c_r, a_w, b_w, v2_eta, d) are the fixed hyperparameters for the prior of \tilde{r}_{tk} = log(r_{tk})
# (L_th0, u2_th0, c_th0, aw_th0, bw_th0, v2_th0, d_th0) are the fixed hyperparameters for the prior of \tilde{alpha}_{0j} = log(alpha_{0j})
# sig2_j \iid IG(a_sig, b_sig)
# lam2_j \iid Ga(a_lam, b_lam)
# M: # of basis points for the process convolution GP model
# u: a set of basis points
# gam2: the fixed parameter for the kernel
# K: a matrix of Z(t_i - u_m)
# tau2_j \iid IG(a_tau, b_tau)


###  MCMC parameters
NN_Burn <- 10  ## No of iterations for Burn-in period
NN_sam <- 100  ## Posterior sample size
Bayes_Reg_Microbime1 =  function (hpara, Dat, n_burn, n_sam) 
{
  # Dat$Y  otu table
  View(Dat$Z_1)
  View(Dat$P)
  View(Dat$Tot_N)
  # Dat$Tot_N number of subjects 
  # Dat$J  number of taxa 
  # Dat$Z_1 in my case Z_1, index, subject and group
  
  ini_sam <- fn.initialize(hpara, Dat$Y, Dat$Z_1, Dat$Tot_N, 
                           Dat$P, Dat$J)
  print(date())
  set.seed(79861)
  ini_sam_1 <- fn.run.MCMC(hpara, Dat, ini_sam, n_burn)
  set.seed(79861)
  print(date())
  Save_MCMC_sam <- fn.run.MCMC.1(hpara, Dat, ini_sam_1, n_sam)
  return(Save_MCMC_sam)
}

names(SIM_dat)
MCMC_sam <- Bayes_Reg_Microbime(hyper, SIM_dat, NN_Burn, NN_sam)
names(MCMC_sam)

class(SIM_dat$Z_1)
## SAVED MCMC SAMPLES
#> names(MCMC_sam)
#[1] "beta"    "phi"     "s"       "r"       "th0"     "th"      "sig2"
#[8] "tau2"    "lam2"    "Z_1_cov"

# beta: P*J*N_sam matrix of posterior samples of size N_sam for beta_{jp}, j=1, ..., J, and p=1, ..., P
# phi: P*J*N_sam matrix of posterior samples of size N_sam for phi_{jp}, j=1, ..., J, and p=1, ..., P
# s: J*N_sam matrix of posterior samples of size N_sam for s_{j}, j=1, ..., J
# r: Tot_N*N_sam matrix of posterior samples of size N_sam for r_{t_i,k}, i=1, ..., n and k=1, ..., K_i (that is, Tot_N many (i, k))
# th0: J*N_sam matrix of posterior samples of size N_sam for alpha_{0j}, j=1, ..., J
# sig2: J*N_sam matrix of posterior samples of size N_sam for sig2_{j}, j=1, ..., J
# tau2: J*N_sam matrix of posterior samples of size N_sam for tau2_{0j}, j=1, ..., J
# lam2: J*N_sam matrix of posterior samples of size N_sam for lam2_{0j}, j=1, ..., J
# Z_1_cov: Tot_N*(# of covaraites with missing values)*N_sam of imputed covariates of size N_sam, X_{t,p} at all K_i







