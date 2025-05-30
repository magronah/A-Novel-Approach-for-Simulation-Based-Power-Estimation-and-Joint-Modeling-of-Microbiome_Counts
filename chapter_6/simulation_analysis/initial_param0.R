seed = 101; nsim  = 300
#######################################################
beta  =  5; n_gt  =  4
#######################################################
betazi   =  qlogis(((0.12) +(0.92))/2)
thetazi  =  log(diff(qlogis(c(0.12, 0.92)))/4)
pho      =  0.7
#######################################################
param_list  =  list(param1 = c(ntaxa = 100, nsubj = 30, ntime = 3),
                    param2 = c(ntaxa = 200, nsubj = 50, ntime = 4))
#########################################################
param  =   param_list[[param_index]]
ntaxa  =   param[["ntaxa"]]
nsubj  =   param[["nsubj"]]
ntime  =   param[["ntime"]]
#######################################################
theta <- c(
  ## AR1 log-sd and transformed cor
  # log(0.1), pho/sqrt(1-pho^2),
  ##########################################
  get_theta_logSD(n_gt, seed =  seed),   #logsd for group*time term
  get_theta_corr(n_gt,n_gt, seed = seed), #correlation for group*time term
  ##########################################
  log(0.1),                                 #logsd for nugget term
  ##########################################
  log(0.2),  #increase this to 0.2 and see
  ##########################################
  get_theta_logSD(ntaxa, seed = seed),   #might be this
  get_theta_corr(ntaxa, nsubj,seed = seed)
)

