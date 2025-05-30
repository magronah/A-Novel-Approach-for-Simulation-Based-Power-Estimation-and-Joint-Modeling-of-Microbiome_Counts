seed = 101; ntaxa = 100; nsubj = 50; beta = 3; betadisp = 0 
nsim =  10 ; n_gt = 2

theta_true =  c(
  get_theta_logSD(n_gt, seed = seed),
  get_theta_corr1(n_gt, n_gt,  seed=seed),
  get_theta_logSD(ntaxa, seed = seed),
  get_theta_corr1(ntaxa, nsubj,  seed=seed)
)
