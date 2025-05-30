###Function to fit GAM
gam_fit <- function(pvalue, effect_size, mean_count,
                    grid_len = 100, alpha_level = 0.05){
  
  pval_reject  =   (!is.na(pvalue) & pvalue < alpha_level)
  
  comb      =   tibble(lmean_count  =  log(mean_count),
                       abs_lfc      =  abs(effect_size),
                       pval_reject  =  as.numeric(pval_reject))
  
  
  #fit scams
  fit_2d       =    mgcv::gam(pval_reject ~ te(lmean_count, abs_lfc),
                         data = comb, family = binomial)
  # fit_2d       =    scam(pval_reject ~ s(lmean_count, abs_lfc, bs="tedmi"),
  #                        data = comb, family = binomial,
  #                        control = list(trace = TRUE, print.warn = TRUE,
  #                                       devtol.fit = 1e-2,
  #                                       steptol.fit = 1e-2,
  #                                       maxit = 1000))
  
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
pvalue_cal =  function(dd){
  apply(dd, 1, function(x){
    2*min(c(mean(x < 0), mean(x > 0)))})
}