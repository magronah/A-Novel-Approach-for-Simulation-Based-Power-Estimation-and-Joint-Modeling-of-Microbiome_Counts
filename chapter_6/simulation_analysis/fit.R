library(RhpcBLASctl)
library(glmmTMB)
source("reproducible/longitudinal/func2.R")
path = paste0(getwd(),"/reproducible/longitudinal/data/")
###########################################################
data       =   readRDS(paste0(path,"dd_long.rds"))

form2      =   count ~1  +  offset(normalizer) +
                      (group*time|taxon) +
                      (1|nugget) + 
                      #(1|subject:taxon)
                      rr(taxon + 0 | subject:time,2) 
##########################################################
forcats::fct_inorder()forcats::fct_inorder()
as.numeric(factor(x, levels = unique(x)))


mod    =   mms(y = ddd, fixed = ~group*time  + offset(normalizer),
               random = ~ 1|subject,
               correlation = corAR1(form = ~ as.numeric(
                 factor(time, levels = unique(time)))| subject),
               zi_fixed = ~1,
               data = meta,
               method = "zinb",
               niter = 100)


# (0 + time|subject:taxa) ~ 1 +time |

# mod    =   mms(y = ddd, fixed = ~group*time  + offset(normalizer),
#                random = list(subject = pdDiag(~time))),
#                data = meta, 
#                method = "nb",
#                niter = 100)

##Fit the model
mod_list  =  list()

par_ctrl <- glmmTMBControl(
  parallel = list(n = 10, autopar = TRUE)
)

for(i in 1:length(data)){
  
  dat = data[[i]]
  
  gprior  <- data.frame(prior = "gamma(2, 2.5)",
                        class = "theta_sd",
                        coef = "")
  
  options(glmmTMB_openmp_debug = TRUE)
  blas_set_num_threads(1)
  
  system.time(
    fit  <-  glmmTMB(form2, data  = dat,
                     family  = poisson, 
                     ziformula  = ~1 + (1|taxon),
                     prior   = gprior,
                     REML    = TRUE,
                     control = par_ctrl
    )
  )
  
  # par_ctrl <- glmmTMBControl(
  #   parallel = list(n = 1, autopar = TRUE)
  # )
  
  saveRDS(fit, file=paste0(path,  "mod",i,".rds"))
  print(i)
  #mod_list[[i]]  =  fit
}

#11  length(theta)
# est2  =  fit$fit$parfull[names(fit$fit$parfull) == "theta"]
# 
# plot(theta[1:11],getME(fit, "theta")[1:11])
# abline(0,1)
# 
# est22       =   est2[seq(n_gt,length(est2),n_gt)]
# length(getME(fit, "theta"))
# 
# saveRDS(mod_list, file=paste0(path,  "mod_list.rds"))
# 
