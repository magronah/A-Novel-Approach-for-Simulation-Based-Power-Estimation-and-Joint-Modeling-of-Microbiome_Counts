library(glmmTMB)
mm <- readRDS("reproducible/rr_mod.rds")

## attempt leverage calculation

## reconstruct data from model frame
mf <- model.frame(mm)
newdata <- with(mf,
                data.frame(count, group, taxon, site,
                           normalizer = `offset(normalizer)`))
dim(newdata)
length(unique(newdata$taxon))

## get starting parameters 
p0 <- with(mm$obj$env, parList(last.par.best[-random]))
## exclude empty (unused) parameters
p0 <- p0[lengths(p0) > 0]


## exclude random effects
p0 <- p0[setdiff(names(p0), "b")]

## test (6 seconds)
system.time(
    lik1 <- with(mm$obj$env, mm$obj$fn(last.par.best[-random]))
)
lik2 <- with(mm$obj$env, mm$obj$fn(unlist(p0)))

## might not want threading if we're doing the leverage calculation in parallel ...
par_ctrl <- glmmTMBControl(
    parallel = list(n = 10, autopar = TRUE)
)

gprior <- data.frame(prior = "gamma(2, 2.5)",
                     class = "theta_sd",
                     coef = "")
## par_ctrl <- glmmTMBControl(optCtrl = list(trace=1))  ## default
par_ctrl <- glmmTMBControl()

options(glmmTMB_openmp_debug = FALSE)

#Run this as jobs. Run this as job
calc_leverage <- function(i, eps = 5) {
    i <- 800; eps <-  0#1e-4#
    newdata$count[i]
    pert_data <- newdata
    pert_data$count[i] <- pert_data$count[i] + eps
    ## 2 minutes to optimize?
    # st <- system.time(
    #     newfit <- update(mm, start = p0, data = pert_data, verbose = FALSE,
    #                      control = par_ctrl, priors = gprior)
    # )
    ## newfit <- update(mm, start = p0, data = pert_data, verbose = FALSE,
    ##                     control = par_ctrl, priors = gprior)

    ## can we do this modularly and skip the slow stuff after the optimization?
    newfit0 <- update(mm, start = p0, data = pert_data, verbose = FALSE,
                      control = par_ctrl, priors = gprior, doFit = FALSE)
    
    system.time(newfit1 <- fitTMB(newfit0, doOptim = FALSE)) ## 1 second
    system.time(newfit2 <- with(newfit1, nlminb(par, fn, gr))) ## 80 seconds
    newfit3 = finalizeTMB(newfit0,newfit1,newfit2)

    ## then we should be
    X <- getME(mm, "X")
    Z <- getME(mm, "Z")
    pp <- with(newfit1$env, parList(last.par.best[-random]))
    ## these bs from  parList(last.par.best[-random]) are the problem
    pp$b  <-  (newfit1$report()$b)
    ## perturbed prediction
    pred <- drop(X[i,] %*% pp[["beta"]] + Z[i,] %*% pp[["b"]])
    ## original prediction
    pred0  <- predict(newfit1$obj, type = "link")[i]
  
    pred0  <- predict(mm, type = "link")[i]
    # pred0 <- fitted(mm)[i]
    (pred - pred0) / eps
}

#epsilon needs to be large because it is going to be converted to the log scale
#################################################################

library(glmmTMB)

# Simulated data
set.seed(123)
data <- data.frame(
  y = rpois(100, lambda = 10),
  x = rnorm(100),
  group = factor(rep(1:10, each = 10))
)

mod <- glmmTMB(y ~ x + (1 | group), data = data, family = poisson())

y_hat <- function(model, data) {
  predict(model, newdata = data, type = "response", re.form = NA)
}

calc_leverage <- function(model, data,resp_name, epsilon = 1e-4) {
  n <- nrow(data)
  H <- matrix(0, n, n)  
  
  y_pred <- y_hat(model, data)  
  
  for (i in 1:n) {
    i = 1
    data_perturb <- data
    
    data_perturb$y[i] <- data$y[i] + epsilon
    model_up <- update(model, data = data_perturb)
    y_pred_up <- y_hat(model_up, data_perturb)
    
    H[, i] <- (y_pred_up - y_pred) / epsilon
  }
  
  return(H)
}

H_matrix <- calc_leverage(mod, data)

diag_H <- diag(H_matrix)
print(range(diag_H))   

sum_H <- sum(diag_H)
num_params <- length(fixef(mod)$cond) 
print(sum_H)
print(num_params)

data_outlier <- data
data_outlier$y[1] <- data$y[1] + 50  # Introduce an outlier
mod_outlier <- update(mod, data = data_outlier)
H_outlier <- calc_leverage(mod_outlier, data_outlier)

print(diag(H_matrix)[1])   # Before outlier
print(diag(H_outlier)[1])  # Should be higher after outlier

# Sanity Check 3: Compare leverage results for different epsilon values
H_eps1 <- calc_leverage(mod, data, epsilon = 1e-5)
H_eps2 <- calc_leverage(mod, data, epsilon = 1e-3)

print(max(abs(H_matrix - H_eps1)))  # Should be small
print(max(abs(H_matrix - H_eps2)))  # Should be small
#################################################################
