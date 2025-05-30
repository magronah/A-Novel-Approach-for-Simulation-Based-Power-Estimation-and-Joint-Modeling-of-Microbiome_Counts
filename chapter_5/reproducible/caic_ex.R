setwd("/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm/")
source("reproducible/leverage_funs.R")

## don't redo this unless necessary (slow)
if (FALSE) {
#/home/agronahm/projects/def-bolker/agronahm/reduced_rank_mm
    install_glmmTMB(pkgdir = "/home/agronahm/projects/def-bolker/agronahm/glmmTMB/glmmTMB",
                    libdir = "/home/agronahm/projects/def-bolker/agronahm/glmmTMB_lev")

   # install_glmmTMB(pkgdir = "~/Documents/glmmTMB/glmmTMB",
   #                 libdir = "~/Documents/glmmTMB_lev")
}


#library(glmmTMB, lib.loc = "/home/agronahm/projects/def-bolker/agronahm/glmmTMB_lev")
library(glmmTMB, lib.loc = "~/Documents/glmmTMB_lev")

## both of these packages *must* be loaded for leverage code to work
library(RTMB)
library(Matrix)

library(lme4)
library(cAIC4)

## testing leverage/cAIC computations on sleepstudy example 
data("sleepstudy", package = "lme4")
fm1 <- glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy)
fm0 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = FALSE)


X <- getME(fm1, "X")
Z <- getME(fm1, "Z")
## compute _conditional_ log likelihood of glmmTMB object?
pp <- fm1$obj$env$last.par.best
ppL <- split(pp, names(pp))

cAIC(fm0)
fm1 = mod1
## this successfully recovers the conditional log-likelihood as given by cAIC.
resp <- model.response(model.frame(fm1))
condlik1 <- sum(dnorm(resp,
           drop(X %*% ppL[["beta"]] + Z%*% ppL[["b"]]), sigma(fm1), log = TRUE))
## so does this
condlik2 <- sum(dnorm(resp, fitted(fm1), sigma(fm1), log = TRUE))

## this could *almost* be automated (with a case/switch statement mapping
## family to a d*() function)

## in general should extract any info we want from the TMB object *before*
## trying to compute leverage ... leverage computation screws up internal
## state of the TMB object

## doesn't quite match, but close
lfm1 <- leverage(fm1)
plot(lfm1, hatvalues(fm0))
abline(a=0, b=1, lty = 2)


cdf <- sum(lfm1)+1
c(clik = condlik2, cdf = cdf, caic = 2*(-condlik2 + cdf))
data  =  sleepstudy

y_hat <- function(model, data) {
  predict(model, newdata = data, type = "response", re.form = NULL)
}

H_matrix <- calc_leverage(fm1, data)
diag_H <- diag(H_matrix)
print(range(diag_H))   

sum_H <- sum(diag_H)
print(sum_H)

y_hat <- function(model, data) {
  predict(model, newdata = data, type = "response", re.form = NULL)
}

calc_leverage <- function(model, data,epsilon = 0.01) {
  n <- nrow(data)
  H <- matrix(0, n, n)  
  
  y_pred <- y_hat(model, data)  
  yname  <- as.character(formula(model)[[2]])
  
  for (i in 1:n) {
    data_perturb <- data
    data_perturb[[yname]][i] <-  data[[yname]][i] + epsilon
    model_up <- update(model, data = data_perturb)
    y_pred_up <- y_hat(model_up, data_perturb)
    
    H[, i] <- (y_pred_up - y_pred) / epsilon
  }
  
  return(H)
}

## not identical, but similar ...
sum(lfm1)
sum(hatvalues(fm0))

## Vaida and Blanchard add 1 to the number of params for the dispersion
cdf <- sum(lfm1)+1
cAIC(fm0)
c(clik = condlik2, cdf = cdf, caic = 2*(-condlik2 + cdf))

mm <- readRDS("reproducible/rr_mod.rds")
mm$fit[c("iterations", "evaluations")]
mf <- model.frame(mm)
length(unique(mf$taxon))  ## 969 taxa
length(unique(mf$site))   ## 8 'subjects'

## do this *before* trying to do leverage computation!
## this is a plausible value for the conditional log-likelihood ...
condlik3 <- sum(dnbinom(model.response(model.frame(mm)),
                        mu = fitted(mm), size = sigma(mm), log = TRUE))
## compare unconditional nll
mm$obj$fn()

conditional_ll <-  function(mod){
  
  condlik3 <- sum(dnbinom(model.response(model.frame(mod1)),
                          mu = fitted(mod1), size = sigma(mod1), 
                          log = TRUE))
  
  system.time(trace_hat <- sum(leverage(mod)))
  cdf3 <- trace_hat+1
  c(clik = condlik3, cdf = cdf3, caic = 2*(-condlik3 + cdf3))
}
  
if (FALSE) {
    ## uses up memory on my 60Gb system, kills R session
    system.time(trace_hat <- sum(leverage(mm)))
    cdf3 <- trace_hat+1
    c(clik = condlik3, cdf = cdf3, caic = 2*(-condlik3 + cdf3))
}

library(peakRAM)
testfun <- function(nsubj = 100, ntax = 100, d = 2, seed = 101) {
    set.seed(seed)
    dd <- expand.grid(subject = factor(seq(nsubj)),
                      taxon = factor(seq(ntax)))
    ## hard-code d=2 for now (issues with scoping/evaluation)
    dd$y <- simulate_new( ~ 1 + rr(taxon | subject, d = 2),
                 family = nbinom2,
                 newdata = dd,
                 newparams = list(beta = 1,
                                  betadisp = 1,
                                  theta = rep(0.1, ntax*d - choose(d,2))))[[1]]
    ## have to run this without parallelization, since we had to turn off
    ## OpenMP for leverage calculations (we could try to load the full
    ## version of glmmTMB with autopar for fitting the model, then
    ## detach and load the glmmTMB_lev
    p1 <- peakRAM(
        mod <- glmmTMB(y ~ 1 + rr(taxon | subject, d = 2),
                       family = nbinom2,
                       data = dd)
    )
    tmpf <- function(x, task = "model_fit") {
        names(x) <- c("task", "time_sec", "total_RAM_Mb", "peak_RAM_Mb")
        x$task <- task
        x <- data.frame(nsubj = nsubj, ntax = ntax, d = d, x)
        return(x)
    }
    p2 <- peakRAM(leverage(mod))
    rbind(tmpf(p1), tmpf(p2, task = "leverage"))
}

## test for a small example
testframe <- expand.grid(nsubj = seq(10, 40, by = 5),
                         ntax = seq(50, 200, by = 25))
res <- list()
for (i in seq(nrow(testframe))) {
    nsubj <- testframe$nsubj[i]
    ntax <- testframe$ntax[i]
    if (nsubj > 10 & ntax > 200) break
    cat(i, nsubj, ntax, "\n")
    res[[i]] <- peakRAM_testfun(nsubj = nsubj, ntax = ntax)
    saveRDS(res, "leverage_timings.rds")
}

library(ggplot2); theme_set(theme_bw())
res <- readRDS("leverage_timings.rds")
res <- do.call(rbind, res)

res_long <- res |>
    tidyr::pivot_longer(c(time_sec, matches("RAM")),
                 names_to = "measure")

gg0 <- ggplot(res_long, aes(nsubj, value, colour= factor(ntax))) +
    facet_grid(measure ~ task, scale = "free") +
    geom_line() + geom_point() +
    scale_x_log10() + scale_y_log10()

print(gg0)

print(gg0 + aes(x=ntax, colour = factor(nsubj)))

