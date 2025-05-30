source("reproducible/leverage_funs.R")

## adjust as required
# glmmTMB_lib <- "./glmmTMB_lev"
glmmTMB_lib <- "~/Documents/glmmTMB_lev"


## install_glmmTMB("~/R/pkgs/glmmTMB/glmmTMB", glmmTMB_lib)


library(lme4)
library(cAIC4)
library(glmmTMB, lib.loc = glmmTMB_lib)
## both of these packages *must* be loaded for leverage code to work
library(RTMB)
library(Matrix)

## 1. compare cAIC results for sleepstudy model:
##   lme4: cAIC4, hatvalues()
##   glmmTMB: 

data("sleepstudy", package = "lme4")
dd <- transform(sleepstudy, nbReact  = ceiling(Reaction))
fm0 <- lmer(Reaction ~ Days + (Days | Subject), dd, REML = FALSE)
fm1 <- glmmTMB(Reaction ~ Days + (Days | Subject), dd)

X <- getME(fm1, "X")
Z <- getME(fm1, "Z")
## compute _conditional_ log likelihood of glmmTMB object?
pp <- fm1$obj$env$last.par.best
ppL <- split(pp, names(pp))

cAIC(fm0)

## this successfully recovers the conditional log-likelihood as given by cAIC.
resp <- model.response(model.frame(fm1))
condlik <- sum(dnorm(resp, fitted(fm1), sigma(fm1), log = TRUE))

## look at leverage vs. epsilon for a few haphazardly selected indices
eps_vec <- 10^seq(-4, 0, length = 51)
lev_vals <- sapply(eps_vec, \(e) leverage_brute(fm1, data = sleepstudy, inds = 1, eps = e))
lev_vals50 <- sapply(eps_vec, \(e) leverage_brute(fm1,  data = sleepstudy, inds = 50, eps = e))
lev_vals100 <- sapply(eps_vec, \(e) leverage_brute(fm1,  data = sleepstudy, inds = 100, eps = e))

cvec <- c(1,2,4)
par(las=1); matplot(eps_vec,
                    cbind(lev_vals, lev_vals50, lev_vals100),
                    pch = 1,
                    col = cvec,
                    type = "p",
                    log = "x", ylim = c(0, 1), xlab = "epsilon", ylab = "leverage")
abline(h=hatvalues(fm0)[c(1,50,100)], col = cvec, lty = 2)

sleepstudy[c(1,50,100),]
## lme4-computed hatvalues for day 9 are identical

## try it for all observations, for different epsilon -- how well do they line up?
## these (and parallel calcs below) take a little while, could benefit from foreach/furrr for parallel calc ...
## might be a memory limitation on fully parallel calcs?
lev_vals2_eps.001 <- leverage_brute(fm1, eps = 1e-2,  data = sleepstudy, progress = TRUE)
lev_vals2_eps.010 <- leverage_brute(fm1, eps = 1e-1,  data = sleepstudy, progress = TRUE)
lev_vals2_eps.100 <- leverage_brute(fm1, eps = 1e0,  data = sleepstudy, progress = TRUE)
## this will probably screw up the object, don't try anything after this ...
lev_vals_TMB <- leverage(fm1)

pfun <- function(x, y, ...) { points(x, y, ...); abline(a=0, b=1, lty = 2) }
all_levs <- cbind(lmer = hatvalues(fm0),
                  eps001 = lev_vals2_eps.001,
                  eps010 = lev_vals2_eps.010,
                  eps100 = lev_vals2_eps.100,
                  TMB = lev_vals_TMB)
pairs(all_levs, gap = 0, panel = pfun)

## conclusions: in this case eps=0.01 matches TMB 'fancy' computation very well (better than lme4/TMB match)
## eps = 1 is terrible (this is weird -- are these real values or is something somehow breaking in some other way?

## now try the same thing with a small nbinom2 reduced rank model.

set.seed(1)
## rr()
form <- y ~ 1 + x + rr(v|subject, d=2); thetavec <- rnorm(4*2-1, mean = 0, sd = 0.1)
## plain
## form <- y ~ 1 + x + (1|subject); thetavec <- 0
dd <- expand.grid(subject = factor(1:30), x = 1:4, v = factor(1:4))
dd$y <- simulate_new(form[-2],
                     newdata = dd,
                     family = nbinom2,
                     newparams = list(beta = c(1, 0.1),
                                      betadisp = 0,
                                      theta = thetavec))[[1]]
summary(dd$y)                     
fm2 <- glmmTMB(form, family = nbinom2, data = dd)
leverage_brute(fm2, data = dd, inds = 1, eps = 0.1)
## re-fit full model: different, but not very different
leverage_brute(fm2, data = dd, inds = 1, eps = 0.1, fit_method = "update")
leverage_brute(fm2, data = dd, inds = 1, eps = 0.1,pred_method = "predict", 
               fit_method = "update")

## results on response scale are also somewhat sensible
leverage_brute(fm2, data = dd, inds = 1, eps = 0.1, scale = "response")
leverage_brute(fm2, data = dd, inds = 1, eps = 0.1, fit_method = "update", scale = "response")

eps_vec2 <- 10^seq(-4, 0, length = 31)
lev_vals_nb <- sapply(eps_vec2, \(e) leverage_brute(fm2, inds = 1, data = dd, eps = e))
## weird, but not varying by a huge amount at least ... 
plot(eps_vec2, lev_vals_nb, log = "x")


evals <- 10^(-3:0)
lev_vals <- lapply(evals,
    \(e) leverage_brute(fm2, eps = e, data = dd, progress = TRUE, scale = "response")
)

lev_vals2_TMB <- leverage(fm2)

e_labs <- trimws(format(evals, scientific = FALSE))
lev_mat <- do.call(cbind, lev_vals)
colnames(lev_mat) <- paste0("eps", e_labs)
all_levs2 <- cbind(lev_mat,
                   TMB = lev_vals2_TMB)
pairs(all_levs2, gap = 0, panel = pfun)
## matches fancy TMB output well, if we use response scale

## we need to 'bootstrap' (in an informal sense): suppose we have

set.seed(101)
samp <- sample(nrow(dd), 20)
## suppose we only have these leverage values
sub_levs <- lev_vals2_TMB[samp]
m_sub <- mean(sub_levs)
sd_sub <- sd(sub_levs)
tt = 2
(nrow(dd)^2*(sd_sub^2))/tt^2
m_sub*nrow(dd)

sum(lev_vals2_TMB)
hist(lev_vals2_TMB, breaks = 40)
####
## tests on (subsets of) real data
m_big <- readRDS("reproducible/rr_mod.rds")
m_big <- up2date(m_big)
m_big$call
## need gprior and par_ctrl in order to refit
par_ctrl <- glmmTMBControl()
gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")


df <- model.frame(m_big) |>
    dplyr::rename(normalizer = "offset(normalizer)")

eps_vec <- c(0.01, 0.03, 0.1, 0.3, 1, 3, 10)
m_big_Lvals <- sapply(eps_vec,
                      \(e) leverage_brute(m_big, data = df, inds = 1, eps = e,
                                          opt_args = list(control = list(trace = 50))))
plot(eps_vec, m_big_Lvals, type = "b", ylim = c(0,1))
eps_vec2 <- seq(0.1, 10, length = 31)

library(parallel)
cl <- makeCluster(12)
clusterExport(cl, c("m_big", "df", "gprior", "par_ctrl", "glmmTMB_lib"))
invisible(clusterEvalQ(cl, source("reproducible/leverage_funs.R")))
invisible(clusterEvalQ(cl, library(glmmTMB, lib.loc = glmmTMB_lib)))
              
m_big_Lvals2 <- parLapply(cl = cl, eps_vec2,
                      \(e) leverage_brute(m_big, data = df, inds = 1, eps = e,
                                          opt_args = list(control = list(trace = 100))))

## all of these leverages are NEGATIVE (impossible???)
## am I doing something fundamentally wrong?
## is this happening because the original model didn't converge properly?

## test with eps = 0; is it -Inf?
leverage_brute(m_big, data = df, inds = 1, eps = 0, return_delta = TRUE) ## -0.4

## even simpler: refit from (supposed) optimal fit
p0 <- with(m_big$obj$env, parList(last.par.best[-random]))
p0 <- p0[lengths(p0) > 0]
p0 <- p0[setdiff(names(p0), "b")]  ## drop 'b' parameters

## full model refit (slow because we have to compute Std Devs etc ...)

detach("package:glmmTMB")
library(glmmTMB) ## we don't need the fancy leverage one, want parallelization
m_big$obj$fn()
logLik(m_big)
## eigen(vcov(m_big, full = TRUE))$values
m_big_u <- update(m_big, start = p0,
                  control = glmmTMBControl(parallel = list(autopar = TRUE, n = 10),
                                           conv_check = "skip",
                                           optCtrl = list(eval.max=100000, iter.max = 10000, trace = 100)))

m_big$obj$fn() - m_big_u$obj$fn()
