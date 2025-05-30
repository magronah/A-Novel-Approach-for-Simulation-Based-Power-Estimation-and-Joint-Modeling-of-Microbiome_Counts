library(tidyverse)
library(glmmTMB)
library(patchwork)
library(purrr)
library(foreach)
library(doParallel)
library(here)
#setwd(here())
theme_set(theme_bw())

#########################################################
custom_theme <- function(n) {
  theme_bw(base_size = n) +
    theme(
      plot.title = element_text(hjust = 0.5,size = n, family = "Roboto",color = "black"),
      text = element_text(size = n, family = "Roboto",color = "black"),
      axis.text.x = element_text(family = "Roboto", size = n, color = "black"),
      axis.text.y = element_text(family = "Roboto", size = n, color = "black")
    )
}

make_pars <- function(pars, ...) {
  ## FIXME: check for name matches, length matches etc.
  L <- list(...)
  for (nm in names(L)) {
    pars[names(pars) == nm] <- L[[nm]]
  }
  return(pars)
}


sfun <-function(object, data, pars, ..., show_pars = FALSE)  {
  form <- object
  form[[3]] <- form[[2]]
  form[[2]] <- quote(..y)
  data[["..y"]] <- if (!identical(list(...)$family, "beta_family"))
    1
  else 0.5
  r1 <- glmmTMB(form, data = data, ..., doFit = FALSE)
  r2 <- fitTMB(r1, doOptim = FALSE)
  return(r2)
}

SimulateFunc <- function(object, nsim = 1, seed = NULL, data, pars, ..., show_pars = FALSE)  {
  if (!is.null(seed)) set.seed(seed)
  r2 <- sfun(object, data, pars, ...)
  if (show_pars)
    return(r2$env$last.par)
  pars <- do.call("make_pars", c(list(r2$env$last.par), pars))
  replicate(nsim, r2$simulate(par = pars), simplify = FALSE)
}
#########################################################
set.seed(101)
ntaxa = 12; nIndiv = 10; ntime = 20
d <- expand.grid(taxa = factor(1:ntaxa),subject = factor(rep(1:nIndiv,ntaxa)))[1:(ntaxa*nIndiv), ]
expdes <- data.frame(subject = factor(1:nIndiv),
                     group=rep(c("control","treatment"), each = nIndiv/2))
## BMB: left_join(0 avoids warnings
dat0 <- left_join(d, expdes, by = "subject")

dd <- do.call("rbind", replicate(n=ntime, dat0, simplify = FALSE))
dd$time = rep(1:ntime, each=ntaxa*nIndiv)
dd$nugget <- factor(1:nrow(dd))
###############################################
################################
form <- ~ar1(factor(time) + 0 |taxa:subject) +
  diag(group*time|taxa) + #deviations of each OTU from fixed Effect
  (1|nugget) +
  (taxa + 0 | subject:time) #  We are assumming that the correlation structure is the same for every subjec

##
beta.param = rep(0,4)
l = length(beta.param)

set.seed(101)
thet = c(
  log(0.1), 0.7/sqrt(1-0.7^2),        # ar1: 2 sd, phi
  log(0.2),log(0.3),log(0.1), log(0.1), # sd: group * times | OTU:
  rnorm(1),                    # nugget
  log(runif(ntaxa)),
  rnorm(ntaxa*(ntaxa + 1)/2)   # covariance: OTU | subject:
)

simulate.data <- function(beta.param,dat,thet) {
  SimulateFunc(form,
               data = dat,
               family = poisson,
               pars = list(
                 theta = thet)
               
  )[[1]]
}

sim=simulate.data(beta.param,dd,thet)
dd$count <- sim$yobs

max(sim$yobs)

dd0= dd[dd$taxa ==1,]
dd0 = (dd0[dd0$subject == c("1","6"),])


nn  =4; lw =1.5
dd0= dd[dd$taxa ==1,]
dd0 = (dd0[dd0$subject == c("1","6"),])
p1 <- dd0 %>% 
  ggplot(aes(x=time, y=count,  colour= group,
             group = group, linetype = group),size = nn)+
  geom_line(linewidth = lw) +
  geom_smooth(method = lm, se=FALSE, linewidth = lw) +
  geom_line(linewidth = 2)+ geom_point(size = nn) + theme_bw() + 
  custom_theme(17) +
  labs(x="time points",y = "count abundance")

getwd()
width =  10; height =  6; dpi = 300
ggsave("fig/Request.png", plot = p1, 
       width = width, 
       height = height, 
       dpi = dpi)

