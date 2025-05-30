source("func2.R")
r = readRDS("/home/michael/Documents/reduced_rank/100_300/GAM/rr.rds")

object = (r$fit_2d)
plot(object)

n  =  50
newdata   = data.frame(lmean_count = rep(5,n),
                       abs_lfc     = seq(1,5,length=n))


library(ggplot2)
pp = power_pred(object, newdata)
ggplot(pp, aes(x=abs_lfc, y = (power))) +
  geom_point() +
  geom_line()  +
  theme_bw()




