library(glmmTMB)

## put us() first to make it easier to find params (first in b vector)
spider_rr <- glmmTMB(abund ~ 1 + us(1+moss|Species) + rr(Species + 0|id, d = 3),
                     data = spider_long)

system.time(p1 <- predict(spider_rr, type = "latent", se.fit = TRUE))
system.time({
    ss <- TMB::sdreport(spider_rr$obj, getJointPrecision = TRUE)
    se_vec <- sqrt(diag(solve(ss$jointPrecision)))
})

nspp <- 12
se1 <- p1$se.fit[1:(2*nspp)]
se2 <- split(se_vec, names(se_vec))[["b"]][1:(2*nspp)]
all.equal(se1, se2)
all.equal(se1, se2, tolerance = 0) ## diff 1e-15
