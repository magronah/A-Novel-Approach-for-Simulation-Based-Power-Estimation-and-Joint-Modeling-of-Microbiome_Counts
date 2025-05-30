library(glmmTMB)
dd <- readRDS("real_data/atlass_data/results/uszi_each_mod.rds")

logLik(dd)
fn <- dd$obj$fn
ee <- dd$obj$env
pp <- ee$last.par.best[-ee$random]
f0 <- fn()
f1 <- fn(pp)
## these are **not** identical (env contains mutable state ...)
stopifnot(all.equal(f0, f1))
## 23868.37

## set tolerance to 0 to report difference
all.equal(f0, f1, tolerance = 0)  ## diff 2.03e-13

## dd$obj$hessian = ee$hessian is FALSE
## (probably means we didn't compute it already)
## full hessian, probably not what we want
H <- ee$spHess()

## default finite-diff eps for optimHess is 0.001; for H2, 1e-4
system.time(H1 <- optimHess(pp, fn, dd$obj$gr))  ## 11 seconds
system.time(H2 <- numDeriv::jacobian(dd$obj$gr, pp, method = "simple"))

## skip this by default -- slow!
if (FALSE) {
    system.time(H3 <- numDeriv::jacobian(dd$obj$gr, pp,
                                         method = "Richardson"))
    ## 39 seconds
}
## close enough ...
all.equal(unname(H1), H2, tolerance = 1e-4)
all.equal(H2, H3, tolerance = 1e-4)  

## non-pos-def in any case
min(eigen(H1)$values)
min(eigen(H2)$values)
min(eigen(H3)$values)

coef(summary(dd))
e1 <- eigen(H1)
## smallest value is *negative* ... and very wide range, anyway
e1$values

## nothing stands out as being of especially large magnitude ...
## can't pinpoint what's wrong here
par(las = 1, bty = "l")
badvec <- e1$vector[,e1$values<=0]
cvec <- c("black", "red")[as.numeric(badvec>0)+1]
plot(abs(badvec), log = "y", type = "h", col = cvec)
points(abs(badvec), pch = 16, col = cvec)

## brute-force fix
H1f_tmp <- Matrix::nearPD(H1)
H1f_tmp$eigenvalues
H1f_tmp$normF  ## Frobenius norm
## (L2 norm of elementwise differences between original matrix
## and 'pos-defified' matrix); I can't tell whether 2.18 is a large
## number or not.

H1f <- as.matrix(H1f_tmp$mat)
sdvec <- sqrt(diag(solve(H1f)))
sdvec0 <- sqrt(diag(solve(H1)))

sdvec1 <- sdvec0
sdvec1[is.nan(sdvec1)] <- 1000
MASS::eqscplot(log10(sdvec1), log10(sdvec),
               xlab = "log10(original Hessian)",
               ylab = "log10(PD-ified Hessian)",
               col = c("black","red")[as.numeric(is.nan(sdvec0))+1],
               main = "SDs of top-level params (beta, theta)",
               cex=3)
abline(a=0, b=1, lty = 2)

sdr <- TMB::sdreport(dd$obj)
## weird. Giving me non-pos-def Hessian warning, but still gets
## sdreport internally calls optimHess in the same way that I do above
## try(chol(hessian.fixed)) gives "the leading minor of order 24 is not positive"
## the reason is that (on my computer) the Hessian is *negative* definite,
## not singular, so we can invert it -- but, as shown above, the results
## aren't necessarily correct ...

try(chol(H1))
ss <- solve(H1)
ee <- eigen(H1)$values
any(abs(ee)<1e-5)  ## no eigenvalues near zero
any(ee<0)          ## but some are negative!
ee[ee<0]

pp <- predict(dd, type = "latent", se.fit = TRUE)
ord <- order(pp$fit)
with(pp, matplot(seq_along(fit), cbind(fit[ord], fit[ord]-se.fit[ord], fit[ord]+se.fit[ord]), type = "p",
     pch = 16, col = c(1,2,2)))
