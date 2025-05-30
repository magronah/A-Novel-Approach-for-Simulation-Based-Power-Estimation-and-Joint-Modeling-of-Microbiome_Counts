form = count ~1 + offset(normalizer) +
    (pregnant*GA_Days|taxon) +  ## overall slopes by taxon
    (1|nugget) +   ## same as taxon:Subect_ID:GA_days
    (1|Subect_ID:taxon) +
    ## could instead/additionally do random slopes by taxon within subject:
    ##  (1 + GA_days | Subect_ID:taxon)
    ## similar to:
    ##  ar1(0 + GA_days | Subect_ID:taxon)
    rr(taxon + 0 | Subect_ID:GA_Days,2)
## we know that pregnant is a between-subjects factor, so we have
## already included the only term for that that's possible
## (in the pregnant*GA_days|taxon term)
## we also can rule out any random effect that *doesn't* include taxon
##  in the grouping structure, because that would imply something
##  that changed all taxa at once
## the big complication that hurts my head, and which we haven't
## included or thought about much, is whether there is room for
## more rr() terms with taxon: i.e. rr(taxon + 0| Subect_ID) ?
## rr(0 + taxon:GA_days | Subect_ID) ? or something somehow combining
## these terms (separable covariance structures?)

## suppose we modeled rr(taxon | subject) + rr(taxon:GA_days | subject)
## ( 0+ terms omitted for brevity)
## these are separate RE terms, so there is no correlation between elements
##  of each of these terms. How would we think about relaxing that?
##
## the analogy that I'm trying to think about is when we have spatial
## correlation C1_(x1, x2) and temporal correlation C2_(t1, t2) we often
## think about *separable* correlation, i.e.
## C({x1, t1}, {x2, t2}) = C1_(x1,x2) * C2_(t1, t2)

library(reformulas)
## find RE terms and make them into character
re_terms <- sprintf("(%s)", sapply(findbars(form), deparse))

## want all possible combinations (1, 2, 3, 4)
## sum(choose(4, i))
model_vec <- character(0)  ## not including '0 terms'

for (i in 1:4) {
    model_vec <- c(model_vec,
                   combn(re_terms, i, paste, collapse = "+"))
}

cat(model_vec, sep = "\n")

## nobars() drops random effects; [[3]] extracts the right-hand side
fixed <- deparse(nobars(form)[[3]])

all_models <- sapply(model_vec,
         function(x) reformulate(c(fixed, x), response = "count"))
names(all_models) <- NULL ## names are ugly

all_fits <- lapply(all_models,
                   function(f) glmmTMB(f, ...)
                   )

sapply(all_fits, AIC)

## backward stepwise is a little more complicated ...


## if we have:
(pregnant*GA_Days|taxon) +  ## (1)
    (1|nugget) +                ## (2)
    (1|Subect_ID:taxon) +       ## (3)
    rr(taxon + 0 | Subect_ID:GA_Days,2)  + ## (4)
    ar1([something])             ## (5)

## if one term is singular, drop it
## if the model 'works' (no warnings, no NaNs) keep it
## keep 1 and 4 no matter what ?
## drop 5, then 3, then 2 ... ???

library(glmmTMB)

vars <- setdiff(names(Salamanders),
                c("site", "cover", "sample", "count", "spp"))
sapply(vars, function(x) isNested(Salamanders[[x]], Salamanders$site))

zipm3 <- glmmTMB(count~1  + (1 |site) + (0 + Wtemp | site),
                 zi=~ 1 + (1|spp) + (1 + Wtemp | site), Salamanders, family="poisson")

performance::check_singularity(zipm3, tol = 1e-3)
performance::check_singularity(zipm3, tol = 1e-16)

## performance:::check_singularity.glmmTMB
##' @param x a glmmTMB fit
##' @param tolerance tolerance for checking equality of determinant to zero
##' @param retval: "all" for a test of whether *any* component is singular;
##' "terms" for a separate test of every random effect term in the model
my_sing <- function (x, tolerance = 1e-06, retval = c("all", "terms")) {

    retval <- match.arg(retval)
    res <- list()
    vv <- VarCorr(x)
    for (component in c("cond", "zi", "disp")) {
        nm <- sapply(findbars(formula(x, component = component)), deparse)
        res[[component]] <- sapply(vv[[component]], function(x) det(x) < tolerance)
        names(res[[component]]) <- nm
    }
    if (retval == "overall") return(any(unlist(res)))
    return(res)
}

my_sing(zipm3, retval = "terms")

## when you get a singular fit,
## makes sense _either_ to drop the singular terms
## _or_ to drop the terms you care least about and try again
