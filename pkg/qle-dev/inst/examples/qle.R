\dontrun{
## Quasi-likelihood simulation based estimation
##
## 1. use criterion `qle` for estimation
## 2. then `mahal` as (generalized) least squares
library(qle)
data(normal)

options(mc.cores=8)

RNGkind("L'Ecuyer-CMRG")
set.seed(1356)

# one step minimization
# no sammpling
x0 <- c(2.5,1.5)
qscoring(qsd, x0, pl=10,verbose=TRUE)

# alternatively use minimization by `nloptr`
searchMinimizer(x0, qsd, method = c("bobyqa"), verbose=TRUE)

# QLE approach:
OPT <- qle(qsd,qsd$sim,nsim=100,
		global.opts=list("maxeval"=50),
		local.opts=list("lam_max"=1e-3,"weights"=0.5),
		pl=3)

# solution (with  details)
print(OPT)

# prepare for a restart 
qsd2 <- OPT$qsd

# and do a pure global search by `ftol_abs=0` 
# sample additional (globally selected) points
# for evaluation by selection criterion `var`
GL <- qle(qsd2, qsd$sim, nsim=100,		
		global.opts = list( "maxiter"=10, "stopval"=0),
		local.opts = list("nextSample"="var","ftol_abs"=0),
		iseed=123, pl = 5)


# Use a least squares approach 
qsd$criterion <- "mahal"
MA <- qle(qsd, qsd$sim, method = c("lbfgs","bobyqa","direct"), 
		global.opts = list("maxiter" = 10, "maxeval"=25),
		iseed=123,pl = 3)

print(MA)
}

