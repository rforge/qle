## Quasi-likelihood simulation based estimation
##
## 1. use criterion `qle` for estimation
## 2. then `mahal` as (generalized) least squares
library(qle)
data(normal)

# setting number of local cores
options(mc.cores=6)

## one step minimization
## no sammpling
# x0 <- c(2.5,1.5)
# qscoring(qsd, x0, pl=10,verbose=TRUE)

## alternatively use minimization by `nloptr`
# searchMinimizer(x0, qsd, method = c("bobyqa"), verbose=TRUE)

# main estimation with new evaluations
# (simulations of the statistical model)
OPT <- qle(qsd,qsd$simfn,nsim=10,
		global.opts=list("maxeval"=10),
		pl=10)

# restart estimation and do a pure global search (setting `ftol_abs=0`), 
# sample additional points for evaluation and select new candidates by
# criterion `var`
GL <- qle(OPT$qsd, qsd$simfn, nsim=10,		
		global.opts = list("maxiter"=10, "stopval"=0),
		local.opts = list("nextSample"="var","ftol_abs"=0),
		pl=10, iseed=1234)
