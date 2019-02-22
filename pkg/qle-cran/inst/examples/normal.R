# Example: apply `qle` to normal model with criterion `score`
# The following code is also part of the vignette
library(qle)

## a local cluster
cl <- makeCluster(8L)
clusterSetRNGStream(cl,1234)

## Multicore parallel processing:
# options(qle.multicore="mclapply")
# options(mc.cores=2L) 

simfunc <- function(pars) {	
	x <- rnorm(10,mean=pars["mu"],sd=pars["sigma"])    
	c("T1"=mean(x),"T2"=var(x))	
}

# box contraints defining the parameter space
lb <- c("mu"=0.5,"sigma"=0.1)
ub <- c("mu"=8.0,"sigma"=5.0)	   

## the (unknown) true parameter
theta0 <- c("mu"=2,"sigma"=1)

# simulate model at a minimum of required design points
sim <- simQLdata(sim=simfunc,nsim=10,N=12,
		method="maximinLHS",lb=lb,ub=ub)

# set number of simulations manually
# since otherwise only `nsim` would be used to 
# calculate sample average variance
attr(sim,"nsim") <- 100

# true and error-free observation
obs <- structure(c("T1"=2,"T2"=1), class="simQL")

# construct QL approximation model
qsd <- getQLmodel(sim,lb,ub,obs,var.type="wcholMean")

# quasi scoring first try
QS <- qscoring(qsd, x0=c("mu"=5,"sigma"=3.0),opts=list("pl"=10))
print(QS)

OPT <- qle(qsd,	simfunc, nsim=10, global.opts=list("maxiter"=10,"maxeval"=25),
		local.opts=list("lam_max"=1e-3), pl=2L, cl=cl)

OPT$final
OPT$why
OPT$ctls
