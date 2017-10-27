## 1D example of QLE:
## Simulate M/M/1 queue using its (noisy) steady
## state distribution, estimate the parameter `rho` by
## the mean number of customers as a statistic 
library(qle)
# setting RNG and seed
RNGkind("L'Ecuyer-CMRG")
set.seed(1356)

## steady state function
# success parameter: 1-rho
# equivalently use `rgeom`
simfn <- function(tet) {
	c(tet/(1-tet) + rnorm(1,0,0.25))
	# mean(rgeom(100,prob=1-tet[1]))		
}

# parameter space
lb <- c("rho"=0.05)
ub <- c("rho"=0.95)
# number of simulations
nsim <- 10	 
# initial sample size
Nsample <- 10

# simulate initial design
X <- multiDimLHS(Nsample,lb,ub,
      method="maximinLHS",type="matrix")

## simulate at design X
sim <- simQLdata(sim=simfn,nsim=nsim,X=X)

## the observed data:
## mean number of customers
tet0 <- c("rho"=0.5)
# ok, this equals 1 anyway
obs <- c("N"=tet0/(1-tet0) )

## check with exact statistic
## simulated and exact mean number of customers
# obs0 <- simQLdata(sim=simfn,X=tet0,nsim=10000,mode="mean")[[1]]
# rbind("obs0"=mean(unlist(obs0)),obs)

## Use kriging of variance of statistic
## (nugget variance fixed to 0.1% of Cholesky decomposed values).
## Note that we use the default `controls`
## parameter although we could choose
## different settings for optimization.
qsd <- getQLmodel(sim, lb, ub, obs, var.type="kriging",
		 set.var=FALSE, var.sim=1e-07, as.nugget=TRUE,
		 var.opts=list("var.sim"=0.001), verbose=TRUE)
 
## powexp
#qsd <- getQLmodel(sim, lb, ub, obs, var.type="kriging",
#		 nfit=1, set.var=FALSE, model="powexp",
#		 var.opts=list("var.sim"=0.001), verbose=TRUE)

## matern
#qsd <- getQLmodel(sim, lb, ub, obs, var.type="kriging",
#		 nfit=1, set.var=TRUE, model="matern", 
#		 var.opts=list("var.sim"=0.001), verbose=TRUE)
 
## sirfk with average approximation of variance 
## no estimate of nugget (fixed), no simulation variance as local nugget
#qsd <- getQLmodel(sim, lb, ub, obs, var.type="wcholMean",trend=2,
#		 set.var=FALSE,var.sim=0,nugget=1e-7,fixed.param="nugget")

## sirfk
## use simulation variances as local nuggets
#qsd <- getQLmodel(sim, lb, ub, obs,
#		var.type="wcholMean",trend=2,nugget=1e-6)
  
### plot statistics
rho <- as.matrix(seq(0.1,0.9,by=0.001))
y <- as.numeric(unlist(simQLdata(sim=simfn,nsim=10,X=rho,mode="mean")))
T <- qsd$qldata[grep("mean.",names(qsd$qldata))]
Y <- predictKM(qsd$covT,rho,X,T,krig.type="var")
# steady state values
y0 <- rho/(1-rho)


## --------------------------------------------------------------------------------------------------
## used only to generate pdf output
# pdf("mm1q.pdf",width = 50, height = 25)
# op <-par(mfrow=c(1, 2), mar=c(5.1, 5.1, 1.1, 1.1),
#		oma=c(5,4,1,1), cex=2.2, cex.axis=2.2, cex.lab=2.2,lwd=0.5,
#		cex.main=2.2, cex.sub=2.2, xaxs='i', yaxs='i')
#----------------------------------------------------------------------------------------------------

## kriging approximation
plot(NULL, type="n", xlab=expression(rho),
		ylab="T",xlim=c(0,1), ylim=c(0,10))

lines(as.numeric(rho),y,col="black",lt=2)
lines(as.numeric(rho),Y,col="blue")
lines(as.numeric(rho),y0,col="red")
legend("topleft", c("Number of customers in the system",
				    "Expected number at steady state",
					"Kriging approximation"),
			cex=2.2, lty=c(2,1,1),col=c("black","red","blue"))

# quasi-deviance plots
p <- seq(lb[1],ub[1],by=0.0001)
QD <- quasiDeviance(X,qsd,value.only=TRUE)
qd <- quasiDeviance(as.matrix(p),qsd)
y <- sapply(qd,"[[","value")
score <- sapply(qd,"[[","score")

## plot quasi-deviance and quasi-score function
plot(NULL, type="n", xlab=expression(rho),
      ylab="quasi-deviance",xlim=c(lb,ub), ylim=c(-10,50))
abline(h=0)
points(X,QD,pch=3)
lines(p,score, type='l',col="blue",lwd=1.5) 
lines(p,y,col="black",lwd=0.8)
legend("top", c("quasi-deviance","quasi-score","sample points", "solution"),
		lty=c(1,1),lwd=c(1.5,1.5,NA,NA),pch=c(NA,NA,3,5),cex=2.2,
		col=c("black","blue","black","magenta"))

# par(op)
# dev.off()

# some options for quasi-scoring
opts <- list("pl"=0, # change to > 1 for output 
		"ftol_stop"=1e-10,
		"slope_tol"=1e-5,
		"score_tol"=1e-6)

## scoring iteration
x0 <- c("rho"=0.25)
# qscoring(qsd,x0,opts)
S0 <- searchMinimizer(x0, qsd, opts=opts,
       method=c("qscoring","bobyqa"), verbose=TRUE)
points(S0$par,S0$val,col="magenta",pch=5)

# check
(QD <- quasiDeviance(S0$par,qsd)[[1]])

# start sampling estimation
OPT <- qle(qsd,
		   simfn,		     	
		   global.opts = list("maxiter" = 20, "NmaxLam"=5,"maxeval"=25),
		   local.opts = list("nextSample"="score","weights"=0.5,"ftol_abs"=1e-4,
				             "lam_max"=1e-5,"useWeights"=FALSE,"eta"=c(0.01,0.1)),
		   method = c("qscoring","bobyqa","direct"),
		   plot=TRUE, pl = 100) 

print(OPT)

# final criterion function results
OPT$final
