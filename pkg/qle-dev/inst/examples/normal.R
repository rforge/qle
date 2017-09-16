\dontrun{
# A pedagogic example of a simple statistical 
# simulation model using normal random numbers
library(qle)	

# init RNG
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# simulation function
simfunc <- function(pars) {	
    x <- rnorm(10,mean=pars["mu"],sd=pars["sigma"])    
    c("T1"=median(x),"T2"=mad(x))	
}

# box contraints defining the parameter space
lb <- c("mu"=0.5,"sigma"=0.1)
ub <- c("mu"=8.0,"sigma"=5.0)	   

## the (unknown) true parameter
theta0 <- c("mu"=2,"sigma"=1)

# simulate model at a minimum of required design points
nsim <- 10
sim <- simQLdata(sim=simfunc,nsim=nsim,N=8,
             method="maximinLHS",lb=lb,ub=ub,
             mode="list")	 

## here generate some test data
#obs <- simQLdata(simfunc,X=theta0,nsim=1)
obs <- structure(c("T1"=2,"T2"=1), class="simQL")

# construct approximation model
qsd <- getQLmodel(sim,lb,ub,obs,var.type="wcholMean")

###########################################################
#	use kriging approximation of variance matrix
###########################################################

## get the model with kriging approximation of
## variance matrix using fixed local nuggets
## and a starting global nugget for reml estimation	

# qsd <- getQLmodel(sim,lb,ub,obs,var.type="kriging",
#          var.opts=list("var.sim"=1e-4,nugget = 1e-6),
#           criterion="qle")

## Construct QL model manually
# qldata <- setQLdata(sim)
# attributes(qldata)

## fit by REML using SIRF-k covariance
## model for all statistics and prepare for
## kriging variance matrix approximation
# mods <- fitSIRFk(qldata, var.type="kriging",		
#            var.opts=list("var.sim"=0.01),verbose=TRUE)	
# qsd <- QLmodel(qldata,lb,ub,obs,mods,var.type="kriging")
	
###########################################################
#	optional: contour plots
###########################################################

## get a device for plotting
## use plot = TRUE in `qle` to show iterates
x <- seq(qsd$lower[1],qsd$upper[1],by=0.05)
y <- seq(qsd$lower[2],qsd$upper[2],by=0.05)
p <- as.matrix(expand.grid(x,y))
X <- as.matrix(qsd$qldata[,1:2])
Tstat <- qsd$qldata[grep("^mean.",names(qsd$qldata))]
Xp <- quasiDeviance(X,qsd,value.only=TRUE)
D <- quasiDeviance(p,qsd,value.only=TRUE)

# library(rgl)
# z <- matrix(D,ncol=length(y))
# open3d()
# persp3d(x,y,z,col="red", alpha=0.3, axes=T)
# cnt <- contourLines(x,y,z,
#      lev=seq(range(z)[1],range(z)[2],
#		  by=dist(range(z))/100))
# for (i in 1:length(cnt))
#   with(cnt[[i]], lines3d(x, y, level, col="darkred"))
# points3d(X[,1],X[,2], Xp, size=3, col="red",add=TRUE)

## contour plot
dev.new()
z1 <- matrix(D,ncol=length(y))
plot(x = 0, y = 0, type = "n", xlim=range(x), ylim=range(y),xlab = "", ylab = "")
contour(x, y, z1, col = "black", lty = "solid",
	nlevels = 50, add = TRUE,vfont = c("sans serif", "plain"))
try(points(X,pch=23,cex=0.8,bg="black"),silent=TRUE)

###########################################################
# 			 apply quasi-scoring iteration
#			 and main estimation function `qle`
###########################################################

# starting point
x0 <- c("mu"=5,"sigma"=3)
# scoring iteration
S0 <- qscoring(qsd, x0, verbose=TRUE)
points(S0$par[1],S0$par[2],col="magenta",pch=23)

#debug(updateCovModels)
# estimate with sampling new points
OPT <- qle(qsd,
          simfunc,		
          nsim=100,
          global.opts=list("maxeval"=100),
          local.opts=list("lam_max"=1e-4,"weights"=0.5),
          plot=TRUE, pl=5)

# get table of stopping conditions
OPT$ctls
# which one stopped iteration
OPT$ctls[OPT$ctls[,"stop"]>0,]
# last local estimation results
local <- attr(OPT,"final")
info <- attr(OPT,"optInfo")

# fit CV models
cvm0 <- prefitCV(qsd, reduce=FALSE)
# CV prediction variance at `theta0` using initial design
crossValTx(qsd, cvm0, type = "cve")
# compare with kriging variance at left out sample points
crossValTx(qsd, cvm0, type = "sigK")
# check design: Kriging variances seem to
# appropriately model the prediction variance
crossValTx(qsd, cvm0, type = "acve")
# given the current sample, we would prefer
# the kriging prediction variance to CV
crossValTx(qsd, cvm0, type = "ascve")
# compare estimated variance with sample variance
obs0 <- simQLdata(simfunc,X=OPT$par,nsim=1000,mode="matrix")[[1]]

# compare simulated and predicted variance matrix
var(obs0)
attr(local,"Sigma")

# MC hypothesis test of OPT$par
# using `local$method`="qscoring"
Stest <- qleTest(OPT,sim=simfunc,nsim=1000,		
	  	    method=c("qscoring","bobyqa"),
			verbose=TRUE)
print(Stest)

## Perform MC test using Mahalanobis criterion,
## first try to improve estimated parameter with final weighting
## matrix `W` at `theta` from the results of function `qle`
## (because of the weighted variance matrix type of approximation).
# info <- attr(OPT,"optInfo")
# MD <- searchMinimizer(x0,OPT$qsd,method=c("lbfgs","bobyqa"),
#		 W=info$W,theta=info$theta,verbose=TRUE)
# Stest.LS <- qleTest(OPT,MD,sim=simfunc,nsim=1000,		
#			    verbose=TRUE)
# print(Stest.LS)


###########################################################
#	optional: contour plots of final approximation
###########################################################

# nmax <- OPT$ctls["maxeval","val"]
# Xnew <- OPT$qsd$qldata[nrow(X):(nrow(X)+nmax),c(1,2)]
#
# qsd2 <- OPT$qsd
# X <- as.matrix(qsd2$qldata[,1:2])
# Xp <- quasiDeviance(X,qsd2,value.only=TRUE)
# D <- quasiDeviance(p,qsd2,value.only=TRUE)
#
# dev.new()
# z1 <- matrix(D,ncol=length(y))
# plot(x = 0, y = 0, type = "n", xlim=,range(x), ylim=range(y),xlab = "", ylab = "")
# contour(x, y, z1, col = "black", lty = "solid",
#		nlevels = 50, add = TRUE,vfont = c("sans serif", "plain"))
# try(points(X,pch=23,cex=0.8,bg="black"),silent=TRUE)


## save results as data set 
## qsd$sim <- simfunc
## save(qsd,file="normal.rda")
}

