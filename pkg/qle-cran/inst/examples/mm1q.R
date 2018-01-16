# TODO: Add comment
# 
# Author: baaske
###############################################################################
library(qle)
library(graphics)

RNGkind("L'Ecuyer-CMRG")
set.seed(1356)

cond <- list("n"=25)
simfn <- function(tet,cond){
	mean(rgeom(cond$n,prob=1-tet[1]))
}

lb <- c("rho"=0.05)
ub <- c("rho"=0.95)

nsim <- 10
X <- multiDimLHS(N=10,lb=lb,ub=ub,
		method="maximinLHS",type="matrix")

sim <- simQLdata(sim=simfn,cond=cond,nsim=nsim,X=X)

qsd <- getQLmodel(sim, lb, ub, obs=c("N"=1),
		var.type="wlogMean",verbose=TRUE)

S0 <- qscoring(qsd,x0=c("rho"=0.8),
		opts=list("pl"=10),verbose=TRUE)

x <- S0$par 
quasiDeviance(x,qsd,verbose=TRUE)[[1]]

# opt
OPT <- qle(qsd,simfn,cond=cond,	     	
		global.opts = list("maxeval"=5, "NmaxLam"=5),
		local.opts = list("nextSample"="score","weights"=0.5,"ftol_abs"=1e-4,
				"lam_max"=1e-5,"test"=TRUE),
		method = c("qscoring","bobyqa","direct"), iseed=1356) 

#plot statistics
op <- par(xaxs='i', yaxs='i')
rho <- as.matrix(seq(0.1,0.9,by=0.001))
y <- as.numeric(unlist(simQLdata(sim=simfn,cond=cond,nsim=nsim,X=rho,mode="mean")))
T <- qsd$qldata[grep("mean.",names(qsd$qldata))]
Y <- predictKM(qsd$covT,rho,X,T,krig.type="var")
# steady state values
y0 <- rho/(1-rho)
plot(NULL, type="n", xlab=expression(rho),
		ylab="y",xlim=c(0,1), ylim=c(0,10))
lines(as.numeric(rho),y,col="black",lt=2,lwd=0.3)
lines(as.numeric(rho),Y,col="blue",lwd=0.3)
lines(as.numeric(rho),y0,col="red",lwd=0.3)
legend("topleft", c("Number of customers in the system",
				"Expected number at steady state","Kriging approximation"),
		lty=c(2,1,1),col=c("black","red","blue"),
		xpd=TRUE,pt.cex=1,cex=1)
par(op)

# next
op <- par(xaxs='i', yaxs='i')
p <- seq(lb,ub,by=0.0001)
QD <- quasiDeviance(X,qsd,value.only=TRUE)
qd <- quasiDeviance(as.matrix(p),qsd)
y <- sapply(qd,"[[","value")
score <- sapply(qd,"[[","score")
## plot quasi-deviance and quasi-score function
plot(NULL, type="n", xlab=expression(rho),
		ylab="",xlim=c(0,1), ylim=c(-10,50))
abline(h=0,col="gray")
points(X,QD,pch=3,cex=1)
lines(p,score, type='l',col="blue",lwd=1.5) 
lines(p,y,col="black",lwd=0.8)
legend("topleft", c("quasi-deviance","quasi-score","sample points", "approximate root","additional samples"),
		lty=c(1,1),lwd=c(1.5,1.5,NA,NA,NA),pch=c(NA,NA,3,5,8),
		col=c("black","blue","black","magenta","green"),pt.cex=1,cex=1)
points(S0$par,S0$val,col="magenta",pch=5,cex=1)
nmax <- OPT$ctls["maxeval","val"]
X <- as.matrix(qsd$qldata[,1])
Xnew <- OPT$qsd$qldata[(nrow(X)+1):(nrow(X)+nmax),1]
points(cbind(Xnew,0),pch=8,cex=2,col="green")
par(op)

# Quasi-deviance
op <-par(xaxs='i', yaxs='i')
qd <- quasiDeviance(as.matrix(p),OPT$qsd)
y <- sapply(qd,"[[","value")
score <- sapply(qd,"[[","score")
## plot quasi-deviance and quasi-score function
plot(NULL, type="n", xlab=expression(rho),
		ylab="",xlim=c(0,1), ylim=c(-10,50))
abline(h=0,col="gray")
lines(p,score, type='l',col="blue",lwd=1.5) 
lines(p,y,col="black",lwd=0.8)
legend("topleft", c("quasi-deviance","quasi-score","sample points", "QL estimate"),
		lty=c(1,1),lwd=c(1,1,NA,NA,NA),pch=c(NA,NA,3,5,8),
		col=c("black","blue","black","magenta","green"),pt.cex=1,cex=1)

X <- as.matrix(OPT$qsd$qldata[,1])
QD <- quasiDeviance(X,OPT$qsd,value.only=TRUE)
points(X,QD,pch=3,cex=1)
points(OPT$par,OPT$val,col="magenta",pch=5)
par(op)