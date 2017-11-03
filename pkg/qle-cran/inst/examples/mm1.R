# TODO: Add comment
# 
# Author: baaske
###############################################################################
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
		#c("N"=tet/(1-tet) + rnorm(1,0,0.25))
		mean(rgeom(25,prob=1-tet[1]))		
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
	
	qsd <- getQLmodel(sim, lb, ub, obs,
			var.type="wlogMean", verbose=TRUE)
	
#qsd <- getQLmodel(sim, lb, ub, obs,
#		var.type="kriging",
#		var.opts = list("var.sim"=1e-3,"nugget"=1e-2), verbose=TRUE)
	
#qsd <- getQLmodel(sim, lb, ub, obs,
#		var.type="kriging", verbose=TRUE)
	
#qsd <- getQLmodel(sim, lb, ub, obs, var.type="kriging",
#		 set.var=FALSE, var.sim=1e-07, as.nugget=TRUE,
#		 var.opts=list("var.sim"=0.001), verbose=TRUE)
# 
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
	
#op <-par(mfrow=c(1, 2),
#		cex=1.7, cex.axis=1.7, cex.lab=1.7,lwd=0.5,
#		cex.main=1.7, cex.sub=1.7,
#		xaxs='i', yaxs='i')
	
	## kriging approximation
	op<-par(mfrow=c(1,2))
	plot(NULL, type="n", xlab=expression(rho),
			ylab="T",xlim=c(0,1), ylim=c(0,10))
	
	lines(as.numeric(rho),y,col="black",lt=2)
	lines(as.numeric(rho),Y,col="blue")
	lines(as.numeric(rho),y0,col="red")
#legend("topleft", c("Number of customers in the system",
#				    "Expected number at steady state",
#					"Kriging approximation"),
#			cex=2.2, lty=c(2,1,1),col=c("black","red","blue"))
	
# quasi-deviance plots
	p <- seq(0,1,by=0.0001)
	QD <- quasiDeviance(X,qsd,value.only=TRUE)
	qd <- quasiDeviance(as.matrix(p),qsd)
	y <- sapply(qd,"[[","value")
	score <- sapply(qd,"[[","score")
	
	## plot quasi-deviance and quasi-score function
	plot(NULL, type="n", xlab=expression(rho),
			ylab="quasi-deviance",xlim=c(0,1), ylim=c(-10,50))
	abline(h=0)
	points(X,QD,pch=3)
	lines(p,score, type='l',col="blue",lwd=1.5) 
	lines(p,y,col="black",lwd=0.8)
#legend("top", c("quasi-deviance","quasi-score","sample points", "solution"),
#		lty=c(1,1),lwd=c(1.5,1.5,NA,NA),pch=c(NA,NA,3,5),cex=2.2,
#		col=c("black","blue","black","magenta"))
	par(op)
	
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
			method=c("qscoring"),verbose=TRUE)
	
	points(S0$par,S0$val,col="magenta",pch=5)
	
# start sampling estimation
#debug(qle)
	OPT <- qle(qsd,
			simfn, 	     	
			global.opts = list("maxeval"=5,"NmaxLam"=5),
			local.opts = list("nextSample"="score","weights"=0.5,"ftol_abs"=1e-4,
					"lam_max"=1e-5,"useWeights"=TRUE),
			method = c("qscoring","bobyqa","direct"),
			plot=TRUE, pl = 100, iseed=1356) 
	
	print(OPT)
	
#debug(checkMultRoot)
	checkMultRoot(OPT,par=S0$par)
	
	Stest <- qleTest(OPT,sim=simfn,nsim=100,
			#method=c("qscoring","bobyqa","direct"),
			iseed=1234, verbose=TRUE)
	
	
# final criterion function results
	OPT$final
	
	# prediction of mean number of customers based on final approximation
	QD <- quasiDeviance(OPT$par,OPT$qsd,verbose=TRUE)[[1]]
	X <- as.matrix(OPT$qsd$qldata[,1])
# values of statistics
	Tstat <- OPT$qsd$qldata[grep("mean.",names(qsd$qldata))]   
	predictKM(OPT$qsd$covT,tet0,X,Tstat)
# variance
	
# variance geom random number
	T <- 25
	simfn2 <- function(tet) rgeom(T,prob=1-tet[1])
	simV <- simQLdata(sim=simfn2,nsim=nsim,X=tet0)
	vars <- var(unlist(simV))
	vars
	tet0/(1-tet0)^2
	
# variance sample mean of geom random number
	T <- 25
	simfn2 <- function(tet) mean(rgeom(T,prob=1-tet[1]))
	simV <- simQLdata(sim=simfn2,nsim=1000,X=tet0)
	vars <- var(unlist(simV))
	vars
	tet0/(1-tet0)^2/T   # variance of sample mean
	
# estimated (average) variance by QLE) vs. exact variance 
	attr(QD,"Sigma")
	(tet <- OPT$par)
	(tet/(1-tet)^2)/T
	sqrt((tet/(1-tet)^2)/T)
	
	QD$score
	ET <- (1-tet0)/(1-tet0)^2
	sF <- function(tet) { 
		x <- rgeom(250,prob=1-tet)
		n*(1/(1-tet) - mean(x)/(tet))
	}
	y <- sF(p)
	lines(p[-c(1,length(p))],y[-c(1,length(p))],type='l',col="magenta",lwd=1.5)
	
	## MLE!
	n <- 25*100 # nur Anz. Simulation des Anfangs-Samples
	x <- rgeom(n,prob=1-tet0)
	y <- sum(x)
	(tet.mle <- y/(n+y)) #1-1/(1+mean(x))
	
			
# or
	1-1/(1+mean(x))
# score
	n*(1/(1-tet.mle) - mean(x)/(tet.mle))	
	(tet.mle*(1-tet.mle)^2)/n
	
# variance N
	(tet/(1-tet)^2)/n
	tet/((1-tet)^2*n)
# rho mle
	(tet0*(1-tet0)^2)/n
	
# small simulation study
	study <- function(tet0,n){
		y <- sum(rgeom(n,prob=1-tet0))
		y/(n+y)
	}
	D <- lapply(1:100,
			function(i) {
				y <- sum(rgeom(n,prob=1-tet0))
				y/(n+y)
			})
	sd(unlist(D)) # -> emp. error MLE
	.MSE(matrix(unlist(D)),tet0)^0.5
	
	sqrt(1/Imle)  # lower bound error: stimmz noch nicht ganz!!!
	Stest		  # QL estimation error: emp. and predicted.
	
# abs deviations
	abs(tet0-tet)
	abs(tet0-tet.mle)
	
# Vgl. estimates and Fisher QI und I (MLE)
	
# Fisher information/quasi-information
	(Imle <- n/((1-tet.mle)^2*tet.mle))
	1/Imle
	1/OPT$final$I > 1/Imle
	OPT$final$I < Imle
		
	tet
	sqrt(1/OPT$final$I)
	sqrt(1/QD$Iobs)
	
	tet.mle
	sqrt(1/Imle)
	sqrt(1/Iobs.mle)
	
	Iobs.mle <- n/(tet.mle*(1-tet.mle)^2)
# or
	(Iobs.mle <- n*(1/(1-tet.mle)^2 + mean(x)/tet.mle^2))
	
# variance rho
	1/Iobs.mle
# or
	(tet.mle*(1-tet.mle)^2)/n
# compare 
	1/QD$Iobs
	1/QD$I
	QD$varS
	
#n <- 1000
#x <- rgeom(n,prob=tet0)
#y <- sum(x)
#tet <- 1/(1+mean(x))
#n*(1/tet - mean(x)/(1-tet))
	
	n*tet/(1-tet)
	n*tet/(1-tet)^2
	
	## plotting
	dev.new()
	op <-par(cex=1.7, cex.axis=1.7, cex.lab=1.7,lwd=0.5,
			cex.main=1.7, cex.sub=1.7,
			xaxs='i', yaxs='i')
	qd <- quasiDeviance(as.matrix(p),OPT$qsd)
	y <- sapply(qd,"[[","value")
	score <- sapply(qd,"[[","score")
	## plot quasi-deviance and quasi-score function
	plot(NULL, type="n", xlab=expression(rho),
			ylab="quasi-deviance",xlim=c(lb,ub), ylim=c(-10,50))
	abline(h=0)
	lines(p,score, type='l',col="blue",lwd=1.5) 
	lines(p,y,col="black",lwd=0.8)
	X <- as.matrix(OPT$qsd$qldata[,1])
	QD <- quasiDeviance(X,OPT$qsd,value.only=TRUE)
	points(X,QD,pch=3,cex=1)
	points(OPT$par,OPT$val,col="magenta",pch=5)
	legend("topleft", c("quasi-deviance","quasi-score","sample points", "QL estimate"),
			lty=c(1,1),lwd=c(1.5,1.5,NA,NA,NA),pch=c(NA,NA,3,5,8),
			col=c("black","blue","black","magenta","green"),pt.cex=1.6,cex=1.6)
	par(op)
	
	
	
#	% TODO: Stest, MLE
#			%	- cond no hinzufügen cond <- list("n"=25)
#			%	-Prediction of mean number of customers based on final approximation
#			%   -Estimated (average) variance by QLE) vs. exact variance 
#			%   -MLE for rho! based on all nsim*N simulations of initial sample set
#			% 	-Compare RMSE (QLE, Stest) with RMSE from simulation study for (MLE) 
#			
	