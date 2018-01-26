## Fit a Matern-Cluster point pattern 
## to the `redwood` data from package spatstat

library(qle)
library(spatstat)

# set options
options(mc.cores=8)
RNGkind("L'Ecuyer-CMRG")
set.seed(297)

simStat <- function(X,cond){
 x <- Kest(X,r=cond$rr,correction="best")
 x <- x[[attr(x,"valu")]]
 x <- x[x>0]
 if(anyNA(x) || any(!is.finite(x))) {
  warning(.makeMessage("`NA`, `NaN` or `Inf` detected.","\n"))
  x <- x[!is.nan(x) & is.finite(x)]}
 return(c(intensity(X),x))	
}

# simulation function
simClust <- function(theta,cond){
 X <- rMatClust(theta["kappa"],theta["R"],theta["mu"],win=cond$win)	
 simStat(X,cond)
}

# plot diagnostics 
plotGraphs <- function(par, nsim = 100) {
	# fit by MC
	fit0.mc <- kppm(redwood, ~1, "MatClust")
	# fit by QLE
	names(par)<-c("kappa","scale","mu")
	fit0.ql <- kppm(redwood,~1,"MatClust",
			improve.type="none",startpar=par[1:2],
			control=list(maxit=0),algorithm="SANN")
	fit0.ql$Fit$mcfit$mu <- fit0.ql$mu <- fit0.ql$modelpar[3]<-par[3]
		
	# plotting
	oldpar <- par(no.readonly = TRUE)
	par(mfrow=c(3,2))	
	plot(envelope(fit0.mc, Kest, nsim = nsim), main="Minimum Contrast")
	plot(envelope(fit0.ql, Kest, nsim = nsim), main="QL estimation")
	plot(envelope(fit0.mc, Gest, nsim = nsim), main="Minimum Contrast")
	plot(envelope(fit0.ql, Gest, nsim = nsim), main="QL estimation")
	plot(envelope(fit0.mc, Fest, nsim = nsim), main="Minimum Contrast")
	plot(envelope(fit0.ql, Fest, nsim = nsim), main="QL estimation")
	par(oldpar)
	
	# statistics
	cat("Fitted by quasi-likelihood: \n\n")
	print(coef(summary(fit0.ql)))
	cat("\n\n")
	cat("Fitted by minimum contrast: \n\n")
	print(coef(summary(fit0.mc)))	
}

# load example data set (spatstat)
data(redwood)

# observation window
win <- owin(c(0, 2), c(0, 2))

# condition object: further options
# needed for the simulation function 
cond <- list(win=win,rr=seq(0,0.3,by=0.05)) 

# quasi-likelihood options for estimation
nsim <- 50
Nsample <- 12

# define parameter space
lb <- c("kappa"=20,"R"=0.01,"mu"=1)
ub <- c("kappa"=35,"R"=0.25,"mu"=5)

# general approach to initialize a (local) cluster object
cl <- makeCluster(8)
clusterSetRNGStream(cl)
clusterCall(cl,fun=function(x) library("spatstat", character.only=TRUE))
clusterExport(cl=cl,varlist=c("simStat","searchMinimizer"), envir=environment())

# simulate design points and statistics
sim <- simQLdata(sim=simClust,cond=cond,nsim=nsim,
		method="randomLHS",lb=lb,ub=ub,N=Nsample,cl=cl)

# check simulations used
attr(sim,"nsim")

# generated random design
X <- attr(sim,"X")

# observed statistics (redwood data)
obs0 <- simStat(redwood,cond)

# set up QL model with kriging approximation of
# covariance matrix estimate with bootstrappingg option
qsd <- getQLmodel(sim,lb,ub,obs0,		
		var.type="kriging",			# kriging variance matrix
		intrinsic=TRUE, Nb = 100,	# use bootstrap option, number of bootstrap samples
		verbose=TRUE)

qsd$covT
qsd$covL

# cross-validation: fitting CV covariance models
cvm <- prefitCV(qsd, reduce=FALSE, verbose=TRUE)

# starting point for local search
x0 <- c("kappa"=24,"R"=0.08,"mu"=2.5)

# use the maximum of kriging and CV-based variances
attr(cvm,"type") <- "max"

# first try quasi-scoring with CV errors (type=max)
(QS0 <- qscoring(qsd,x0,
		 opts=list("pl"=10,
				   "ftol_rel"=1e-7,
				   "xtol_rel"=1e-7,
				   "score_tol"=1e-7,
				   "slope_tol"=1e-7),
		 cvm=cvm,pl=10,verbose=TRUE))
 
QS0$Qnorm
quasiDeviance(QS0$par,qsd,cvm=cvm,verbose=TRUE)[[1]]

# multistart version of finding a root
# ... passed to searchMinimizer
method <- c("qscoring","bobyqa","cobyla")
opts <-list("ftol_stop"=1e-10,"ftol_rel"=1e-7,
		    "xtol_rel"=1e-10,"score_tol"=1e-3,
			"slope_tol"=1e-6) 

S0 <- multiSearch(xstart=x0,qsd,method,opts,check=FALSE,
		cvm=cvm,nstart=25,optInfo=TRUE,cl=cl,verbose=TRUE)

# best found root
roots <- attr(S0,"roots")
(id <- attr(roots,"id"))
stopifnot(!is.na(id))
attr(roots,"par")
roots

## inspect CV errors vs. kriging variances
# no significant bias in predicting the statistics 
crossValTx(qsd, cvm, type = "acve")

# compare magnitudes of predictions:
# here: first statistic (intensity) is more sensitive to
# leave out a single sample point of the initial design
# than the others.
crossValTx(qsd, cvm, type = "mse")

# adequacy of the prediction models
crossValTx(qsd, cvm, type = "ascve")
# T2,T4,T5 -> kriging variance underestimates the actual prediction error (measured by CV error)
# T1,T3,T6,T7 -> the actual CV error seems to not sufficiently reflect the predicted
# error by the kriging variance. Strategy: use the maximum of prediction errors

# compute the kriging variance at the sample points i=1,...,n
# leaving out the ith each time 
crossValTx(qsd, cvm, type = "sigK")

## could compare these to the kriging variance:
# dx <- attr(qsd$qldata,"xdim")
# T <- qsd$qldata[(dx+1):(dx+length(qsd$covT))]
# varKM(qsd$covT,X,X,T)

# start main estimation using selection
# criterion `score` (see vignette) and
# the maximum of CV errors and kriging variances
# in order to accouont for the prediction uncertainty
# of sample means of the statistics

OPT <- qle(qsd, simClust, cond=cond,  
		qscore.opts=list("pl"=0,					# >=10 show quasi-scoring iterations
				         "xtol_rel"=1e-7,
						 "ftol_rel"=1e-7,
						 "ftol_abs"=1e-6,
						 "score_tol"=1e-3),
		global.opts = list("maxiter"=5,
				           "maxeval" = 10,
				           "weights"=c(50,10,5,1,0.1),
						   "NmaxQI"=3, 
						   "nstart"=25),			# default number of starting points
		local.opts = list("lam_max"=1e-2,
				          "nobs"=50,				# number of (bootstrap) observations for testing local minimizer
				          "nextSample"="score",		# sample criterion
				          "ftol_abs"=1e-7,			# upper bound on criterion value
						  "weights"=c(0.55),		# constant weight factor
						  "eta"=c(0.025,0.075),	    # ignored, automatic adjustment of weights
						  "test"=FALSE),			# testing is enabled
		method = c("qscoring","bobyqa","direct"),		
		errType="max", iseed=297, multistart=TRUE,
		cl=cl, pl=10)								# pl=10 also show failed local minimizations

print(OPT,pl=10)

# extract information of parameter estimation
local <- OPT$final
info <- attr(OPT,"optInfo")
track <- attr(OPT,"tracklist")
# last message from local minimization
local$message
# history of roots
do.call(rbind,lapply(track,function(x) x$S0$score))

# do a global search with final QL model
# and compare with the following local results
xt <- track[[5]]$Snext$par
D0 <- searchMinimizer(xt, OPT$qsd,
		opts=list("pl"=10,"xtol_rel"=1e-6),
		method="qscoring",cvm=OPT$cvm,
		verbose=TRUE)


S0 <- searchMinimizer(OPT$par, OPT$qsd,
			method="qscoring",cvm=OPT$cvm,
			verbose=TRUE)

# quas-scoring again more precise results
QS <- qscoring(OPT$qsd,OPT$par,
		opts=list("score_tol"=1e-3),
		cvm=OPT$cvm,pl=10)

# compare the different estimates
checkMultRoot(OPT,par=rbind("QS"=QS$par,"S0"=S0$par))

# MC hypothesis testing 
Stest <- qleTest(OPT,QS,sim=simClust,cond=cond,nsim=200,
		  method=c("qscoring","bobyqa","direct"),
		  cl=cl, verbose=TRUE)

print(Stest)

# show envelopes for K,G,F function
plotGraphs(OPT$par,nsim=1000)

## fit model by Minimum Contrast
## and compare with QLE
#data(redwood)
#fitM <- kppm(redwood, ~1, "MatClust")
#fitM$modelpar
#QS$par

# do not forget
stopCluster(cl)


# ------------------------- ONLY FOR THE VIGNETTE ---------------------------

## save results for vignette
#matclust <- list("qsd"=qsd,"cvm"=cvm,"OPT"=OPT,"Stest"=Stest)
#save(matclust,file="matclust.rda")

## plot and store envelopes
#pdf("Kfunc.pdf",width = 8, height = 10)
#plotGraphs(OPT$par,nsim=1000)
#dev.off()
