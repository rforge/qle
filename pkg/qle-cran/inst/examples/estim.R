# (1) Kriging prediction of sample means of statistics
# (2) Estimation of variance matrix by sample average approximation
library(qle)
data(normal)

# points to estimate statistics
p <- multiDimLHS(N=1000,qsd$lower,qsd$upper,
		method="randomLHS",type="list")

# kriging mean
X <- as.matrix(qsd$qldata[,1:2])
# values of statistics
Tstat <- qsd$qldata[grep("^mean[.]",names(qsd$qldata))] 
# Cholesky decompostions of variance matrices
Lstat <- qsd$qldata[grep("L+",names(qsd$qldata))]		   

# kriging prediction (low level functions)
(est.m1 <- estim(qsd$covT,p,X,Tstat,krig.type="var"))
(est.m2 <- estim(qsd$covT[[1]],p,X,Tstat[1],krig.type="var"))

stopifnot(
  all.equal(as.numeric(extract(est.m1,type="mean")[,1]),
			as.numeric(extract(est.m2,type="mean")))
)

#for(i in 1:2) qsd$covT[[i]]$fix.nugget <- c(qsd$covT[[i]]$fix.nugget,1e-4) 
#M <- varKM(qsd$covT,p,X2,Tstat)
#QD0 <- quasiDeviance(c(2,1),qsd)[[1]]
#L <- mclapply(1:NROW(M), function(i) det(QD0$jac%*%solve(attr(QD0,"Sigma")+diag(M[i,]))%*%t(QD0$jac)))
#L <- mclapply(1:NROW(M), function(i) det(QD0$jac%*%solve(diag(M[i,]))%*%t(QD0$jac)))
#which.max(rowSums(M))
#which.max(unlist(L))
#(QD0$jac)%*%solve(attr(QD0,"Sigma"))%*%t(QD0$jac)
#QD0$varS
#L2 <- mclapply(1:NROW(M), function(i) det(diag(M[i,])))
#which.max(unlist(L2))

# estimate derivatives
jacobian(qsd$covT,p,X,Tstat)

# average variance matrix interpolation
covarTx(qsd,theta=X)

# predict and extract values 
predictKM(qsd$covT,p,X,Tstat,krig.type="var")

# calculate kriging prediction variances
varKM(qsd$covT,p,X2,Tstat)

# update with new found root as minimum of QD
#qsd <- updateQLdata(list(S0),qsd,fit=TRUE,verbose=verbose)																																	
## check results of kriging
#if(.isError(qsd)){
#	msg <- .makeMessage("Cannot update QL data model at approximater root.")				 			    
#	err <- attr(qsd,"error")
#	if(!is.null(err) && inherits(err,"error"))
#		msg <- c(msg, conditionMessage(err))
#	message(msg)
#	break
#}
#X <- as.matrix(qsd$qldata[seq(xdim)])
### TODO: funktioniert cvm nach 'updateQLdata'?											
#QD <- quasiDeviance(Y,qsd,NULL,cvm=cvm,check=FALSE,W=I,theta=xt,
#		criterion="qle",cl=if(isTRUE(use.cluster)) cl)
#if(.isError(QD))
#	stop("Failed to compute quasi-deviance in 'traceQIsamp'.")
#ok <- which(sapply(QD, function(x) !(.isError(x) || anyNA(unlist(x)))))	
#if(length(ok) == 0L)					
#	stop("Quasi-deviance computations have errors.")				
#
#RES <-
#		do.call(doInParallel,
#				c(list(X=QD,
#								FUN=function(qd,qsd,qd0) {					
#									options(mc.cores=1L)
#									qsd <- updateQLdata(list(qd),qsd,fit=FALSE,cl=NULL)
#									if(.isError(qsd))
#										return(qsd)	
#									# Information-based criterion													
#									QDnew <- quasiDeviance(qd0$par,qsd,check=FALSE)[[1]]
#									if(.isError(QDnew))
#										return(.qleError(message="eval error",call=match.call(),error=QDnew))												 	
#									try(det(QDnew$varS),silent=TRUE)												 
#								}, cl=cl, qsd=qsd, qd0=S0)))					
#
#ok <- which(sapply(RES,function(x) !.isError(x) ))			
#if(length(ok) == 0L){
#	msg <- paste("All evaluations of minimum prediction variances failed.")
#	message(msg)
#	return(.qleError(message="Errors in evaluations",call=match.call(),error=RES))	   
#} else if(length(ok) < length(RES)){
#	msg <- .makeMessage("A total of ",length(ok), " minimum prediction variances failed. Check attribute `optRes`.")
#	message(msg)
#}											
#which.max(unlist(RES[ok]))	