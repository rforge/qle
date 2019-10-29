## Copyright (C) 2018 Markus Baaske. All Rights Reserved.
# This code is published under the GPL (>=3).
#
# File: 	qsOpt.R
# Date:  	14/03/2018
# Author: 	Markus Baaske
#
# Functions for simulated quasi-likelihood estimation,
# sampling candidates for evaluation, cross-validation,
# implements local and global search methods,
# quasi-scoring method 

.qleError <- function(subclass = NULL,
		        message = "", call = as.list(sys.call(-1))[[1]],
		 		error=NULL, ...) {
	
	 args <- list(...)
	 if(length(args) > 0L && any(is.null(names(args))))
	   stop("Additional arguments must have names.")
	 structure(
		list(message = .makeMessage(message), call = call),
	 	 ...,
	 class = c(subclass, c("error","condition")), error = error)
}

.isError <- function(x) {	
	( "error" %in% class(x) ||
	  !is.null(attr(x,"error")) || isTRUE(attr(x,"error")) ||
	  inherits(x, "error") || inherits(x, "try-error") 
	)
}

# internal function to check the arguments `args` with of function `fun` 
.checkfun <- function(fun, args, hide.args = NULL, remove = NULL, check.default=TRUE) {	
	funname <- deparse(substitute(fun)) 
	if( !is.function(fun) )
	  stop(paste0(funname, " must be a function\n"))
	flist <- formals(fun)
	# remove `...`
	if("..." %in% names(formals(fun))) {
	  flist <- flist[-which(names(formals(fun))=="...")]
  	  if("..." %in% names(args))
	    args <- args[-which(names(args)=="...")]
	}	
    if ( length(flist) > 0L ) {		
		fnms  <- if(!is.null(remove)) names(flist)[-remove] else names(flist) 
		rnms  <- names(args)				
		m1 <- match(fnms, rnms)
		if(length(m1) == 0L && length(rnms) > 0L) {
		  for(i in 1:length(rnms)) {
			stop(paste0("Argument `",rnms, "` passed but not required in function `",
				funname,"`.\n"))
			}
		} 
		if(anyNA(m1)) {
			mx1 <- which(is.na(m1))			
			if(!check.default) 
				mx1 <- mx1[!(mx1 %in% which(nzchar(flist)))]
			if(length(mx1) > 0L && !is.null(hide.args))
		      mx1 <- mx1[-which(mx1==pmatch(hide.args,fnms))]
		    if(length(mx1) > 0L) {
			 for( i in 1:length(mx1)){
				stop(paste0(funname, " requires argument `",
						fnms[mx1[i]], "` which has not been passed.\n"))
			 }
		    }
		}
		m2 <- match(rnms, fnms)
		if(anyNA(m2)){
			mx2 <- which(is.na(m2))
			for( i in 1:length(mx2)){
				stop(paste0("Argument `",rnms[mx2[i]], "` passed but not required in function `",
						funname,"`.\n"))
			}
		}
	}
	return(0)
}

.checkOptions <- function(optlist, opts) {
	if(is.null(names(opts)))
		stop("Options should be a list of named arguments.")
	if (!is.list(opts) || "" %in% names(opts))
		stop("Argument 'opts' must be a list of named (character) elents.")
	optnames <- (names(opts) %in% names(optlist))
	if (!all(optnames)) {
		unames <- as.list(names(opts)[!(optnames)])
		stop(paste(c("Unknown arguments in 'opts': ",do.call("paste", c(unames, sep = ", "))), collapse=" "))
	}
	return (0)
}

.qsOpts <- function(options = list(), xdim = 1L) {
	opts <- .addQscoreOptions(xdim)	
	if(length(options) > 0L) {
	 .checkOptions(opts,options)
	 namc <- match.arg(names(options), choices = names(opts), several.ok = TRUE)
	 if (!is.null(namc))		 
		opts[namc] <- options[namc]
	}
	# invert scaling constants
	txid <- which(opts$xscale != 1)
	if(length(txid)>0L)
	 opts$xscale[txid] <- 1/opts$xscale[txid] 
	tfid <- which(opts$fscale != 1)
	if(length(tfid)>0L)
	 opts$fscale[tfid] <- 1/opts$fscale[tfid]
 
	return(opts)
}

.addQscoreOptions <- function(xdim) {
	list( "ftol_stop" = .Machine$double.eps,	  
		  "step_tol"  = 1e-13,
		  "xtol_rel"  = 1e-11,								# see also steptol (Dennis & Schnabel)
		  "grad_tol"  = 1e-4,  
		  "ftol_abs"  = 1e-6,								# used if stepmin or grad_tol reached 
		  "ltol_rel"  = 1e-4,								# relative step length tolerance
		  "score_tol" = 1e-5,								# also used to select best roots		 
		  "slope_tol" = 1e-9,
		  "maxiter"   = 100,
		  "xscale" 	  = rep(1,xdim),						# scaling independent variables, e.i. parameter theta
		  "fscale" 	  = rep(1,xdim),						# scaling quasi-score components for 0.5*norm^2 of quasi-score only 
		  "restart"	  = TRUE,
		  "pl" = 0L)
}

.getDefaultoptions <- function(xdim) {
	list("stopval"=0.0,"lam_rel"=1e-2,"lam_max"=1e-3,"xtol_rel"=1e-3,"ftol_rel" = 1e-4,
		 "ftol_abs"= 1e-6,"var_qd"=0.1,"var_qd_abs"=1e-2,"sampleTol" = 1e-4,"maxiter" = 100,"pmin"=0.05,
		 "weights" = c(1,2,4,6,12,16,20,50), "nsample" = 1e3*xdim,"npsample" = 1e2*xdim,"ntotal"=(1e2*xdim)-10L,
		 "dtol"=c(0.5,0.05), "useWeights"=FALSE, "eta" = c(1,2),"nfail" = 1, "nsucc" = 1, "tau"=2,
		 "Nstopval"=1,"Nxtol_rel" = 3,"Nftol_rel" = 5,"Nlam_rel" = 3,"NsampleTol" = 3,"Nlam_max" = 3,"Nvar_qd" = 2,
		 "nstart" = 10*(xdim+1L), "ntest"=2, "doStest"=TRUE, "alpha"=0.05)		
}

.setControls <- function(optlist) {		
	defaults <- c("lam_rel","lam_max", "var_qd", "xtol_rel", "ftol_rel",
			"stopval", "sampleTol", "maxiter", "nfail","nsucc")   
	namc <- match.arg(names(optlist), choices = defaults, several.ok = TRUE)
	ctls <- data.frame(
			 cbind("val" = 0,	
				   "cond" = unlist(optlist[namc]),
				   "tmp" = 0,				   
				   "stop" = 0,
  				   "count" = c(unlist(optlist[c(paste0("N",namc))]),1,1,1)),
			 row.names = namc, check.names = FALSE)	
	# init some controls	
	ctls["sampleTol","val"] <- Inf
	ctls[c("lam_rel","lam_max","var_qd"),"val"] <- rep(1,3)	
	return (ctls)
}

#' @name getDefaultOptions
#' 
#' @title Print default options for optimization
#' 
#' @description Print default options for global and local optimization routines
#' for function \code{\link{qle}}
#' 
#' @param xdim 		dimension of the unknown model parameter
#' 
#' @return List of options.
#' 
#' @details The function returns a list of available options
#'  for functions \code{\link{qscoring}} and \code{\link{qle}}.
#' 
#' @examples
#' getDefaultOptions(xdim=2)
#'  
#' @author M. Baaske
#' @rdname getDefaultOptions
#' @export
getDefaultOptions <- function(xdim) {
	if(!is.numeric(xdim))
	  stop("`xdim` mus be a numeric value.")  
	list("qscoring" = .addQscoreOptions(xdim),		 
		 "qle_opts" = .getDefaultoptions(xdim))		 
}

## Internal, not exported
## 'points' is best chosen as matrix
cverrorTx <- function(points, Xs, dataT, cvm, Y, type, useMax = FALSE, cl = NULL) {	
	# kriging predictions of full models	
	dfx <- as.data.frame(extract(Y,type="mean"))
	# kriging variances of full models
	dfs2 <- if(useMax || type == "cvmix") 
	 { as.data.frame(extract(Y,type="sigma2")) } else NULL 
	
    np <- nrow(dfx)	 # number of fitted cov models equals 
	n <- length(cvm) # number of blocks for jackknife variance
		
	# prediction function for CV
	.statsCV <- function(mods,points,Xs,dataT,type) {				
		id <- attr(mods$covT,"id")
		if(type == "cvmix") {
		 .COL2LIST(varKM(mods$covT,points,Xs[-id,,drop=FALSE],dataT[-id,]))				
		} else .COL2LIST(predictKM(mods$covT,points,Xs[-id,,drop=FALSE],dataT[-id,],"var"))		
	} 	
	L <- tryCatch(
			doInParallel(cvm, .statsCV, points=points, Xs=Xs, dataT=dataT, type=type, cl=cl)
			 ,error = function(e) {
				msg <- .makeMessage("Cross-validation prediction failed: ",
						conditionMessage(e))
				message(msg)
				.qleError(message=msg,call=match.call(),error=e)
			} )
	# on error
	if(.isError(L))
	 return(L)
	sigL <- NULL
	if(type == "nmsd" || type=="sig2LOO") {
		.varCV <- function(mods,points,Xs,dataT,type) {
			id <- attr(mods$covT,"id")
			.COL2LIST(varKM(mods$covT,points,Xs[-id,,drop=FALSE],dataT[-id,]))
		}
		# kriging variances of LOO-CV models
		sigL <- tryCatch(
			doInParallel(cvm, .varCV, points=points, Xs=Xs, dataT=dataT, type=type, cl=cl)
			,error = function(e) {
				msg <- .makeMessage("Cross-validation prediction failed: ",
						conditionMessage(e))
				message(msg)
				.qleError(message=msg,call=match.call(),error=e)
			} )
	}
	if(.isError(sigL))
	 return(sigL)
	
	do.call(cbind,
		lapply(1:length(L[[1]]), function(i) {
			##  index i (i=1,...,p) is over statistics
			##  index k is over prefitted covariance models with exactly (n-1) points
			y <- do.call(cbind,lapply(L,function(mod) mod[[i]]))					
			switch(type,
				"cve" =	{
					m.cvjack <- n*dfx[,i]-(n-1)*y
					cv <- apply(m.cvjack,1,
							function(x) {
								xn <- mean(x)
								sum((x-xn)^2)/((n-1)*n) 
							}
					)
					if(useMax) {
						pmax(cv,dfs2[,i])
					} else cv		
				},
				"scve" = { 										# standardized (by kriging variance) cross-validation error				 
					if(n!=np)
						stop("Standardized jackknife variance calculation only if number of samples equals number of models.")				   	
					sigK <- sapply(1:length(cvm),
							function(k) {
								id <- attr(cvm[[k]]$covT,"id")
								varKM(cvm[[k]]$covT[[i]],points[k,], Xs[-id,,drop=FALSE],dataT[-id,i])								 
							}
					)	 		
					(dfx[,i]-diag(y))^2/sigK 		  	 
				},			 
				"msd" = { rowMeans((y - dfx[,i])^2) },				# CV based mean squared deviation (prediction uncertainty)
				"nmsd" = {
					# extract kriging variances of LOO-CV 
					sig2 <- do.call(cbind,lapply(sigL,function(mod) mod[[i]]))
					rowMeans(((y-dfx[,i])^2)/sig2)									
				},
				"sig2LOO" = {
					sig2 <- do.call(cbind,lapply(sigL,function(mod) mod[[i]]))	
					rowMeans(sig2)
				},
				"aae" = { rowMeans(abs(y - dfx[,i])) },
				"rmsd" = { sqrt(rowMeans((y - dfx[,i])^2)) },		# CV based root mean squared deviation (prediction uncertainty)
				"acve" = {											# average CV errors at sample points (model validity) to assess the bias in estimation
					if(n!=np)
						stop("Average cross-validation error calculation only available if the number of sample points equals the number of CV models.")			  	 
					# should be approximately zero for all statisitcs i=1,...,p
					# for no systematic over- or under estimation of statistics
					mean(diag(y)-dfx[,i])
				},
				"mse" = {											# average mean squared CV errors at sample points (model validity)
					if(n!=np)
						stop("Cross-validation MSE calculation can be computed only if number of samples equals the number of CV models.")		  	 	
					mean((diag(y)-dfx[,i])^2)
				},
				"ascve" = {											# standardized CV based mse by kriging variances 
					if(n!=np)									    # at sample points (model validity)
						stop("Standardized cross-validation MSE can be computed only if number of samples equals number of CV models.")			     	 
					sigK <- sapply(1:length(cvm),
							function(k) {
								id <- attr(cvm[[k]]$covT,"id")
								varKM(cvm[[k]]$covT[[i]],points[k,], Xs[-id,,drop=FALSE],dataT[-id,i])										
							})	 
					mean((dfx[,i]-diag(y))^2/sigK) 			  	 
				},
				"sigK" = { 										# standardized (by kriging variance) cross-validation error				 
					if(n!=np)
						stop("Leave-one-out kriging variance is only available if the number of sample points equals the number of CV models.")				   	
					sapply(1:length(cvm),
							function(k) {								
								id <- attr(cvm[[k]]$covT,"id")
								varKM(cvm[[k]]$covT[[i]],points[k,], Xs[-id,,drop=FALSE],dataT[-id,i])								 
							}
					)				 		  	 
				},
				"cvmix" = {
					# full model kriging variances are stored in dfs2[,i]
					try(sqrt(rowMeans(sqrt(y))*sqrt(dfs2[,i])),silent=TRUE)
					#try(rowMeans(sqrt(y)),silent=TRUE)
					#try(rowMeans(sqrt(dfs2[,i])),silent=TRUE)
				}
			)		
		}))	
}

#' @name crossValTx 
#'
#' @title Prediction variances by cross-validation
#'
#' @description The function estimates the prediction variances by a cross-validation approach (see vignette)
#'  applied to each sample means of summary statistics. 
#'
#' @param qsd   	object of class \code{\link{QLmodel}}
#' @param cvm		list of prefitted covariance models obtained from function \code{\link{prefitCV}}
#' @param theta		optional, default \code{NULL}, list or matrix of points where to estimate prediction variances
#' @param type		name of type of prediction variance 
#' @param cl	    cluster object, \code{NULL} (default), of class \code{MPIcluster}, \code{SOCKcluster}, \code{cluster}
#' 						 
#' 	
#' @return A matrix of estimated prediction variances for each point given by the argument \code{theta} (as rows)
#'  and for each statistic (as columns).  
#'
#' @details	Other than the kriging prediction variance, which solely depends on interdistances of sample points
#'  and estimated covariance parameters of some assumed to be known spatial covariance model, the cross-validation
#'  based approach (see [4] and the vignette) even takes into account the predicted values at `\code{theta}` and thus can be thought of
#'  a more robust measure of variability between different spatial locations. By default, `\code{theta}` equals the current sampling set 
#'  stored in the object `\code{qsd}`.
#' 
#'  If we set the type of measure `\code{type}` equal to "\code{cve}", the impact on the level of accuracy (predicting at unsampled
#'  points) is measured by a \emph{delete-k jackknifed variance} of prediction errors. This approach does not require further
#'  simulation runs as a measure of uncertainty for predicting the sample means of statistics at new candidate points accross the parameter space.
#'  If \code{attr(cvm,"type")} equals "\code{max}", then the maximum of kriging and CV-based prediction variances is returned. 
#' 
#'  In addition, other measures of prediction uncertainty are available such as the \emph{root mean square deviation}
#'  (\code{rmsd}) and \emph{mean square deviation} (\code{msd}) or the \emph{standardized cross-validation error}
#'  (\code{scve}). The details are explained in the vignette. In order to assess the predictive quality of possibly
#'  different covariance models (also depending on the initial sample size), including the comparison of different
#'  sizes of initial sampling designs, the following measures [8] are available for covariance model validation and adapted
#'  to the cross-validation approach here by using an \emph{average cross-validation error} (\code{acve}), the \emph{mean square error} (\code{mse})
#'  or the \emph{average standardized cross-validation error} (\code{ascve}). These last measures can only be computed in case the total number
#'  of sample points equals the number of leave-one-out covariance models. This requires to fit each cross-validation
#'  covariance model by \code{\link{prefitCV}} using the option `\code{reduce}`=\code{FALSE} which is then based on exactly
#'  one left-out point. Also, we can calculate the kriging variance at the left-out sample points if we set the option `\code{type}`
#'  equal to "\code{sigK}". 
#'
#' @examples
#' data(normal)
#' 
#' # design matrix and statistics
#' X <- as.matrix(qsd$qldata[,1:2])
#' Tstat <- qsd$qldata[grep("^mean[.]",names(qsd$qldata))]
#' 
#' # construct but do not re-estimate
#' # covariance parameters by REML for CV models
#' qsd$cv.fit <- FALSE
#' cvm <- prefitCV(qsd)
#' theta0 <- c("mu"=2,"sd"=1)
#' 
#' # get mean squared deviation using cross-validation at theta0 
#' crossValTx(qsd, cvm, theta0, type = "msd")
#' 
#' # and kriging variance  
#' varKM(qsd$covT,theta0,X,Tstat) 	 
#' 
#' 
#' @seealso \code{\link{prefitCV}}
#'
#' @author M. Baaske
#' @rdname crossValTx
#' @export
crossValTx <- function(qsd, cvm, theta = NULL, 
		type = c("rmsd","msd","cve","scve","acve","mse","ascve","sigK","cvmix","nmsd","aae","sig2LOO"),
		cl = NULL)
{		
	type <- match.arg(type)
	stopifnot(class(cvm) %in% c("cv","cvfull"))	
	dx <- attr(qsd$qldata,"xdim")
	Xs <- as.matrix(qsd$qldata[seq(dx)])
	
	# set sample points as default points to predict the CV error
	if(is.null(theta)) theta <- Xs
	# dataT has to be list (of class data.frame)
	dataT <- qsd$qldata[(dx+1):(dx+length(qsd$covT))]
	
	tryCatch({
		Y <- estim(qsd$covT,theta,Xs,dataT,krig.type="var")		
		# cross-validation variance/RMSE of statistics
		cv <- cverrorTx(theta,Xs,dataT,cvm,Y,type,FALSE,cl)
		structure(cv,dimnames=list(NULL,names(dataT)))			
	 }, error = function(e) {
	 		msg <- .makeMessage("Could not calculate cross-validation variance: ",
					conditionMessage(e))
			message(msg)
			return(.qleError(message=msg,call=match.call(),error=e))
		}
	)
}

## Internal
## COMMENT: 
##	  i is numeric vector of indices of left out points
updateCV <- function(i, qsd, fit, N, opts, stats.only = TRUE) {	
	xdim <- attr(qsd$qldata,"xdim")
	Xs <- as.matrix(qsd$qldata[seq(xdim)])
	stid <- (xdim+1):(xdim+length(qsd$covT))
	fitit <- (fit && !( (N-length(i)) %% qsd$nfit))
	
	tryCatch({
		 ## CV update of statistics
		 attr(qsd$covT,"id") <- i
		 for(j in 1:length(qsd$covT)) {			 					
			 qsd$covT[[j]]$start <- qsd$covT[[j]]$param[qsd$covT[[j]]$free]				
			 if(!is.null(qsd$covT[[j]]$fix.nugget))
				 qsd$covT[[j]]$fix.nugget <- qsd$covT[[j]]$fix.nugget[-i]
			 if(fitit){ 
				 qsd$covT[[j]]$dataT <- as.numeric(qsd$qldata[-i,stid][[j]])			  
			  res <- doREMLfit(qsd$covT[[j]],Xs[-i,],opts)
			  if(!inherits(res,"error")) {
				qsd$covT[[j]]$param <- res$param    
			  } else {		
				msg <- message("Could not update covariance parameters because `REML` failed.")				
				attr(qsd$covT[[j]],"error") <- .qleError(message=msg,error=res,id=i)				  				
			  } 
		  	}			
		 } 
		 ## CV update of variance kriging models 
		 if(qsd$var.type %in% c("kriging","logKrig") && !stats.only) {
		   if(is.null(qsd$covL)) 
		    stop("Covariance models for kriging the variance matrix not set.")
		   attr(qsd$covL,"id") <- i   
		   for(j in 1:length(qsd$covL)) {			 					
			  qsd$covL[[j]]$start <- qsd$covL[[j]]$param[qsd$covL[[j]]$free]				
			  if(!is.null(qsd$covL[[j]]$fix.nugget))
				  qsd$covL[[j]]$fix.nugget <- qsd$covL[[j]]$fix.nugget[-i]
			  if(fitit){ 
				  qsd$covL[[j]]$dataT <- qsd$qldata[-i,grep("^L[^b]",names(qsd$qldata))[j]]
				  res <- doREMLfit( qsd$covL[[j]],Xs[-i,],opts)
				  if(!inherits(res,"error")) {
					  qsd$covL[[j]]$param <- res$param	    
				  } else {		
					  msg <- message("Could not update covariance parameters because `REML` failed.")					 
					  attr(qsd$covL[[j]],"error") <- .qleError(message=msg,error=res,id=i)				  				
				  } 
			  }			  
		  }
	  	}
	  }, error = function(e) {	
		  msg <- .makeMessage("Cross-validation update of covariance models failed: ", conditionMessage(e))
		  message(msg)
		  return(.qleError(message=msg,call=match.call(),error=e))
	  })
	  if(stats.only)
	   return( structure(list("covT"=qsd$covT)) )
	 
	  qsd$qldata <- qsd$qldata[-i,]
	  return(qsd)  	  			    	# for CV error of quasi-score	  
}

#' @name prefitCV 
#'
#' @title Covariance parameter estimation for cross-validation 
#'
#' @description The function constructs a list of covariance models of statistics in order to estimate the prediction error
#'  variances by a cross-validation approach at unsampled points. 
#'
#' @param qsd   	  object of class \code{\link{QLmodel}}
#' @param reduce	  if \code{TRUE} (default), reduce the number of covariance models to refit
#' @param type		  type of prediction variances, "\code{cv}" (default) or "\code{cvfull}", see \code{\link{crossValTx}}
#' @param stats.only  logical, \code{TRUE} (default), only return list of updated cross-validation covariance models of statistics 
#' @param return.qsd  logical, \code{FALSE} (default), return list of QL models for cross-validation 
#' @param control	  control arguments for REML estimation, passed to \code{\link[nloptr]{nloptr}}  	
#' @param cl	      cluster object, \code{NULL} (default), of class \code{MPIcluster}, \code{SOCKcluster}, \code{cluster}
#' @param verbose	  if \code{TRUE}, print intermediate output
#'
#' @return A list of certain length depending on the current sample size (number of evaluated points).
#'  Each list element corresponds to a (possibly reduced) number of sample points with at most \eqn{k} points
#'  (see details) left-out for fitting the corresponding covariance models. 
#'
#' @details Using the cross-validation approach (see vignette) for estimating the prediction variances 
#' 	might require a refit of covariance parameters of each statistic based on the remaining sample points.
#'  The covariance models are refitted if `\code{fit}` equals \code{TRUE} and otherwise simply updated without fitting which
#'  saves some computational resources. The number of points left-out, if applicable, is dynamically adjusted depending on the number
#'  of sample points in order to prevent the main estimation algorithm to fit as many models as there are points already evaluated.  
#' 
#'  The number \eqn{n_c} of covariance models still to fit, that is, the number of partitioning sets of sample points, is limited by
#'  \eqn{n_c\leq n}, with maximum \eqn{k} sampling points deleted from the full sample set with overall \eqn{n} sample points such that
#'  \eqn{n=n_c k} (see also the vignette for further details). 
#' 
#' @examples 
#'   data(normal)
#'   
#'   # without re-estimation of covariance parameters, default is TRUE
#'   qsd$cv.fit <- FALSE  
#'   cvm <- prefitCV(qsd)
#'   
#' @seealso \code{\link{QLmodel}}
#' 
#' @author M. Baaske
#' @rdname prefitCV
#' @export
prefitCV <- function(qsd, reduce = TRUE, type = c("cv","cvfull"),
		       control = list(), cl = NULL, verbose = FALSE)
{		
	# do not exclude corner points from CV because of bas extrapolation properties
	#rng <- apply(X, 2, range)
	#id <- as.numeric(sapply(1:ncol(X),function(i) which(X[,i] %in% rng[,i])))
	
	xdim <- attr(qsd$qldata,"xdim")
	X <- as.matrix(qsd$qldata[seq(xdim)])
	N <- nrow(X)
	Ni <- seq_len(N)
	# number of bins for LOOCV grows like square root of number of samples
	nb <-
	 if(reduce) {
		if((N-2^xdim) > qsd$minN){
			pts <- as.matrix(expand.grid(qsd$lower,qsd$upper))
			id <- .min.distXYIdx(X,pts)
			Ni <- as.numeric(rownames(X[-id,,drop=FALSE]))
			N <- length(Ni)
		}
		try(floor(2*qsd$minN+sqrt(max(N-qsd$minN,0))),silent=TRUE) 
	} else N	
	# block size
	k <- if(!is.numeric(nb) || inherits(nb,"try-error")){
		message("Block size of Leave-(k)-out is invalid. Using k=1.")  
		1L
	} else {
		if(nb>0){ ceiling(N/nb)	}
		else { message("Cannot determine block size. Using k=1."); 1L}
	}
	S <-
	 if((N-k) >= qsd$minN) {	 			
		split(Ni, sort(Ni%%nb))
	} else {
	 stop(paste0("Total number of points must be at least ",qsd$minN," for cross-validation."))
	}	
	type <- match.arg(type)	
	fit <- isTRUE(qsd$cv.fit)		
	# Leave-(k)-Out CV of statistics only
	tryCatch({			 
		if(length(control) > 0L) {		
			opts <- nloptr::nl.opts()
			opts[names(control)] <- control
		} else {
			opts <- attr(qsd,"opts")		
		}	
		structure(doInParallel(S, updateCV, qsd=qsd, fit=fit, N=N, opts=opts,
					stats.only=(type == "cv"), cl=cl), class = type)

		## debugging
		#debug(updateCV)
		#res <- updateCV(1,qsd,fit=fit, N=N,opts=opts, stats.only=(type == "cv"),pl=0,verbose=FALSE)

	 },error = function(e) {
		msg <- paste0("Prefitting covariance models failed.\n")
		message(msg)
		return(.qleError(message=msg,call=match.call(),error=e))
	 }
	)	
}

# internal: prepare variance matrix of statistics
varMatrix <- function(qsd, Sigma = NULL, ...) {				
	if(!(qsd$var.type %in% c("kriging","logKrig","const"))){		
		covTx <- covarTx(qsd,...)
		if(.isError(covTx)){
			msg <- paste0("Could not compute variance matrix.")
			message(msg)
			return(.qleError(message=msg,error=covTx))	
		}
		Sigma <- covTx[[1]]$VTX
	} else if(qsd$var.type == "const"){
		if(is.null(Sigma))
			stop("`Sigma` cannot be NULL if used as a constant variance matrix.")
		# Only for constant Sigma, which is used as is and inverted by default
		Sigma <- try(gsiInv(Sigma),silent=TRUE)
		if(inherits(Sigma,"try-error")) {
			msg <- paste0("Inversion of constant variance matrix failed.")
			message(msg)
			return(.qleError(message=msg,error=Sigma))
		}
	}
	return(Sigma)		
}

# internal, alloc C structure
# Also for QL, a pre-set (inverse) variance matrix can be supplied by VTX
# No predictions variances here (done at C level), theta is only needed
# for the weighted version of avergage variance approximation
.qdAlloc <- function(qsd, Sigma = NULL, ..., cvm = NULL) {	
	xdim <- attr(qsd$qldata,"xdim")	
	Sigma <- varMatrix(qsd,Sigma,...,cvm=cvm)
	if(.isError(Sigma)){
		msg <- paste0("Could not compute variance matrix for CV models.")
		message(msg)
		return(.qleError(message=msg,error=Sigma))
	}	
	qlopts <- list("varType"=qsd$var.type,"useCV"=!is.null(cvm)  && class(cvm)=="cv")
	args.qsd <- list("qlm"=qsd,"VTX"=Sigma,"X"=as.matrix(qsd$qldata[seq(xdim)]))
	
	# LOO CV models patch
	if(!is.null(cvm) && class(cvm) == "cvfull") {
		cvm <- doInParallel(cvm,
				function(qlm,...) {				
					Sigma <- varMatrix(qlm,...)
					if(.isError(Sigma))						  
					 return(.qleError(message=.makeMessage("Could not compute variance matrix for CV models."),error=Sigma))						
					list("qlm"=qlm,"VTX"=Sigma,"X"=as.matrix(qlm$qldata[seq(attr(qlm$qldata,"xdim"))]))				
				},Sigma=Sigma,...) 
		hasErr <- which(sapply(cvm,function(x) .isError(x)))
		if(length(hasErr)>0) {
			msg <- .makeMessage("Failed to compute variance matrix for cross-validation models.")
			message(msg)
			return(.qleError(message=msg,error=cvm))
		}		  
	}
	# return TRUE for success otherwise signal error
	try(.Call(C_initQL,args.qsd,qlopts,cvm))	
}

# internal, free memory
.qdDealloc <- function() {
   try(try(.Call(C_finalizeQL),silent=TRUE))	
}

.checkArguments <- function(qsd, x0=NULL, Sigma = NULL, ...) {
	if(class(qsd) != "QLmodel"){
	   stop("`qsd` object must be of class `QLmodel`.")
    }
   	if(!is.null(x0)) {
	    if(!is.numeric(x0) || anyNA(x0))
		  stop("Starting point must be numeric vector.")		
		# bounds checks
		if( length(qsd$lower)!=length(x0) || length(qsd$upper)!=length(x0))
			stop("Length of 'x0' does not match 'lower' or 'upper' bounds length.")	
		if(any(x0<qsd$lower) || any(x0>qsd$upper))
			stop("At least one element in 'x0' does not match bound constraints. Please check!")
	}
	if(!is.null(Sigma)){
		stopifnot(is.matrix(Sigma))			  	  	  
		if(nrow(Sigma)!=length(qsd$covT) )
		 stop("Dimensions of `Sigma` must match the number of statistics.\n")
		if(qsd$var.type == "const" && qsd$criterion == "qle")
			warning("'Sigma' is chosen to be a constant variance matrix while criterion 'qle' is applied.")			
				
	} else if(qsd$var.type == "kriging" && is.null(qsd$covL))
		stop("Covariance models for kriging variance matrix must be given, see function `setQLdata`.")	
	  else if(qsd$var.type == "const") 
		stop("`Sigma` must not be NULL for `const` variance matrix approximation.")
	  	
}

#' @name searchMinimizer
#'
#' @title Minimize a criterion function 
#'
#' @description The function either searches for a root of the quasi-score or minimizes one of the criterion functions.
#'
#' @param x0		  (named) numeric vector, the starting point
#' @param qsd   	  object of class \code{\link{QLmodel}}
#' @param method	  names of possible minimization routines (see details) 
#' @param opts		  list of control arguments for quasi-scoring iteration, see \code{\link{qscoring}}
#' @param control 	  list of control arguments passed to \code{\link[nloptr]{nloptr}} routines defined in \code{method}
#' @param ...		  further arguments passed to \code{\link{covarTx}}, \code{\link{qscoring}} and \code{\link{quasiDeviance}} 
#' @param obs		  numeric vector of observed statistics, overwrites `\code{qsd$obs}`
#' @param optInfo	  logical, \code{FALSE} (default), not yet used argument (ignored)
#' @param check		  logical, \code{TRUE} (default), whether to check input arguments
#' @param minimize    logical, \code{TRUE} (defualt), unless \code{method} is not "\code{qscoring}", the posterior density is maximized 
#' @param restart 	  logical, \code{TRUE} (default), whether to restart optimization in case of non convergence
#' @param pl		  numeric value (>=0), the print level 
#' @param verbose	  if \code{TRUE} (default), print intermediate output
#'
#' @details The function provides an interface to local and global numerical minimization routines
#'  using the approximate quasi-deviance (QD) or Mahalanobis distance (MD) as an objective (monitor) function.
#'  
#'  The function does not require additional simulations to find an approximate minimizer or root of the quasi-score.
#'  The numerical iterations always take place on the fast to evaluate criterion function approximations.
#'  The main purpose is to provide an entry point for minimization without the need of sampling new candidate points for evaluation.
#'  This is particularly useful if we search for a "first-shot" minimizer or to re-iterate a few further steps after estimation of the
#'  model parameter.  
#' 
#'  The criterion function is treated as a deterministic (non-random) function during minimization
#'  (or root finding) whose surface depends on the sample points and chosen covariance models. Because of the typical nonconvex
#'  nature of the criterion functions one cannot expect a global minimizer by applying any local search method like, for example,
#'  the scoring iteration \code{\link{qscoring}}. Therfore, if the quasi-scoring iteration or some other available method gets stuck
#'  in a local minimum of the criterion function showing at least some kind of numerical convergence we use such minimizer as it is and
#'  finish the search, possibly being unlucky, having not found an approximate root of the quasi-score vector (or minimum of the Mahalanobis distance).
#'  If there is no obvious convergence or any error, the search is restarted by switching to the next user supplied minimization routine defined
#'  in the vector of method names `\code{method}`. 
#' 
#'  \subsection{Choice of local minimization routines}{  
#'  Besides the quasi-scoring method, `\code{method}` equal to "\code{qscoring}", the following
#'  (derivative-free) routines from the \code{\link[nloptr]{nloptr}} package are available for minimizing
#'  both criterion functions:
#'  
#' 	\itemize{
#' 	  \item{}{ \code{\link[nloptr]{bobyqa}}, \code{\link[nloptr]{cobyla}} and \code{\link[nloptr]{neldermead}}}
#'    \item{}{ \code{\link[nloptr]{direct}}, global search with a locally biased version named \code{directL}}
#' 	  \item{}{ \code{\link[nloptr]{lbfgs}},  for minimizing the MD with constant `\code{Sigma}` only}
#' 	  \item{}{ \code{\link[nloptr]{nloptr}}, as the general optimizer, which allows to use further methods}
#'  }
#'    
#'  Using quasi-scoring first, which is only valid for minimizing the QD, is always a good idea since we might have done
#'  a good guess already being close to an approximate root. If this fails we switch to any of the above alternative methods
#'  (e.g. \code{\link[nloptr]{bobyqa}} as the default method) or eventually - in some real hard situations - to the
#'  method `\code{direct}`, if given, or its locally biased version `\code{directL}`. The order of processing is determined
#'  by the order of appearance of the names in the argument `\code{method}`. Any method available from package `\code{nloptr}` can be
#'  chosen. In particular, setting \code{method="nloptr"} and defining `\code{control}` in an appropriate way allows to choose a multistart
#'  algorithm such as \code{\link[nloptr]{mlsl}}, see also \code{\link{multiSearch}} for an alternative solution.
#' 
#'  Only if there are reasonable arguments against quasi-scoring, such as expecting local minima rather than a root first or an available
#'  limited computational budget, we can always apply the direct search method `\code{direct}` leading to a globally exhaustive search.
#'  Note that we must always supply a starting point `\code{x0}`, which could be any vector valued parameter of the parameter space unless
#'  method `\code{direct}` is chosen. Then `\code{x0}` is still required but ignored as a starting point since it uses the "center point" of
#'  the (hyper)box constraints internally. In addition, if a list of cross-validation models either `\code{cvm}` of the covariance models of
#'  the statistics or of the full quasi-score are given, then these cross-validation prediction variances
#'  are inherently used instead of kriging prediction variances for estimation of the variance of the quasi-score during consecutive iterations
#'  of all optimization methods. This results in some additional computational effort due to the repeated computations of the statistics or
#'  quasi-score vector, respectively, to calculate these variances during each new iteration.  
#' }
#' 
#' @return A list as follows
#' 	  \item{par}{solution vector}
#' 	  \item{value}{objective value}
#' 	  \item{method}{applied method}
#' 	  \item{convergence}{termination code}
#' 	  \item{score}{if applicable, quasi-score vector (or gradient of MD)}
#' 
#' @examples
#' data(normal)
#' searchMinimizer(c("mu"=2.5,"sd"=0.2),qsd,method=c("qscoring","bobyqa"),verbose=TRUE) 
#' 
#' @seealso \code{\link{prefitCV}}, \code{\link[nloptr]{nloptr}}, \code{\link{qscoring}}, \code{\link{multiSearch}}
#' 			
#' @rdname searchMinimizer
#' @author M. Baaske
#' @export
#' @importFrom nloptr direct directL cobyla bobyqa lbfgs neldermead
searchMinimizer <- function(x0, qsd, method = c("qscoring","bobyqa","direct"),
					 opts = list(), control = list(), ...,  
					   obs = NULL, optInfo = FALSE, check = TRUE, 
					    minimize = TRUE, restart = TRUE, pl = 0L, verbose = FALSE)
{
	args <- list(...)
	x0 <- if(is.matrix(x0))
		structure(as.numeric(x0),names=colnames(x0))	
	else unlist(x0)
	if(check){
	 .checkArguments(qsd,x0,...)
	 stopifnot(is.numeric(pl) && pl >= 0L )
    }  
    fun.name <- ""
	nms <- names(x0)	
	scoring <- isTRUE("qscoring" %in% method)
	# to be sure
	x0 <- .PROJMED(x0,qsd$lower,qsd$upper)
	# current sample points
	xdim <- attr(qsd$qldata,"xdim")
	if(xdim != length(x0))
	 stop("Dimension of starting point 'x0' does not match the problem dimension.")
 	
 	# init tracklist for error tracking
 	tracklist <- structure(list(),class="QDtrack") 
	# may overwrite (observed) statistics	
	if(!is.null(obs)) {
		obs <- unlist(obs)
		if(anyNA(obs) | any(!is.finite(obs)))
			warning("`NA`, `NaN` or `Inf` values detected in argument `obs`.")
		if(!is.numeric(obs) || length(obs)!=length(qsd$covT))
		  stop("Object `obs` must be a (named) numeric vector or list of length equal to the number of given statistics in `qsd`.")
		qsd$obs <- obs
		# TODO
		if(!is.null(args$cvm)){
		 for(i in 1:length(cvm)) args$cvm[[i]]$obs <- obs				
		}
  	} 

	S0 <-
	 if(qsd$criterion != "qle"){
		m1 <- pmatch("qscoring",method)
		if(!is.na(m1)) {
		  method <- method[-m1]
		  message(.makeMessage("Scoring not available for criterion `mahal`, using `",method,"` instead.\n"))
	  	}
	    if(length(method) == 0L)
		  stop("Only a single local search method is specified: ")
	 	fun.name <- method[1]
	    NULL
	 } else {
		fun.name <- "qscoring"
		m1 <- pmatch(fun.name,method)
		if(!is.na(m1)){
		 if(m1!=1)
		  method <- c("qscoring",method[-m1])		
		 tryCatch({			
		    qscoring(qsd,x0,opts,...,check=FALSE,verbose=verbose)			
		   }, error = function(e) {	e }
  		 )
		} else NULL
	}
	
    if(!is.null(S0) && (.isError(S0) || S0$convergence < 0L)){	   
	   msg <- .makeMessage("Minimization by `",fun.name,"` did not converge: ")
	   if(!is.null(S0$convergence))
 	    msg <- c(msg, paste0(" (status = ",S0$convergence,")") )
	   if(inherits(S0,"error"))
		msg <- c(msg, conditionMessage(S0))
	   if(verbose)  message(msg)	   	   
	   method <- method[-1]
	   if(is.na(method[1]) || !restart){
		message(.makeMessage("No convergence or restart required and only one local method supplied."))
		return(S0)	
	   }
	   tracklist <- c(tracklist,list(S0))
    }
	
	if(is.null(S0) || (restart && S0$convergence < 0L)) {	  	
	  S0 <- 
		tryCatch({			
			# allocation at C level
			if(.isError(.qdAlloc(qsd,...)))
			 stop("Allocation error: cannot request C memory.")
			
			fn <- 
			 if(minimize) {
			  switch(qsd$criterion,
				"qle" = function(x) { .Call(C_qDValue,x) },
				"mahal" = { function(x) .Call(C_mahalValue,x) }				 
			  )
			 } else {
			  switch(qsd$criterion,
				"qle" = function(x) { -.Call(C_qDValue,x) },
				"mahal" = { function(x) -.Call(C_mahalValue,x) }	)
			 }
		 	 repeat {
				if(!is.na(method[1])) {
					if(verbose)
					  cat(paste0("Using method: ",method[1],"...\n"))
					fun.name <- method[1]					
				} else {
					return(
					 .qleError(message = "No convergence and only one method supplied: ",
							   call = sys.call(),
							   error = if(inherits(S0,"error")) conditionMessage(S0) else NULL,
					   		S0 = S0, method = method[1], tracklist = tracklist)
					)	
				}				
			 	S0 <-
					tryCatch({
						switch(fun.name,
								"direct" = {
									direct(fn, lower=qsd$lower, upper=qsd$upper, control=control)
								},
								"directL" = {
									directL(fn, lower=qsd$lower, upper=qsd$upper, control=control)
								},
								"lbfgs" = {									
									if(qsd$criterion != "mahal" || qsd$var.type != "const")
									  stop("`lbfgs` only for criterion `mahal` using a constant `Sigma`. Note that, in this case only the quasi-score vector is the gradient of the criterion function 'mahal'.")
									lbfgs(x0,
										  fn = function(x) {
												 val <- fn(x)
										 		 return(
												   list("objective" = val, 
														"gradient" = -attr(val,"score")))
									   	},
										lower=qsd$lower, upper=qsd$upper, control=control)
							   	},
								"nloptr" = {
									fn <- function(x,cvm){									
										val <- .Call(C_qDValue,x)
										as.numeric(qDPress(as.matrix(x),val,cvm,fun="lapply")$mean)										
									} 									
									if(is.null(control$algorithm)){
										control["algorithm"] <- "NLOPT_LN_BOBYQA"
										if(verbose)
										 message(paste0("Using default derivative-free method: ",control$algorithm))									
									}									
									ret <- do.call(nloptr::nloptr,
											list(x0, eval_f=fn, lb=qsd$lower, ub=qsd$upper, cvm=args$cvm, opts=control))
									
									structure(list("par"=ret$solution,
												   "value"=ret$objective,
												   "iter"=ret$iterations,
												   "convergence"=ret$status,	# 0 for success
												   "message"=ret$message))									
								},
								{
									fun <- try(get(fun.name),silent=TRUE)
									if(inherits(fun,"try-error") || !is.function(fun))
									   stop(paste0("Unknown function call: ",fun.name,".\n"))									
								    # default call to `nloptr`
								    do.call(fun, list(x0, fn, lower=qsd$lower, upper=qsd$upper, control=control))								
								}
						)		 
					  }, error = function(e) { e })
			    
				if(!inherits(S0,"error") && S0$convergence >= 0L) {					
					attr(S0,"restarted") <- TRUE
					S0$bounds <- which(S0$par >= qsd$upper | S0$par <= qsd$lower)
					break
				} else {
					msg <- .makeMessage("'Nloptr' failed or no convergence by: ",fun.name,".")
					message(msg, if(inherits(S0,"error")) conditionMessage(S0) else "", sep=" ")
				  	method <- method[-1]
					tracklist <- c(tracklist,list(S0))
				}
			}										
			S0  # success: nloptr			
		}, error = function(e) {
			 msg <- .makeMessage("Restarted 'nloptr' minimization failed: ", conditionMessage(e))
			 message(msg)
			 return(.qleError(message = msg,call = sys.call(), error = e,
					   method = fun.name, tracklist = tracklist))			
		}, finally = { 
			 if(!.qdDealloc())
			   stop("Allocation error in C memory management.")
		})	
	}
	if(.isError(S0))
	  return(S0)		
	if(!is.null(nms))
 	  names(S0$par) <- nms     
 	
    if(class(S0) != "QSResult") {		
	  	if(!minimize && is.numeric(S0$value)) S0$value <- -S0$value
		S0 <- structure(
	    	    c(S0,list("method"=fun.name,				   	  
						  "criterion"=qsd$criterion,						 
				 		  "start"=x0)),
	  		   restarted=attr(S0,"restarted"), class="QSResult")
	
		qd <- tryCatch({				
			    quasiDeviance(S0$par,qsd,...,check=FALSE,verbose=verbose)									
			}, error = function(e) {
				 msg <- .makeMessage("Error in criterion function: ",conditionMessage(e))
				 message(msg)
				.qleError(message=msg,call=sys.call(),error=e)		
		  })
		if(!.isError(qd)){			
	 		S0 <- structure(
					  c(S0,qd[[1]][which(!(names(qd[[1]]) %in% names(S0)))],
					     "Qnorm"=0.5*sum(qd[[1]]$score^2)),
					 Sigma = attr(qd[[1]],"Sigma"), restarted = attr(S0,"restarted"),				
				   class = "QSResult")			
				
	 	} else { 
			message(qd$message)
			return( structure(S0, tracklist = tracklist, error = qd) )
		}	  
    }	
		
	if(verbose){
	  cat(paste0("Successful minimization by: ",fun.name,
		if(isTRUE(attr(S0,"restarted"))) " [restarted]", " (status = ",S0$convergence,")","\n"))
	}
	if(pl >= 5L){
	  cat("\n")
	  print(S0,pl=pl)
	  cat("\n")
	}
   
	if(length(tracklist) > 0L)
	  attr(S0,"tracklist") <- tracklist
    return(S0)   
}
#' @name multiSearch
#'
#' @title A multistart version of local searches for parameter estimation
#'
#' @description  The function is a multistart version of \code{\link{searchMinimizer}} which selects the best
#' 	root of the quasi-score (if there is any) or a local minimum from all found minima according to the criteria described
#'  in the vignette.
#' 
#' @param x0	 	  numeric, \code{NULL} (default), list,  vector or matrix of starting parameters
#' @param qsd		  object of class \code{\link{QLmodel}}
#' @param ...    	  arguments passed to \code{\link{searchMinimizer}} 
#' @param nstart 	  number of random samples from which to start local searches (if `\code{x0}`=\code{NULL}, then ignored)
#' @param xstart 	  list of starting points for local minimization or root finding by quasi-scoring   
#' @param optInfo 	  logical, \code{FALSE} (default), whether to store original local search results
#' @param multi.start logical, \code{FALSE} (default), whether to perform a multistart local search always otherwise only if first local search did not converge 
#' @param cores		  integer, number of local CPU cores used, default is \code{options(mc.cores,0L)} and if \code{=0} then take all available cores 
#' @param cl 	 	  cluster object, \code{NULL} (default), of class \code{MPIcluster}, \code{SOCKcluster}, \code{cluster}
#' @param pl		  print level, use \code{pl}>0 to print intermediate results
#' @param verbose	  if \code{TRUE} (default), print intermediate output
#' 
#' @details The function performs a number of local searches depending which local method `\code{method}` was passed to
#'  \code{\link{searchMinimizer}}. Either the starting points are given by `\code{x0}` or are generated as an augmented 
#'  design based on the sample set stored in `\code{qsd}`. The function evaluates all found solutions and selects the one which 
#'  is best according to the criteria defined in the vignette. If none of the criteria match, then the parameter for which the lowest value
#'  of the criterion function was found is returned. Multistart searches can be done using a cluster object. Then for each generated/given obervation
#'  a number of \code{cores>1} multistart searches is performed in parallel if \code{fun="mclapply"} using the local cores of each cluster node. 
#' 
#' @return Object of class \code{QSResult} and attribute `\code{roots}`, i.e. the matrix of estimated parameters for which any of
#'  the available minimization methods has been successfully applied. If `code{optInfo}` is \code{TRUE}, then the originally estimtation reuslts
#'  are also returned. The best solution is stored as an attribute named `\code{par}` if found any.
#' 
#' @seealso \code{\link{checkMultRoot}}
#'  
#' @examples 
#'  data(normal)
#'  x0 <- c("mu"=3.5,"sigma"=1.5)
#'  S0 <- multiSearch(x0=x0,qsd,method=c("qscoring","bobyqa"),
#'            opts=list("ftol_stop"=1e-9,"score_tol"=1e-3),nstart=4,
#'             optInfo=TRUE,verbose=TRUE)
#' 
#'  roots <- attr(S0,"roots")
#'  id <- attr(roots,"id")
#'  stopifnot(!is.na(id)) 
#'  id  # index of best root found in matrix roots
#'  attr(roots,"par")  # the final parameter estimate w.r.t. id
#'  
#' @rdname multiSearch
#' @author M. Baaske 
#' @export 
multiSearch <- function(x0 = NULL, qsd, ..., nstart = 10, xstart = NULL, optInfo = FALSE,
		 					multi.start = FALSE, cores = getOption("mc.cores",0L),
								cl = NULL, pl = 0L, verbose = FALSE)
{	 	
	if(!(nstart > 0L))
	 stop("Number of multistart points must be greater than zero!")
 	args <- list(...)
	
	S0 <- 
	 if(!is.null(x0)){
	   if(!is.list(x0))
		 x0 <- .ROW2LIST(x0)
	   # use only first provided method, usually `qscoring`
	   # if non convergence then do a multistart search if enabled
	   # otherwise use a restart with some nloptr minimization routine				
	   if(verbose)
		cat("Starting first local search...","\n")  
	   do.call(searchMinimizer,	c(list(x0=x0[[1]],qsd=qsd,optInfo=optInfo,pl=pl,verbose=(pl>=3L)),args))
	 } else NULL
	
	if(!is.null(S0)){
		if(.isError(S0))
		 message("First local search has errors.")
	    else if(S0$convergence < 0L || S0$convergence == 10) {
		 if(verbose)
		  cat(paste0("First local search did not converge ( status = ",S0$convergence),")\n")						  
		} else {
 		  mroot <- try(.evalRoots(list(S0),opts=args$opts),silent=TRUE)
		  if(.isError(mroot))
			message(.makeMessage("Failed to check first local solution."))					
		  else if(!attr(mroot,"valid")) {
			multi.start <- TRUE							# start multistart search if not a valid root
		}
	  }
	}	
	
    if(is.null(S0) && !multi.start)
		stop("No starting `x0` given. Argument `multi.start` should be set TRUE.")
	RES <- 
	 if(.isError(S0) || multi.start ||
		S0$convergence < 0L || S0$convergence == 10) {	# do not accept convergence by `xtol_rel`
		 if(isTRUE(args$method[1]=="qscoring")) {
			 if(is.null(xstart)){			 
				# generate random LHS starting points for quasi-scoring
				X <- as.matrix(qsd$qldata[seq(attr(qsd$qldata,"xdim"))])
				xstart <- try(multiDimLHS(N=nstart,qsd$lower,qsd$upper,X=X,
								method="augmentLHS",type="list"),silent=TRUE)		 
				if(inherits(xstart,"try-error")) {
					msg <- .makeMessage("Could not generate random starting points in function `multiDimLHS`.")
					if(!is.null(x0)){
						warning(msg)
						return ( structure(S0, message=msg, error=xstart) )
					} else return(.qleError(message=msg,call=match.call(),error=xstart))
				}
				if(verbose)
				 cat(paste("Multistart search using",nstart,"random starting points (LHS)."),"\n")		
			} 
			do.call(doInParallel,
			 c(list(X=xstart,
					FUN=function(x,...){
						searchMinimizer(x,...)						# including a restart by default
					}, cl=cl, cores=cores, fun="mclapply", 
					qsd=qsd, optInfo=optInfo, pl=pl, verbose=(pl>=3L)
			 ), args) )	
		 } else {
			 if(verbose) cat("Starting mulisearch...","\n")
			 if(is.null(xstart)){			 	
			    x0 <-
				 if(is.null(x0))
				  (qsd$lower + qsd$upper)/2.0
				  else if(is.list(x0) || is.vector(x0)) {		
				   unlist(x0)
				  } else { stop("Starting vector 'x0' must be list or a numeric vector.") }
				 if(length(args$control) > 0L) {		
					  opts <- nloptr::nl.opts()
					  opts[names(args$control)] <- args$control
				 } else {
					 opts <- list("algorithm" = "NLOPT_LN_BOBYQA",
						"ftol_abs"=1e-12,"ftol_rel"=1e-7,"xtol_rel"=1.0e-6,"maxeval"=100)
				 }
				 args$control <- list("algorithm"="NLOPT_GN_MLSL", "ftol_rel"=1e-5,
						           "xtol_rel"=1.0e-4, "maxeval"=100, "local_opts" = opts)	  		 
		  		 # no restart
			     args$restart <- FALSE	 
				 # default method then mlsl
			     args$method <- "nloptr"  		 
				 list(do.call(searchMinimizer,c(list(x0=x0,qsd=qsd,optInfo=optInfo,pl=pl,verbose=(pl>=3L)),args)))
			 } else {				 
				 do.call(doInParallel,
					 c(list(X=xstart,										# list of starting points if provided
							 FUN=function(x,...){
								 searchMinimizer(x,...)						# including a restart by default
							 }, cl=cl, cores=cores, fun="mclapply", 
							 qsd=qsd, optInfo=optInfo, pl=pl, verbose=(pl>=3L)
					 ), args) )	
			 }
		 }   	
	 } else { list(S0) }
	# check results
	if(.isError(RES))	
	 return(RES)    
 	# no evaluation for just a single parameter
    if(length(RES) == 1L){
		if(.isError(RES[[1]]) || RES[[1]]$convergence < 0L){
			msg <- .makeMessage("Local search did not converge.")
			message(msg)
			return(.qleError(message=msg,call=match.call(),error=RES))
		}
		return (RES[[1]])
 	}
		
	# check results again
	ok <- which(sapply(RES,function(x) !.isError(x) && x$convergence >= 0L))	
	if(length(ok) == 0L){
		msg <- .makeMessage("All local searches did not converge.")
		message(msg)
		return(.qleError(message=msg,call=match.call(),error=RES))							
	} else if(length(ok) < length(RES)){
		message(paste0("A total of ",length(RES)-length(ok)," local searches did not converge."))							
	}
	
	hasError <- which(!(1:length(RES) %in% ok))
	# get the best roots
	roots <- try(.evalRoots(RES[ok],opts=args$opts),silent=TRUE)
	if(.isError(roots)) {
		msg <- .makeMessage("Could not evaluate best results of local searches")
		message(msg)
		attr(roots,"optRes") <- RES
		attr(roots,"hasError") <- hasError
		return(.qleError(message=msg,call=match.call(),error=roots))	   
	}
	id <- attr(roots,"id")
	if(anyNA(id)){
		msg <- .makeMessage("Could not find any root.")
		message(msg)
		attr(roots,"optRes") <- RES
		attr(roots,"hasError") <- hasError
		return(.qleError(message=msg,call=match.call(),error=roots))
 	}

	structure(RES[[ok[id]]],												# best choice
		"roots"=if(optInfo) roots else NULL,        						# successful optimizations
		"optRes"=if(optInfo) c(list(S0),RES[hasError]) else NULL,			# first result and failed optimization results
		"hasError"=hasError) 	
}


#' @name maxNextSample
#'
#' @title Find new global evaluation point 
#'
#' @description  The function finds the maximum value of the global selection criterion, which basically
#'  is a distance weighted version of the quasi-deviance over the whole parameter domain. 
#' 
#' @param x0		  (named) numeric vector, the starting point
#' @param qsd   	  object of class \code{\link{QLmodel}}
#' @param w 		  density potence
#' @param method	  names of possible minimization routines (see details) 
#' @param control 	  list of control arguments passed to \code{\link[nloptr]{nloptr}} for maximization of selection criterion
#' @param ...		  further arguments passed to \code{\link{covarTx}} and \code{\link{quasiDeviance}}
#' @param obs		  numeric vector of observed statistics, overwrites `\code{qsd$obs}`
#' @param optInfo	  logical, \code{FALSE} (default), not yet used argument (ignored)
#' @param check		  logical, \code{TRUE} (default), whether to check input arguments
#' @param pl		  numeric value (>=0), the print level
#' @param verbose	  if \code{TRUE} (default), print intermediate output 
#' 
#' @return list which contains the new evaluation point and the value of the criterion function
#' @seealso \code{\link{searchMinimizer}} 
#' @export 
maxNextSample <- function(x0, qsd, w = 1, method = c("nloptr","direct"), 
					control = list(), ..., obs = NULL, optInfo = FALSE, check = TRUE,
						pl = 0L, verbose = FALSE) {
	args <- list(...)
	x0 <- if(is.matrix(x0))
				structure(as.numeric(x0),names=colnames(x0))	
			else unlist(x0)
	nms <- names(x0)	
	if(check) {
		.checkArguments(qsd,x0,...)
		stopifnot(is.numeric(pl) && pl >= 0L )	
	}
	# to be sure project starting point
	x0 <- .PROJMED(x0,qsd$lower,qsd$upper)
	# check dimension
	xdim <- attr(qsd$qldata,"xdim")
	if(xdim != length(x0))
		stop("Dimension of starting point 'x0' does not match the problem dimension.")
	# current sampled points (design)
	X <- data.matrix(qsd$qldata[seq(xdim)])
	if(length(control) > 0L) {		
		nlopts <- nloptr::nl.opts()
		nlopts[names(control)] <- control
	} else {
		#nlopts <-
		#	list("algorithm" = "NLOPT_GN_MLSL",
		#		"local_opts" = list("algorithm" = "NLOPT_LN_BOBYQA",
		#				"ftol_abs"=1e-12,"ftol_rel"=1e-7,"xtol_rel" = 1.0e-6,"maxeval" = 1000),
		#		"maxeval" = 100,"ftol_rel"=1.0e-4,"xtol_rel" = 1.0e-4)
		nlopts <- list("algorithm"="NLOPT_GN_CRS2_LM", "xtol_rel" = 1.0e-6,"maxeval" = 1000)
	}	
	ret <- 
		tryCatch({			
			# allocation at C level
			if(.isError(.qdAlloc(qsd,...)))
				stop("Allocation error: cannot request C memory.")			
			fn <-		
		  	 switch(qsd$criterion,
				"qle" = function(x,w,X) { 
					dmin <- .min.distXY(X,rbind(x))										
					qd <- .Call(C_internQD,x)
					cv <- as.numeric(sum(sqrt(qd$sig2)))
					-(dmin*cv*exp(-w*0.5*qd$value))
					#qd <- .Call(C_qDValue,x)
					#cv <- crossValTx(qsd,cvm,theta=x,type="nmsd")
					#rowSums(cv) * exp(-w*0.5*qd)
				},
				"mahal" = function(x,w,X) {
					dmin <- .min.distXY(X,rbind(x))
					#qd <- .Call(C_mahalValue,x)	
					qd <- .Call(C_internMD,x)
					cv <- as.numeric(sum(sqrt(qd$sig2)))
					-(dmin*cv*exp(-w*0.5*qd$value))
				}
			)	
			# maximize selection criterion
			if(verbose)
			 cat("Starting optimization to find next sample point.","\n")
			do.call(nloptr::nloptr,
					list(x0=x0, eval_f=fn, lb=qsd$lower, ub=qsd$upper, opts=nlopts,
						 w=w, X=X))
			
		}, error = function(e) {
			msg <- .makeMessage("Global search by maximizing criterion failed: ", conditionMessage(e))
			message(msg)
			return(.qleError(message = msg,call = sys.call(), error = e))			
		}, finally = { 
			if(!.qdDealloc())
			 stop("Allocation error in C memory management.")
		}
	)
	# check results
	if(.isError(ret) || !is.numeric(ret$solution) || ret$status < 0L) {					
		msg <- .makeMessage("'Nloptr' no convergence by algorithm `",nlopts$algorithm,"`.",
				 if(inherits(ret,"error")) conditionMessage(ret) else "", sep=" ")
		message(msg)
		return(	.qleError(message = msg, call = sys.call(), error = ret, method = "nloptr") )		
	}
	if(verbose)
	 cat(paste0("Found new global evaluation point.","\n"))
	if(!is.null(nms))
	 names(ret$solution) <- nms 
	# quasi-deviance at found solution
	qd <- tryCatch({				
				quasiDeviance(ret$solution,qsd,...,check=FALSE,verbose=verbose)									
			}, error = function(e) {
				msg <- .makeMessage("Error in criterion function: ",conditionMessage(e))
				message(msg)
				.qleError(message=msg,call=sys.call(),error=e)		
			})
	if(.isError(qd))
	 return(structure(qd, "optRes" = if(optInfo) ret else NULL) )
		
	S0 <- structure(c(qd[[1]],
			list("convergence"=ret$status,
				 "message"=ret$message,
				 "iter"=ret$iterations,
				 "start"=x0,
				 "bounds"=which(qd[[1]]$par >= qsd$upper | qd[[1]]$par <= qsd$lower),
				 "Qnorm"=0.5*sum(qd[[1]]$score^2)),
		 		 "method"=method,
				 "criterion"=qsd$criterion),
		 	 Sigma = attr(qd[[1]],"Sigma"), class = "QSResult")
    structure(S0, "optRes"=if(optInfo) ret)	
}

#' @name qle
#'
#' @title Simulated quasi-likelihood parameter estimation
#'
#' @description  This is the main estimation function of the simulated quasi-likelihood estimation approach. 
#' 
#' @param qsd			object of class \code{\link{QLmodel}}
#' @param sim		    user-defined simulation function, see details
#' @param ...			further arguments passed to `\code{sim}` 
#' @param nsim			numeric, number of (initial) simulation replications for each new sample point
#' @param fnsim 		optional, a call returning the number of simulation replications applied to a new
#' 						sample point with the current environment of calling function \code{qle},
#' 						default is the initial value `\code{qsd$nsim}`, respectively `\code{nsim}`
#' @param obs			optional, numeric vector of observed statistics, overwrites `\code{qsd$obs}` if given
#' @param Sigma			optional, constant variance matrix estimate of statistics (see details) 
#' @param qle.opts		options for local search phase
#' @param method		vector of names of local search methods which are applied in consecutive order	
#' @param qscore.opts   list of control arguments passed to \code{\link{qscoring}}
#' @param control		list of control arguments passed to any of the routines defined in `\code{method}` 
#' @param errType		type of prediction variances used for error estimation of the quasi-score, one of "\code{kv}"
#' 					   (kriging prediction variance), "\code{cv}" (cross-validation variances of the statistics),
#' 						"\code{max}" (using the maximum of both) or "\code{cvfull}" (cross-validation of quasi-score) (see details)  
#' @param pl			print level, use \code{pl}>0 to print intermediate results
#' @param verbose  		logical, \code{TRUE} show additional status information
#' @param use.cluster   logical, \code{FALSE} (default), whether to use the cluster environment `\code{cl}` for computations other than model simulations or
#'   a multicore forking which requires to set \code{options(qle.multicore="mclapply")} using at least a number of cpus 
#' 	 cores \code{>1}, e.g. \code{options(mc.cores=2L)}.
#' @param cl			cluster object, \code{NULL} (default), of class "\code{MPIcluster}", "\code{SOCKcluster}", "\code{cluster}" 
#' @param iseed			integer, seed number, \code{NULL} (default) for default seeding of the random number generator (RNG) stream for each worker in the cluster or
#' 						  for parallel processing by "\code{mclapply}", if available on non windows platforms. Note that only using the RNG L'Ecuyer-CMRG
#' 						  leads to reproducible results. Only for \code{iseed} different from \code{NULL} a seed is set including any cluster worker.
#' @param plot 			if \code{TRUE}, plot newly sampled points (for 2D-parameter estimation problems only)
#'
#' @return List of the following objects:
#' 	  \item{par}{ final parameter estimate}
#' 	  \item{value}{ value of criterion function}
#'    \item{ctls}{ a data frame with values of the stopping conditions}
#'    \item{qsd}{ final \code{\link{QLmodel}} object, including all sample points and covariance models}
#' 	  \item{cvm}{ CV fitted covariance models}
#'    \item{why}{ names of matched stopping conditions}
#'	  \item{final}{ final local minimization results of the criterion function, see \code{\link{searchMinimizer}} }
#'	  \item{score}{ quasi-score vector or the gradient of the Mahalanobis distance}
#' 	  \item{convergence}{ logical, whether the QL estimation has converged, see details} 	  
#' 
#'  Attributes: 	 
#'  
#'  \item{tracklist}{ an object (list) of class \code{QDtrack} containing local minimization results,
#'     evaluated sample points and the status of the corresponding iterations}    
#'  \item{optInfo}{ a list of arguments related to the estimation procedure:}
#'   \itemize{
#'    \item{x0:}{ starting parameter vector}
#' 	  \item{W:}{ final weight matrix, equal to quasi-information matrix at \code{theta}, used for both the variance
#' 			 average approximation, if applicable, and as the predicted variance for the (local) sampling of new candidate points
#' 			 according to a multivariate normal distribution with this variance and the current root as its mean parameter.}
#'    \item{theta:}{ the parameter corresponding to \code{W}, typically an approximate root or local minimzer of the criterion function} 
#' 	  \item{status:}{ status of `qle` optimization, whether last iteration was global, local and minimzed}
#' 	  \item{minimum:}{ whether last local minimization was successful}
#' 	  \item{useCV:}{ logical, whether the cross-validation approach was applied}
#' 	  \item{method:}{ name of final search method applied}
#'    \item{nsim:}{ number of simulation replications at each evaluation point}
#' 	  \item{iseed}{ the seed to initialize the RNG}
#'   }
#'     
#' @details
#'  The function sequentially estimates the unknown model parameter. Basically, the user supplies a simulation function `\code{sim}`
#'  which must return a vector of summary statistics (as the outcome of model simulations) and expects a vector of parameters
#'  as its first input argument. Further arguments can be passed to the simulation function by the `\code{\ldots}` argument. The object
#'  `\code{qsd}` aggregates the type of variance matrix approximation, the data frame of simulation runs, the
#'  initial sample points and the covariance models of the involved statistics (see \code{\link{QLmodel}}). In addition, it sets
#'  the criterion function by `\code{qsd$criterion}`, which is either used to monitor the sampling process or minimized itself. The user
#'  also has the option to choose among different types of prediction variances: either "\code{kv}" (kriging variances), "\code{cv}"
#'  (cross-validation-based variances) or "\code{cvfull}" using cross-validation of the quasi-score to estimate the prediction variance
#'  of the quasi-score induced by the kriging models, are available.
#' 
#'  \subsection{Criterion functions}{
#'  The QD as a criterion function follows the quasi-likelihood estimation principle (see vignette)
#'  and seeks a solution to the quasi-score equation. Besides, the Mahalanobis distance (MD) as an alternative 
#'  criterion function has a more direct interpretation. It can be seen as a (weighted or generalized) least squares criterion
#'  depending on the employed type of variance matrix approximation. For this reason, we support several types of variance matrix
#'  approximations. In particular, given `\code{Sigma}` and setting `\code{qsd$var.type}` equal to "\code{const}" treats `\code{Sigma}`
#'  as a constant estimate throughout the whole estimation procedure. Secondly, if `\code{Sigma}` is supplied and used as
#'  an average variance approximation (see \code{\link{covarTx}}), it is considered an initial variance matrix approximation and
#'  recomputed each time an approximate (local) minimizer of the criterion function is found. This is commonly known as an iterative update
#'  strategy of the variance matrix. Opposed to this, setting `\code{qsd$var.type}` equal to "\code{kriging}" corresponds to continuously
#'  updating the variance matrix each time a new criterion function value is required at any point of the parameter space. In this way the
#'  algorithm can also be seen as a simulated version of a least squares approach or even as a special case of the \emph{simulated method of moments}
#'  (see, e.g. [3]). Note that some input combinations concerning the variance approximation types are not applicable since the criterion "\code{qle}",
#'  which uses the QD criterion function, is not applicable to a constant variance matrix approximation.
#'  }
#'       
#'  \subsection{Monte Carlo (MC) hypothesis testing}{
#'  The algorithm sequentially evaluates promising local minimizers of the criterion function during the local phase in order to assess the plausibility
#'  of being an approximate root of the corresponding quasi-score vector. We use essentially the same MC test procedure as in \code{\link{qleTest}}.
#'  First, having found a local minimum of the test statistic, i.e. the criterion function, given the data, new observations are simulated w.r.t. to the
#'  local minimizer and the algorithm re-estimates the approximate roots for each observation independently. If the current minimizer is accepted as an
#'  approximate root at some significance level `\code{local.opts$alpha}`, then the algorithm stays in its local phase and continues sampling around the
#'  current minimizer accoring to its asymptotic variance (measured by the inverse of the predicted quasi-information) and uses the additional simulations
#'  to improve the current kriging approximations. Otherwise we switch to the global phase and do not consider the current minimizer as an approximate root.
#' 
#'  This procedure also allows for a stopping condition derived from the reults of the MC test. We can compare the estimated mean squared error (MSE) with the
#'  predicted error of the approximate root by its relative difference and terminate in case this value drops below a user-defined bound `\code{perr_tol}`
#'  (either as a scalar value or numeric vector of length equal to the dimension of the unknown parameter). A value close to zero suggests a good match of both
#'  error measures. The testing procedure is disabled by default. Set `\code{qle.opts$test=TRUE}` for testing an approximate root during the estimation.
#'  In this case, if for some value of the criterion function its value exceeds `\code{qle.opts$ftol_abs}`, then the corresponding minimizer is tested for an approximate root.
#'  Otherwise the last evaluation point is used as a starting point for next local searches  by a multistart approach when the algorithm is in its global phase.
#'  Note that this approach has the potential to escape regions where the criterion function value is quite low but, however, is not considered trustworthy as an
#'  approximate root. If testing is disabled, then a local minimizer whose criterion function drops below the upper bound `\code{qle.opts$ftol_abs}` is considered
#'  as a root of the quasi-score and the algorithm stays in or switches to its local phase. The same holds, if the quasi-score vector matches the numerical tolerance for 
#'  being zero. A multistart approach for searching a minimizer during the local phase is only enabled if any of the local search methods fail to converge since most of the
#'  time the iteration is started from a point which is already near a root. The re-estimation of parameters during the test procedure also uses
#'  this type of multistart approach and cannot be changed by the user. 
#'  
#'  If one of the other termination criteria is met in conjunction with a neglectable value of the criterion function, we
#'  say that the algorithm successfully terminated and converged to a local minimizer of the criterion function which could be an approximate root of the quasi-score
#'  vector. This can be further investigated by a goodness-of-fit test in order to assess its plausibility (see \code{\link{qleTest}}) and quantify the empirical and
#'  predicted estimation error of the parameter. If we wish to improve the final estimate the algorithm allows for a simple warm start strategy though not yet as an fully
#'  automated procedure. A restart can be based on the final result of the preceeding run. We only need to extract the object
#'  `\code{OPT$qsd}` and pass it as an input argument to function \code{\link{qle}}. 
#'  }
#' 
#'  \subsection{Sampling new points}{
#'  The QL estimation approach dynamically switches from a \emph{local} to a \emph{global search phase} and vise versa for
#'  sampling new promising candidates for evaluation, that is, performing new simulations of the statistical model. Depending on the current value of
#'  the criterion function three different sampling criteria are used to select next evaluation points which aim on potentially improving the quasi-score
#'  or criterion function approximation. If a local minimizer of the criterion function has been accepted as an approximate root, then a local search
#'  tries to improve its accuracy. The next evaluation point is either selected according to a weighted minimum-distance criterion (see [2] and vignette),
#'  for the choice `\code{nextSample}` equal to "\code{score}". In all other cases, for example, if identifiable roots of the QS could not be found
#'  or the (numerical) convergence of the local solvers failed, the global phase of the algorithm is invoked and selects new potential
#'  candidates accross the whole search space based on a weighted selection criterion. This assigns large weights to candidates
#'  with low criterion function values and vise versa. During both phases the cycling between local and global candidates is
#'  controlled by the weights `\code{global.opts$weights}` and `\code{qle.opts.opts$weights}`, respectively. Besides this, the smaller
#'  the weights the more the candidates tend to be globally selected and vice versa during the global phase. Within the local phase,
#'  weights approaching one result in selecting candidates close to the current minimizer of the criterion
#'  function. Weights approximately zero maximize the minimum distance between candidates and previously sampled points and
#'  thus densely sample the search space almost everywhere if the algorithm is allowed to run infinitely. The choice of weights
#'  is somewhat ad hoc but may reflect the users preference on guiding the whole estimation more biased towards either a local
#'  or global search. In addition the local weights can be dynamically adjusted if `\code{useWeights}` is \code{FALSE}
#'  depending on the current progress of estimation. In this case the first weight given by `\code{qle.opts.opts$weights}` is 
#'  initially used for this kind of adjustment. Make sure to export all functions to the cluster environment `\code{cl}` beforehand,
#'  loading required packages on each cluster/compute node, which are used in the model simulation function (see \code{\link{clusterExport}}
#'  and \code{\link{clusterApply}}). 
#'  }
#' 
#'  \subsection{Parallel options}{
#'  Parallel processing of all computations is supported either by the multicore approach (spawning/forking tasks by "\code{mclapply}" for non Windows-based systems) using the parallel
#'  package or a (user-defined) cluster object, e.g. created by \code{\link{makeCluster}}, currently of class \code{"MPIcluster","SOCKcluster","cluster"} which supports calling
#'  the function \code{\link{parLapply}}. By default there is no parallel processing. A cluster object will be automatically created for functions which also could take such object as an argument
#'  (if \code{NULL}) locally on a single host and finalized after function termination.
#'  In this case, the cluster is set up based on a local forking (under Linux) or as a local \code{PSOCK} connection with a number of CPUs available for other operating systems. The option
#'  "\code{mc.cores}", e.g. \code{options(mc.cores=2L}, sets the number of cores used for both types of parallel computations. One can also set an integer value `\code{iseed}` as a seed to
#'  initialize each worker/thread, see also \code{\link{clusterSetRNGStream}}, for reproducible results of any function involving random outputs if \code{mc.cores>1} and RNG kind L'Ecuyer-CMRG. 
#'  Seeding is either done via setting \code{iseed} once the user calls the function \code{\link{qle}} or beforehand using a cluster defined elsewhere
#'  in the user calling program. Then the user should set \code{iseed=NULL} in order to not reinitialize the seed. The seed passed is stored in the attribute `\code{optInfo}$iseed` of the
#'  final return value. If no cluster object is provided, the user sets \code{mc.cores>1L} and \code{options(qle.multicore="mclapply")}, then most computations and model simulations are performed
#'  by function "\code{mclapply}" on a single host. In particular, the user-defined stochastic model simulations can be performed in a cluster environment `\code{cl}` whereas all other
#'  computations are spawned to the cpus of a single host setting \code{use.cluster=FALSE}. We recommend to specify a cluster object once in the user calling program and pass it to functions which take 
#'  an cluster object as their argument, e.g. \code{\link{qle}}. Auxiliary computations, e.g. local optimizations or root findings by \code{\link{multiSearch}}, should be run in parallel on a single host
#'  with many CPUs, e.g. setting \code{options(qle.multicore="mclapply")} and \code{mc.cores>1}. For parallel computations without using a cluster object on a local host simply set both options to "\code{mclapply}"
#'  and \code{mc.cores>1}.  
#'  }
#'  
#'  For a 2D parameter estimation problem the function can visualize the sampling and selection process, which
#'  requires an active 2D graphical device in advance.    
#' 
#'  The following controls `\code{qle.opts}` for the local search phase are available:
#'   \itemize{
#'   \item{\code{lam_max}:}{ upper bound on the maximum eigenvalue of the generalized eigenvalue decomposition of
#' 		the quasi-information matrix and estimated interpolation error (variance) of quasi-score.
#'  	This stops the main iteration sampling new locations following the idea that in this case
#' 		the quasi-score interpolation error has dropped below the estimated precision at best measured by
#' 		quasi-information matrix for `\code{global.opts$NmaxLam}` consecutive iterates.}
#' 	 \item{\code{pmin}:}{ minimum required probability that a new random candidate point is inside the parameter
#'                space. Dropping below this value triggers the global phase. This might indicate
#' 				  that the inverse of the quasi-information matrix does not sufficiently reflect the variance
#' 				  of the current parameter estimate due to a sparse sample or the (hyper)box constraints of the
#' 				  parameter space could be too restrictive.}
#' 	 \item{\code{nsample}:}{ sampling size of candidate locations at the local phase}
#' 	 \item{\code{ntotal}:}{ number of  random samples for evaluation of the integral of total prediction variance of the quasi-score vector at
#'        the local phase with criterion "\code{VQS}"}
#' 	 \item{\code{weights}:}{ vector of weights, \eqn{0\leq\code{weights}\leq 1}, for weighted error sensitivity sampling}
#'    \item{\code{dtol}:}{ scalar value, factor of the mean value of all minimum distances between sample points as a cluster threshold }
#'	 \item{\code{ftol_abs}:}{ upper bound on the function criterion: values smaller trigger the local phase
#'    treating the current minimzer as an approximate root otherwise forces the algorithm to switch to the global phase and vice versa
#'    unless testing is enabled in which case, depending on the result of testing, could even stay in the local phase.}
#'   \item{\code{stopval}:}{ stopping value related to the criterion function value, the main iteration terminates
#' 				     as soon as the criterion function value drops below this value. This might be preferable to a time consuming
#' 					 sampling procedure if one whishes to simply minimize the criterion function or find a first
#' 					 approximation to the unknown model parameter.}
#'   \item{\code{lam_rel}:}{ upper bound on the relative change of `\code{lam_max}`, see `\code{qle.opts}`
#'         The algorithm terminates	its value drops below after a number of `\code{global.opts$NmaxLamRel}` 
#'         consecutive iterations.}
#' 	 \item{\code{xtol_rel}:}{ relative change of found minimizer of the criterion function or root of quasi-score.}
#' 	 \item{\code{ftol_rel}:}{ relative change of the criterion function value.}
#'   \item{\code{maxiter}:}{ maximum allowed global phase iterations }
#' 	 \item{\code{sampleTol}:}{ minimum allowed distance between sampled locations at global phase}	
#' 	 \item{\code{weights}:}{ vector of \eqn{\code{weights}>0} for global sampling}
#'   \item{\code{nsample}:}{ sampling size of candidate locations at the global phase}
#'   \item{\code{NmaxRel}:}{ maximum number of consecutive iterates until stopping according to `\code{xtol_rel}`}
#'   \item{\code{NmaxLamRel}:}{ maximum number of consecutive iterates until stopping according to `\code{lam_rel}`}
#'   \item{\code{NmaxSample}:}{ maximum number of consecutive iterations until stopping according to `\code{sampleTol}`}
#'   \item{\code{NmaxLam}:}{ maximum number of consecutive iterations until stopping for which the generalized eigenvalue of the variance
#' 		 of the quasi-score vector within the kriging approximation model and its total variance measured by the quasi-information matrix
#'       at some estimated parameter drops below the upper bound `\code{qle.opts$lam_max}` }
#'   \item{\code{nstart}:}{ number of random starting points, \code{10*(xdim+1)} (default), if local searches do not converge}
#'  }
#' 
#' 
#' @seealso \code{\link{quasiDeviance}}, \code{\link{qleTest}} 
#' 
#' @examples
#' data(normal)
#'  
#' # main estimation with new evaluations
#' # (simulations of the statistical model)
#' OPT <- qle(qsd,qsd$simfn,nsim=10,qle.opts=list("maxiter"=1)) 
#' 
#' @author M. Baaske
#' @rdname qle 
#' @useDynLib qle, .registration = TRUE, .fixes = "C_"
#' @export  
#' @import parallel stats
#' @importFrom graphics points
qle <- function(qsd, sim, ..., nsim, fnsim = NULL, obs = NULL,
		        Sigma = NULL, qle.opts = list(),
				  method = c("qscoring","bobyqa","direct"), qscore.opts = list(),
				   control = list(), errType = "cv", pl = 0L, verbose = TRUE, 
				    use.cluster = FALSE, cl = NULL, iseed = NULL, plot=FALSE)
{		
	# print information 	
	.printInfo = function(){		
		if(pl > 0L) {
			cat("\n")			
		    cat("Total evaluations...",niter,"\n")			
			cat("Criterion value.....",formatC(ft, digits=6, format="e", big.mark=","),"\n")
			cat("\n")
		}		
	}

	.showConditions = function() {
		if(pl > 0L) {		   
		   cat("Iterations..................",paste0(niter,"\n"))
		   cat("Sampling status.............",paste0(if(status[["global"]]>1L) "global" else paste0(if(status[["global"]]>0L) "(tested)" else "local","\n")))
		   cat("Local method................",paste0(if(!.isError(S0)) {if(any(S0$bounds>0L)) paste0("`",S0$method, "` (success at bounds)") else paste0("`",S0$method,"` (success)") } else "failed"))			
		   if(!.isError(S0) && isTRUE(attr(S0,"restarted"))) cat(" [restarted]","\n") else cat("\n")				
		   cat("Number of replications......",nsim,"\n")			
		   cat("Weight factors..............",paste0("w=",w," [nfail: ",ctls["nfail","val"]," nsucc: ",ctls["nsucc","val"], "]\n"))
	       cat("\n")
			df <- as.data.frame(
					cbind(
					 formatC(signif(as.numeric(xs),digits=6),digits=6,format="e", flag="#"),
					 formatC(signif(as.numeric(xt),digits=6),digits=6,format="e", flag="#"),
					 if(.isError(Snext))
						 formatC(signif(as.numeric(rbind(rep(NA,xdim))[1,,drop=FALSE]),digits=6),digits=6,format="e", flag="#")
   					 else formatC(signif(as.numeric(rbind(Snext$par)[1,,drop=FALSE]),digits=6),digits=6,format="e", flag="#")))
	       		
			dimnames(df) <- list(names(x0),c("Start","Estimate", "Sample"))
			dfv <- as.data.frame(
					cbind(
						formatC(signif(fs,digits=6),digits=6,format="e"),
						formatC(signif(ft,digits=6),digits=6,format="e"),
						if(is.null(Snext) || .isError(Snext)) NA
						else formatC(signif(Snext$value, digits=6),digits=6,format="e")))
			dimnames(dfv) <- list("value",c("","",""))
		
			if(!.isError(S0)){
				df <- cbind(df,formatC(signif(as.numeric(S0$par),digits=6),digits=6,format="e", flag="#"))
				dimnames(df)[[2]] <- c("Start","Estimate", "Sample", "Local")
				dfv <- cbind(dfv,formatC(signif(as.numeric(S0$value),digits=6),digits=6,format="e", flag="#"))
				dimnames(dfv)[[2]] <- c("", "", "", "")
			}					
			print(format(df, digits=6),	print.gap = 2, right=TRUE, quote = FALSE)		    			
			print(format(dfv,digits=6) ,print.gap = 2, right=TRUE, quote = FALSE)
			cat("\n\n")		
		}
		
		if(pl > 1L) {
			# other conditions
			# max of quasi-score depends on whether criterion was minimized (local) or not
			cat("Current stopping conditions: \n\n")			
			cond <- 
			 if(!.isError(S0) && status[["minimum"]]) {
			   	c("|score_max|" = max(abs(S0$score)))				
			 } else c("|score_max|" = NA)
	 
			cond <-
			 c(cond,
			   "lam_max"=unlist(ctls["lam_max","val"]),
			   "lam_rel"=unlist(ctls["lam_rel","val"]),
			   "var_qd"=unlist(ctls["var_qd","val"]),
			   "xtol_rel"=unlist(ctls["xtol_rel","val"]),
			   "ftol_rel"=unlist(ctls["ftol_rel","val"]),
			   "sampleTol"=unlist(ctls["sampleTol","val"]))	
			
		 	print.default(format(cond, digits = 4, justify="left"),	print.gap = 2, quote = FALSE)
		}	
		
		if(pl >= 2L) {
			if(!.isError(S0)){
			 cat("\n\n")
			 print(S0,pl=pl)				 
		    }			
			cat("\n")
		}
		if(pl >= 3) {
			if(!is.null(Stest) && !.isError(Stest)){
				cat("\n\n")
				cat("MC testing: \n\n")
				print(Stest)
			}
			cat("\n")
		}			
		if(pl > 0L) cat("----------------------------------------------------------------------\n\n")	  	
	}	
	
	args <- list(...)
	if(!is.numeric(pl) || pl < 0L)
	  stop("Print level `pl` must be some positive numeric value.")		
  	
  	# simulation increase function `nsim` 
    if(missing(nsim))
	  nsim <- attr(qsd$qldata,"nsim")  	  
    Fnsim <-
     if(is.null(fnsim)){
	  as.call(list(function(n) n, quote(nsim)))
	 } else if(is.call(fnsim)) {		 
		 fnsim[[1]] <- match.fun(fnsim[[1]])		 
		 # make current environment available in function call 
		 environment(fnsim[[1]]) <- environment()		 
		 fnsim
	 } else {
	  stop("Expected numeric value `nsim` or an object of class `call` in argument `fnsim`.")
	 }
	
	# overwrite observed statistics	
	if(!is.null(obs)) {
	  obs <- unlist(obs)
	  if(anyNA(obs) | any(!is.finite(obs)))
		  warning("`NA`, `NaN` or `Inf` values detected in argument `obs`.")
	  if(!is.numeric(obs) || length(obs)!=length(qsd$covT))
		  stop("Object `obs` must be a (named) numeric vector or list of length equal to the number of given statistics in `qsd`.")
	  qsd$obs <- obs
	}
	# verbose if any print level is set
	verbose <- pl>0L  
    # wrapping simulator function
    stopifnot(is.function(sim))
	sim <- match.fun(sim)	
	# silently remove not required	
	.checkfun(sim,args,remove=TRUE)
	simFun <-
	  structure(
		 function(x) {
			 try(do.call(sim,c(list(x),args)))			
		 },
	     class=c("simObj","function")
	  )
		
    # check here 
	var.type <- qsd$var.type    
	.checkArguments(qsd,NULL,Sigma)	
		
	weighting <- (var.type %in% c("wcholMean","wlogMean"))
	# clean or invert Sigma if supplied
	if(!is.null(Sigma)){
		if(var.type == "kriging"){
			Sigma <- NULL	# to be safe
			message("Constant variance matrix `Sigma` is ignored because using the kriging approximation is set.")
		} else if(var.type == "const") {
			# if constant, then invert ones for all iterations
			Sigma <- try(gsiInv(Sigma),silent=TRUE)
			if(inherits(Sigma,"try-error"))
			 stop("Failed to invert initial estimate `Sigma` as a constant variance matrix.")		
		}
	}
	xdim <- attr(qsd$qldata,"xdim")
	# available local optimization method(s) to choose
	nlopt.fun <- c("cobyla","bobyqa","neldermead","direct","directL","lbfgs","nloptr")		   
	all.local <- 
	 if(qsd$criterion == "qle") {		
		c("qscoring",nlopt.fun)
	 } else nlopt.fun	
	# quasi-score options
	qscore.opts <- .qsOpts(qscore.opts,xdim)

	loc <- pmatch(method,all.local)
	if(length(loc) == 0L || anyNA(loc)) {
	 stop(paste0("Invalid local method(s) for criterion `mahal`. Choose one or more of: ", paste(all.local, collapse = ", ")))
	}
	mid <- pmatch("qscoring",method)
	if(!is.na(mid) ) {
		# if using 'qscoring' method always try it first
		if(mid!=1)
		 method <- c("qscoring",method[-mid])		
	}
	# get initial design points
	X <- data.matrix(qsd$qldata[seq(xdim)])
	nstart <- nrow(X)
	xnames <- colnames(X)
	
	# list of consistency checks
	tmplist <- NULL
	tracklist <- structure(list(),class="QDtrack")
			
	## set 'qle.opts'
	local.opts <- .getDefaultoptions(xdim)
	if(length(qle.opts) > 0L) {
		.checkOptions(local.opts,qle.opts)		
		local.opts[names(qle.opts)] <- qle.opts	
	}		
	# setting controls as data frame	
	ctls <- .setControls(local.opts)
				
	# init weights for update location parameter		
	w <- local.opts$weights[1]	
	if(!is.numeric(w) || any(w<0))
	 stop("Sampling weights must be positive.")	
	mWeights <- length(local.opts$weights)
	ntest <- local.opts$ntest
	
	# parallel options for POSIX OS`s only
	tryCatch({
		noCluster <- is.null(cl)				
		if(noCluster) {
			if(.Platform$OS.type != "windows") {
				if(getOption("qle.multicore","mclapply") == "mclapply"){
					if(getOption("mc.cores",0L) == 0L) {
						if(!inherits(np <- try(detectCores(),silent=TRUE),"try-error"))
							options("mc.cores"=np)	
					}						
					# force L'Ecuyer-CMRG RNG
					if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE) || RNGkind()[1] != "L'Ecuyer-CMRG"){
						RNGkind("L'Ecuyer-CMRG")
						if(!is.null(iseed)) 
						 set.seed(iseed)
				 	}													
				} else { message("You have no parallel processing enabled!") }				
			} else {
				message("Please specifiy an external cluster object as an argument to function `qle` for parallel processing!")
			}	
		} else if(any(class(cl) %in% c("MPIcluster","SOCKcluster","cluster"))){
			if(!is.null(iseed)) 
			 parallel::clusterSetRNGStream(cl,iseed)		 				
		} else if(!is.null(iseed)) set.seed(iseed)
		if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)){
			if(!is.null(iseed)) 
			 set.seed(iseed)
			else runif(1)
		}			
	},error = function(e)  {		
	     message(.makeMessage("Could not initialize or manage the cluster environment."))
	}, finally = {
		if(noCluster && !is.null(cl)) {
			if(inherits(try(parallel::stopCluster(cl),silent=TRUE),"try-error"))
				message(.makeMessage("Failed to stop cluster object."))
			cl <- NULL			
			invisible(gc())
		}
	} )	
		   	
 	# select criterion function	
	# W=W,theta=theta,cvm=cvm,
	criterionFun <- function(x,...)
	 quasiDeviance(x,qsd,Sigma,W=W,theta=theta,cvm=cvm,...,check=FALSE, cl=if(isTRUE(use.cluster)) cl else NULL)
 
	.nextCandidate <- function(that=NULL, n=1e3, w=1, QS=NULL, verbose=TRUE) {	
		Y <- try(multiDimLHS(N=n,qsd$lower,qsd$upper,X=X,
					method="augmentLHS",type="matrix"),silent=TRUE)
		#Y <- nextLOCsample(NULL,that,n,qsd$lower,qsd$upper,X,eps,0.05)				     
		if(.isError(Y))
		 stop(message("Uniform sampling failed for posterior density sampling."))	
		if(verbose)
		 cat("Compute cross-validation error sensitivity...\n")
	 
		dists <- .min.distXY(X,Y)
		# average minimum distance from cluster threshold 	 			
		ds <- local.opts$dtol[1]*mean(distC)
		# distance constraint on new sample points
		idx <- which(dists < ds)	
		if(length(idx) > 0L){
		  if((nrow(Y)-length(idx)) < floor(0.1*length(dists)))											
		   message(.makeMessage("Number of sampled candidates too small and less than 10% of sample size:, ", nrow(Y)))				
		  else Y <- Y[-idx,,drop=FALSE]		
		}	 	
		cv <- crossValTx(qsd,cvm,theta=as.matrix(Y),type="nmsd")
		if(.isError(cv) || anyNA(cv))
		 return(.qleError(message="Failed to cross-validation error sensitivities.",call=match.call(), error=cv))		
	 	QD <- criterionFun(as.matrix(Y),value.only=FALSE)	 	
		if(.isError(QD))
		 return(.qleError(message="Failed to compute quasi-deviance.",call=match.call(), error=QD))		
	 	
		qd <-
		  if(!is.null(QS)) {
			 m <- QS$score
			 S <- gsiSolve(QS$varS,use.solve=FALSE)
			 unlist(mclapply(QD, function(q) { as.numeric(crossprod(q$score-m),S%*%as.matrix(q$score-m)) }))
		} else as.numeric(sapply(QD,"[[","value")) 	
				
		id <- try(which.max(rowSums(cv) *exp(-w*qd)),silent=TRUE)		
		if(inherits(id,"try-error") || !is.finite(id) || length(id) == 0L)
			return(.qleError(message="Could not get index of next sample point.",call=match.call(), error=id))
		QD <- criterionFun(rbind(Y[id,]),value.only=FALSE)
		if(.isError(QD))
		 .qleError(message="Failed to compute quasi-deviance at next candidate.",call=match.call(), error=QD)
		structure(QD[[1]], id=id)
	}
	
	# prior sampling function as fallback solution
	.posteriorSamp <- function(that = NULL, n=1e3, nk=1e4, verbose=FALSE) {		
		Y <- nextLOCsample(NULL,that,n,qsd$lower,qsd$upper,X,eps,0.05)				     
		if(.isError(Y))
		 stop(message("Uniform sampling failed for posterior density sampling."))		
		qd <- criterionFun(as.matrix(Y),value.only=TRUE)		
		if(.isError(qd))
		 return(.qleError(message="Failed to compute quasi-deviance.",call=match.call(), error=qd))
		fp <- (max(qd)-qd)
		smp <- try(dsample(fp, Y, n=n, nk=nk),silent=TRUE)
		if(inherits(smp,"try-error")) {
			msg <- message("Could not sample posterior density.")
			return(.qleError(message=msg,call=match.call(),error=smp))			
		}		
		# evaluate at best point (first mode)		
		QD <- criterionFun(unlist(cbind(smp$X)[1,]),value.only=FALSE)
		if(.isError(QD))
		 return(.qleError(message="Failed to compute quasi-deviance value at mode.",call=match.call(), error=QD))
		return( structure(QD[[1]], Y=as.matrix(smp$X)))
	}
	
	## loop init	
	niter <- 1L
	EPS <- .Machine$double.eps
	maxIter <- local.opts$maxiter
	# miniumum allowable distance of sampling candidates to all others
	eps <- 1e-3*prod(abs(qsd$upper-qsd$lower))^(1/xdim)
	# cluster distance threshold (to be computed)
	distC <- 0.0
	## record iteration status
	status <- list("global"=0L, "last"=0L, "minimum"=FALSE, "eval"=FALSE, "reset"=TRUE)			
		  
    # initialize	
	S0 <- Stest <- NULL
	msg <- cvm <- W <- theta <-  NULL
	
	errId <- pmatch(errType,c("kv","cv","cvfull"))
	if(anyNA(errId) || length(errId)!=1L)
		stop("Invalid argument `errType`. Please choose one of `kv`, `cv`, `max`")
	useCV <- errId > 1L
		
	# pre-fit CV models either full or statistics only	
	if(useCV) {
	 if(verbose)
	  cat("Update cross-validation covariance models...\n")
 	 cvm <- try(prefitCV(qsd,TRUE,errType,list(),if(isTRUE(use.cluster)) cl else NULL,verbose),silent=TRUE) 
	  if(.isError(cvm)) {						
		msg <- .makeMessage("Fitting CV models failed. Cannot continue. Please see the error message for details.")
		message(msg)
		return(.qleError(message=msg,call=match.call(), error=cvm)) 
	  }
	}
	# init		 	
	ft <- fs <- fold <- Inf
	Snext <- varS <- qdP <- NULL
	x0 <- (qsd$lower + qsd$upper)/2.0
	xt <- xs <- xold <- x0
		
	ret <- 
	 tryCatch({	 
	   repeat{			   
		   # find new sample point
		   Snext <- 
		    tryCatch({			   											
			   # Start a (multistart) local search	and an additional multistart
			   # search only at local phase in case optimization from `x` failed
			   # Otherwise at global phase always do a multistart search if any of
			   # the methods given in `method` does not converge or has errors										
			   S0 <- multiSearch(xs, qsd=qsd, method=method, opts=qscore.opts,
					   control=control, Sigma=Sigma, W=W, theta=theta, cvm=cvm,
					    check=FALSE, nstart=max(local.opts$nstart,(xdim+1L)*nrow(X)),
						 multi.start=status[["global"]]>1L, cl=if(isTRUE(use.cluster)) cl,
						  pl=pl, verbose=verbose)
			   															
			   # store local minimization results and # last status
			   tmplist <- list("S0"=S0)							   
			   status[["last"]] <- status[["global"]]
			   
			   # Set current iterate to last best mode of quasi-deviance in case of no convergence												
			   if(.isError(S0) || S0$convergence < 0L) {		
				   msg <- .makeMessage("Local minimzation of criterion function did not convergence.")
				   # discretize quasi-deviance and sample modes as starting points
				   if(verbose)
					 cat("Continue by sampling modes of quasi-deviance...\n")
				   S0 <- .posteriorSamp(xs, local.opts$npsample, nk=local.opts$ntotal)
				   if(.isError(S0)) 			  
					 stop(.makeMessage("Error in local minimization and sampling density also failed."))
				   status[["global"]] <- 2L
				   status[["minimum"]] <- FALSE				   				   
			   }
			   else {
				   status[["global"]] <- 0L	
				   status[["minimum"]] <- TRUE
			   }
			   I <- S0$I
			   xt <- S0$par
			   ft <- S0$value					 
			   varS <- S0$varS				
			   # determine the current status: found root or not
			   if(max(abs(S0$score)) > qscore.opts$score_tol || ft > local.opts$ftol_abs) {					    													
					   # check global stopping value (use with care!)
					   if(ft < ctls["stopval","cond"]) { 
						   ctls["stopval","stop"] <- 1L														
					   }
					   # convergence in ftol_rel
					   ctls["ftol_rel","val"] <- abs(ft-fold)/abs(ft)
					   if( ctls["ftol_rel","val"] < ctls["ftol_rel","cond"]) {
						   ctls["ftol_rel","stop"] <- ctls["ftol_rel","stop"] + 1L					
					   } else { ctls["ftol_rel","stop"] <- 0L }						  
			   }											
			   # compute stopping rule maximum generalized eigenvalue				  
			   ctls["lam_max","tmp"] <- ctls["lam_max","val"]
			   ev <- try(geneigen(varS,I,only.values=TRUE),silent=TRUE)
			   if(!inherits(ev,"try-error") && all(is.numeric(ev))){
				   ctls["lam_max","val"] <- max(ev-1.0, na.rm=TRUE)
				   # maximum generalized eigenvalue
				   if( ctls["lam_max","val"] < ctls["lam_max","cond"]) {
					   ctls["lam_max","stop"] <- ctls["lam_max","stop"] + 1L									
				   } else { ctls["lam_max","stop"] <- 0L }
				   # relative change in maximum generalized eigenvalue 				  				  
				   ctls["lam_rel","val"] <- abs(ctls["lam_max","val"]-ctls["lam_max","tmp"])/ctls["lam_max","val"]
				   if(ctls["lam_rel","val"] < ctls["lam_rel","cond"]) {
					   ctls["lam_rel","stop"] <- ctls["lam_rel","stop"] + 1L
				   } else { ctls["lam_rel","stop"] <- 0L }	
			   } else {
				   message("Could not compute generalized eigenvalues so discard this as a stopping condition.")
			   }							  
			 
			   
			   # convergence in xtol_rel
			   ctls["xtol_rel","val"] <- max(abs(xt-xold)/pmax(abs(xt),1.0/qscore.opts$xscale))
			   if( ctls["xtol_rel","val"] < ctls["xtol_rel","cond"]) {
				   ctls["xtol_rel","stop"] <- ctls["xtol_rel","stop"] + 1L					
			   } else { ctls["xtol_rel","stop"] <- 0L }
			   
			   # compute jacknife variance for relative comparisson 
			   if(verbose)
				 cat("compute jacknife variance of quasi-deviance...","\n")
			   qdP <- qDPress(xt,ft,verbose=TRUE)
			   if(!.isError(qdP)){			
				   # storing the sign also, accept only better values as increase of accuracy
				   ctls["var_qd","val"] <- (qdP$var-ctls["var_qd","tmp"])/max(qdP$var,ctls["var_qd","tmp"])				   
			   } else {
				   msg <- .makeMessage("Jacknife variance calculation of quasi-deviance failed.")
				   message(msg)
			   }			   	
			   		   
			   # reached minimum sampling distance ?
			   dm <- attr(.check.distAll(X,xTol=ctls["sampleTol","cond"]),"min")
			   if(dm < ctls["sampleTol","val"])
				   ctls["sampleTol","val"] <- dm
			   if(dm < ctls["sampleTol","cond"]) {
				 	 ctls["sampleTol","stop"] <- ctls["sampleTol","stop"] + 1L				
			   } else { ctls["sampleTol","stop"] <- 0L }					  
			  
			   # maximum iterations/evaluations				
			   # If stopping at global phase, then the current sample point
			   # at maximum weight corresponds to a sampled minimum of the
			   # criterion function if not locally minimized.
			   
			   if(niter >= maxIter){
				   if(status[["global"]] > 0L){					   
					   ctls["maxiter","stop"] <- 1L						  	
				   } else ctls["maxiter","stop"] <- 1L 
			   } 
			   		   
			   # stop before sampling new candidate point
			   if(any(ctls[,"stop"] >= ctls[,"count"])) {
				  # show info
				  .printInfo()		
				  # print stopping conditions				
				  .showConditions()
				  break
			   }		   
			  
			   if(local.opts$useWeights) {										
				   # use user defined weights for cycling if wanted
				   w <- local.opts$weights[(niter %% mWeights)+1L]				   
			   } else {				   
				   # remember successful and failed iteration 
				   if(ctls["var_qd","val"] > 0 && ctls["var_qd","val"] > ctls["var_qd","cond"] ||
					 (ft > local.opts$ftol_abs && ft < (0.9*fold)) ) # some improvement in quasi-deviance if numerically not zero  													 
				   {
					   ctls["nfail","val"] <- 0L
					   ctls["nsucc","val"] <- ctls["nsucc","val"] + 1L												
				   } else if(abs(ctls["var_qd","val"]) > local.opts$tau*ctls["var_qd","cond"]){
					   ctls["nsucc","val"] <- 0L
					   ctls["nfail","val"] <- ctls["nfail","val"] + 1L												
				   } else if( (abs(qdP$var-ctls["var_qd","tmp"])/ctls["var_qd","tmp"]) < local.opts$var_qd_abs){
					   ctls["nfail","val"] <- 0L
				   	   ctls["nfail","val"] <- ctls["nfail","val"] + 1L
					   #ctls["nsucc","val"] <- 0L
					   #ctls["nsucc","val"] <- ctls["nsucc","val"] + 1L					   
				   }
				   if(status$reset){
					   status$reset <- FALSE
					   w <- local.opts$weights[1]					   
				   } else {					   
					   # update weights									
					   if(ctls["nfail","val"] > 0L && 
						  !(ctls["nfail","val"] %% ctls["nfail","cond"])){
						  
						   w <- min(w+log(local.opts$eta[2]*ctls["nfail","val"]),100)
						   #w <- min(local.opts$eta[2]*w,100)
					   } else if(ctls["nsucc","val"] > 0L &&
							   !(ctls["nsucc","val"] %% ctls["nsucc","cond"])){
						  
						   w <- max(w-local.opts$eta[1]*ctls["nsucc","val"]^(-0.5),1)
						   #w <- max(local.opts$eta[1]*w,1)
					   }
					   ctls["var_qd","tmp"] <- qdP$var
			       }
			   }			   
			   # cluster threshold computation			   			   
			   distC <- .min.pair.dist(X)
			   dxt <- .min.distX(X,rbind(xt))
			   
			   # accept roots that are closer to any existing evaluation
			   # point than this is the case for infill sampling (less dense)
			   if(local.opts$doStest && dxt > local.opts$dtol[2]*mean(distC)) {
			    # testing 'ntest' previous minimizers or approximate roots
			    # which signals clustering of roots hence sampling at current root is prefered
			    if(status[["minimum"]]) {
				  # start again with infill sampling after root sampling
				  if(status[["last"]] != 1L && niter > ntest){					  
					  ntrack <- length(tracklist)
					  mroots <- lapply(tracklist[ntrack:(ntrack-ntest)],
							      function(x){
									  if(!.isError(x$S0) &&
									     isTRUE(x$status$minimum && !x$status$eval) &&
										 is.numeric(x$S0$par) ) 
										x$S0$par
									  else NULL  
								  } )
					  okroot <- sapply(mroots,function(x) !is.null(x))
					  if(all(okroot)) {
						  # previous roots must not be too far from current one
						  dxt <- .min.distX(do.call(rbind,mroots),rbind(xt))
						  if(dxt < mean(distC)) {
							  # include current minimzer for testing which is accepted for sure
							  # if it is an approximate (numeric) root of quasi-score but could be
							  # rejected if it is only a local minimizer
							  Stest <- rootCVtest(S0,c(list(xt),mroots),qsd,cvm,opts=qscore.opts,
									    control=control,Sigma=Sigma,W=W,theta=theta,alpha=local.opts$alpha,
									     check=FALSE,info=FALSE,pl=pl,verbose=verbose)				
							  tmplist <- c(tmplist,list("Stest"=Stest))
							  if(.isError(Stest)) {
								  msg <- paste(c("Could not test previous minimizer or root: \n\t ",
										  format(xold, digits=6, justify="right")),collapse = " ")
								  message(msg)
								  Stest <- NULL					  
							  } else {
								  passed <- sapply(Stest,function(x) attr(x$test,"passed"))
								  # all ntest previous roots/minimizers are required to be 
								  # roots given the current approximation QL model 
								  if(all(passed)) {	status[["global"]] <- 1L }
								  else message(.makeMessage("Last CV candidates rejected."))
							  }
					   } else { message(.makeMessage("Previous minima are too far away.")) }
					  } else { message(.makeMessage("Last CV candidate minimum has already been evaluated.")) }					    
				   }				  				
			     } else { message(.makeMessage("Do not test an approximate sampled minimum.")) }
			   } else {
				   # disregard minimum otherwise hence sampling a minimum is postponed
				   # until distance has decreased to an allowable value by infill-sampling first
				   message(.makeMessage("Testing: Root is not feasible for evaluation being too close."))
			   }
			   
			   # check if it is reasonable to sample minimum			   
			   status[["eval"]] <-
				if(status[["global"]] == 1L){
			       # sample at best guess so far,
				   # i.e. the minimum of quasi-deviance or root of quasi-score
				   Stmp <- S0
				   status[["reset"]] <- TRUE
				   TRUE
			    } else {
				   # Infill sampling approach
				   Stmp <- .nextCandidate(NULL,local.opts$nsample,w,S0,verbose=TRUE)
				   
				   # TODO: coninuous approach			   
				   #   Stmp <- maxNextSample(x0, qsd, w, method = c("nloptr","direct"), 
				   #		control=control, cvm=cvm, Sigma=Sigma, W=W, theta=theta,
				   #			check=FALSE,  pl=pl, verbose=verbose)
				   #
				   #   #try to compute evaluation criterion based on a uniform random sample
				   #   if(.isError(Stmp) || Stmp$convergence < 0L) {				   				   
				   #      message(paste0("Trying to discretize sampling criterion for next evaluation point.","\n"))
				   #		Stmp <- .nextCandidate(NULL,local.opts$nsample,w,S0,verbose=TRUE)				   
			   	   #   }
				   FALSE
			   } 
			   # criterion at selected point `Snext`
			   tracklist <- c(tracklist, list(c(tmplist,"Snext"=list(Stmp),"qdP"=qdP,"status"=list(status))))				
			   if(.isError(Stmp))					
			    stop(.makeMessage("Could sample next evaluation point."))		      
			   Stmp
			   									   
			  }, error = function(e) {						
					msg <- .makeMessage("Sampling new candidates failed: ",conditionMessage(e))
					message(msg)
					.qleError(message=msg,call=match.call(),error=e)
			 	 }
		 )   		
		 	 
		 # show info
		 .printInfo()		
		 # print stopping conditions				
		 .showConditions()
		 
		 # stop on error: no next sampling point found
		 if(.isError(Snext)) { break }
		 
		 # next sampling location (optional: plot iterates (2D)) 
		 if(plot && xdim < 3L) {
			 p <- rbind("xs"=Snext$par,xt) 
			 if(xdim == 1L) 
				 p <- cbind(p,c(0,0))					 
			 cols <- if(status[["global"]]>1L) "blue" else "green"
			 try(points(rbind(p[1,,drop=FALSE]),pch=8,cex=1,col=cols,bg=cols))
			 try(points(rbind(p[2,,drop=FALSE]),pch=21,cex=0.5,col="magenta",bg="magenta"))					
		 }	 		 
		 	
		 # run new simulations and update QL/CV models		
		 tryCatch({
			   # next number of simulations: define rule to increase `nsim` 
			   environment(Fnsim) <- environment()
			   nsim <- try(eval(Fnsim),silent=TRUE)
			   if(inherits(nsim,"try-error") || !is.numeric(nsim)){
				   stop("Number of model simulations is not numeric. Stop now.")						 
			   }
			   # simulate at new locations (use number of simulations qsd$nsim by default)
			   newSim <- simQLdata(simFun, nsim=nsim, X=Snext$par, cl=cl, verbose=verbose)
			   if(.isError(newSim)){
				   tracklist[[length(tracklist)]]$newSim <- newSim
				   stop(paste(c("Cannot simulate observations at new candidate point: ",
					 format(Snext$par, digits=6, justify="right")),collapse = " "),"\n\t(",newSim$message,")")								  								
			   }
			   # try to update current QL model
			   qsd <- .updateQLmodel(qsd, Snext$par, newSim, fit=TRUE,
					     cl=if(isTRUE(use.cluster)) cl else NULL, verbose=verbose)
			   
			   # some really severe error stopping main iteration
			   if(.isError(qsd)){
				   stop(paste(c("Cannot update QL model at candidate points.: ",
					  format(Snext$par, digits=6, justify="right")),collapse = " "))
			   }
			   # refit CV models if applicable
			  if(useCV) {
			    if(verbose)
				  cat("Update cross-validation covariance models...\n")
			    cvm <- try(prefitCV(qsd,TRUE,errType,list(),if(isTRUE(use.cluster)) cl else NULL,verbose),silent=TRUE) 
			    if(.isError(cvm))						
				  stop(.makeMessage("Fitting CV models failed. Cannot continue. Please see the error message for details."))				   		
		 	  }
		   }, error = function(e) {
			   stop(.makeMessage("Failed to simulate and update at new sample point: ", conditionMessage(e)))			  	
			  }
		)									        
		# update iterates and next starting points
		Stest <- NULL
		fold <- criterionFun(xold,value.only=TRUE)[[1]]
		if(ft < fold){
		 xs <- xt										
		 fs <- ft		 	
		}	
		if(weighting) {									# update weighting matrix only if
			W <- varS									# local minimizer has been found    						
			theta <- xt			  					
		}
		xold <- xt										 
		fold <- ft	
		niter <- niter+1L				
		# update current sample matrix
		X <- as.matrix(qsd$qldata[seq(xdim)])			
	}  # end main loop			
	NULL			
	}, error = function(e) {
		 msg <- .makeMessage("Current iteration stopped unexpectedly: ", conditionMessage(e))
		.qleError(message=msg,call=match.call(),error=e)							
	}, finally = {
		  if(noCluster && !is.null(cl)) {
			if(inherits(try(parallel::stopCluster(cl),silent=TRUE),"try-error"))
			   message(.makeMessage("Failed to stop cluster object."))
		  	cl <- NULL			
			invisible(gc())
		  }
	} ) # end outer tryCatch	
	
	err <- 
	if(.isError(ret)){
	 msg <- c(msg,ret$message)
 	 message(msg)
	 ret
 	} else NULL
	
	# only for estimte theta=(xt,ft)	
	ctls["stopval",c(1,3)] <- c(ft,ft < ctls["stopval","cond"])	
	ctls["maxiter",c(1,3)] <- c(niter,niter >= maxIter)		
	
	ctls <-
	 if(.isError(S0)) {
	   message("Last search results have errors. Please see list element `\"final\"`.")	   
	   ctls[1:8,-3]
	 } else {		
	   val <- if(!is.null(S0)) max(abs(S0$score)) else Inf
	   ctls <- rbind(ctls[1:8,-3],   										# remove `tmp` column
		    	as.data.frame(cbind("cond" = qscore.opts$score_tol,
			 				        "val" = val,
									"stop" = as.integer(val < qscore.opts$score_tol),
									"count" = 1L),
		   		     row.names = ifelse(qsd$criterion == "qle","score","grad"),
	   			check.names = FALSE))							  	
	}
	# names of active stopping conditions		
	arg.names <- row.names(ctls[which(ctls[,"stop"] >= ctls[,"count"]),])
	
	structure(
	    list("par"=xt,
			 "value"=ft,
			 "ctls"=ctls,			
		 	 "qsd"=qsd,
			 "cvm"=cvm,
			 "why"=arg.names,
			 "final"=S0,			 		
			 "convergence"=(status[["minimum"]] && length(arg.names) > 0L),
			 "message"=msg),	 	
		tracklist = tracklist,		
		optInfo = list("W"=W,
					   "theta"=theta,						   
			    	   "status"=status,
				  	   "minimum"=status[["minimum"]],
				  	   "useCV"=useCV,
				  	   "method"=S0$method,
				  	   "nsim"=nsim,
					   "iseed"=iseed),
		 class = "qle", call = match.call(), error = err)  	 	 
}

#' @name print.qle
#' 
#' @title print results of class \code{qle}
#' 
#' @description S3 method to print the results from \code{\link{qle}}.
#' 
#' @param x      object of class \code{qle} from \code{\link{qle}}
#' @param pl 	 numeric (positive) value, the print level (higher values give more information)
#' @param digits number of digits to display
#' @param format format character(s), see \code{\link{formatC}}
#' @param ... 	 ignored, additional arguments
#' 
#' @rdname print.qle
#' @method print qle
#' @export 
print.qle <- function(x, pl = 1L, digits = 4, format="e",...){	
	if(.isError(x)) 
	  stop("Estimation result has errors.")	
	if(!inherits(x,"qle"))
	  stop("Method is only for objects of class `qle`.")
    if(!is.numeric(pl) || pl < 0L )
	  stop("Print level must be a positive numeric value.")
  
	if(pl >= 0L) {
	  	if(x$qsd$criterion == "mahal")
		 cat("Mahalanobis distance: \n\n",x$value,"\n\n")
		else			
		 cat("Quasi-deviance: \n\n",x$value,"\n\n")
	  	
		cat("Estimate:\n\n")
		print.default(formatC(signif(x$par, digits = digits), digits = digits, format=format, flag="#"),
				print.gap = 4, quote = FALSE)
		
		cat("\n")
		if(x$convergence >= 0L) {
			by <- x$ctls[x$why,"val"]
			names(by) <- x$why
			cat("Convergence: ")
			if("maxeval" %in% x$why)
			  cat("maximum evaluations reached.\n\n")
		    else cat("\n\n")
			print.default(formatC(by, format=format, digits=digits), right=FALSE, print.gap=2,	quote=FALSE)
		} else cat("(None of convergence criteria matched.)\n\n")
		 
	}
  	if(pl >= 2L) {
		cat("\n")	
		nsim <- unlist(x$ctls["maxeval","val"])*attr(x,"optInfo")$nsim	
		cat("Evaluations: ",unlist(x$ctls["maxeval","val"])," (simulations: ",nsim,")\n\n")
		cat("Variance approximation type: ",x$qsd$var.type,"\n")	
	}
	if(pl >= 10L) {
		if(!.isError(x$final)) {
			cat("\n\n ***  Final results *** \n\n\n")			
			print(x$final)
	 	}
		if(x$qsd$var.type != "kriging") {
			W <- attr(x,"optInfo")$W
			if(!is.null(W)) {
				cat("Weighting matrix: \n\n W = \n\n")
				print(W)
				cat("\n")
			}
		}
	}
	invisible(x)
}

#' @name print.QSResult
#' 
#' @title print results of class \code{QSResult}
#'
#' @description S3 method to print the results from \code{\link{qscoring}}.
#' 
#' @param x  		object of type \code{QSResult} from \code{\link{qscoring}}
#' @param pl		numeric positive value, the print level (higher values give more information) 
#' @param digits 	number of digits to display
#' @param format 	format character(s), see \code{\link{formatC}}
#' @param ... 	    ignored, additional arguments
#' 
#' @rdname print.QSResult
#' @method print QSResult
#' @export 
print.QSResult <- function(x, pl = 1L, digits = 4, format="e", ...) {	
	if(.isError(x)) 
	  stop("Quasi-scoring iteration had errors.")
	if(!inherits(x,"QSResult"))
	  stop("Method is only for objects of class 'QSResult'")	
    if(!is.numeric(pl) || pl < 0L )
	 stop("Print level must be a positive numeric value.")
 	cat("\n")
	msg <- unlist(strsplit(paste(x$message, "\n" ),':'))	
	cat(paste0("Local method:\n\n `",x$method,"`",if(isTRUE(attr(x,"restarted"))) "(restarted)"),"\n")
	cat("\n")	
	if(x$criterion == "qle"){
	  cat("Quasi-deviance:\n\n")
	} else cat("Mahalanobis-distance:\n\n")
	cat(formatC(signif(x$value,digits=digits), digits=digits, format=format, big.mark=","),"\n")
	df <- as.data.frame(
		   cbind(
			"Start"=formatC(signif(as.numeric(x$start),digits=digits),digits=digits,format=format, flag="#"),
			"Estimate"=formatC(signif(as.numeric(x$par),digits=digits),digits=digits,format=format, flag="#")			
		   ))	
	if(!is.null(x$score)){		
	  df.score <- cbind(df,as.data.frame(cbind("Quasi-score"=
		  formatC(signif(as.numeric(x$score),digits=digits),digits=digits,format=format, flag="#"))))	  
	}
	cat("\n")
	print(format(df.score, digits=digits), print.gap = 2, right=TRUE, quote = FALSE)
	cat("\n")
	cat("Iterations............",x$iter,"\n")
	cat("Status................",x$convergence,"(",msg[1],")\n\n")	
	cat(msg[2],"\n")	
	cat("\n")	
	if(pl >= 10L) {
		Sigma <- attr(x,"Sigma")		
		if(!is.null(Sigma)) {	
			cat("\nApproximation of variance matrix: \n\n Sigma = \n\n")
			print(formatC(signif(Sigma,digits=digits),digits=digits,format=format, flag="#"),quote = FALSE)
			cat("\n\n")
		}
		if(!is.null(attr(x,"call"))){
			cat("\nCall:\n\n")
			cat(paste(deparse(attr(x,"call")), sep="\n", collapse = "\n"), "\n\n", sep="")
		}
	}	
	invisible(x)
}

# Rejection sampling using `rmvnorm` to draw
# from the multivariate normal distribution and
# reject those samples which do not fall into the domain
.rejectSampling <- function(N, p, x, S, lb, ub, maxTry = 1e+6){
	doN <- N
	ntry <- 0L
	nMax <- 1e+8
	Y <- matrix(double(0L),ncol=length(x))
	colnames(Y) <- names(x)
	
	while(N > 0L){			  
		nNew <- ifelse (N/p > nMax, N, ceiling(max(N/p,100)))
		X <- mvtnorm::rmvnorm(nNew, mean=x, sigma=S)		
		id <- apply(X,1,function(x) all(x>lb & x<ub))
		naccept <- sum(id)		
		if (naccept == 0L) {		  
		  if(ntry > maxTry) {
		  	warning(paste0("Maximum iterations reached in rejection sampling. Could not sample ",N," points."))
			break;
		  }
		  ntry <- ntry + 1L
		  next				
	 	}
		naccept <- min(naccept, N)		
		Y <- rbind(Y,X[which(id)[1:naccept],,drop = FALSE])		
		N <- N - naccept 
	}
	return( structure(Y, success = NROW(Y) > 0.9*doN, ntry = ntry) )
}

#' @name nextLOCsample
#'
#' @title Generate a random sample of points
#'
#' @description Generate a random sample of points as a set of candidate points for evaluation 
#'
#' @param S			    variance matrix, \code{NULL} (default), of sample points 
#' @param x			    numeric vector, \code{NULL} (default), multivariate mean vector of the normal distribution from
#' 						which the random sample is going to be generated if \code{S} is given#' 						
#' @param n				number of points to sample
#' @param lb			vector of lower bounds of the (hyper)box
#' @param ub			vector of upper bounds of the (hyper)box
#' @param X				matrix of surrogate design points, default \code{NULL}, no minimum distance checks used
#' @param eps			minimum required distance of two neighboring locations
#' @param pmin			minimum required ratio of points falling inside the hyperbox of parameters
#' @param invert	    optional, \code{invert=FALSE} (default) for no inversion of `\code{S}`
#'
#' @return Matrix of sampled locations.
#'
#' @details If `\code{S}` is given, then the function generates a random sample of points with mean and
#'   variance given by `\code{x}` and `\code{S}`, respectively, according to a (truncated) multivariate
#'   normal distribution (using rejection sampling) to match the parameter space given by the lower and
#'   upper bound vectors. Otherwise random points are generated according to a uniform distribution on
#'   the parameter space. 	 
#'
#' @examples
#'  X <- nextLOCsample(matrix(c(1,0,0,1),nr=2), c(0,0), 10, c(-0.5,-0.5), c(0.5,0.5))
#'  
#' @author M. Baaske
#' @rdname nextLOCsample
#' 
#' @importFrom mvtnorm pmvnorm rmvnorm
#' @export
nextLOCsample <- function(S = NULL, x = NULL, n, lb, ub, X = NULL, eps = 1e-5, pmin = 0.05, invert = FALSE) {
	if(is.null(S)) {
		Y <- try(sapply(seq_len(length(x)),function(i) runif(n,lb[i],ub[i])),silent=TRUE)
		if(.isError(Y) || !is.matrix(Y)) {
			msg <- .makeMessage("Uniform sampling failed.")
			message(msg)
			return(.qleError(message=msg,call=match.call(),error=Y))
		}
		p <- 1.0
		colnames(Y) <- names(x)	
	} else {
		if(!is.numeric(S) || anyNA(S) || !is.matrix(S))
			stop("Variance matrix has `NaN` values. A matrix?")		
		if(invert) 
			S <- try(gsiInv(S),silent=TRUE)
		if(.isError(S)) {
			msg <- .makeMessage("Failed to invert matrix.")
			message(msg)
			return(.qleError(message=msg,call=match.call(),error=S))		
		}
		if(anyNA(S) || !is.matrix(S))
			warning("Variance matrix has `NaN` values. A matrix?")
		if( rcond(S) < sqrt(.Machine$double.eps)) {
			warning(" Variance matrix is ill-conditioned.")
			if( (is.pos = .isPosDef(S)) != 0L )
				return(.qleError(message=.makeMessage("Variance matrix not positive definite: ",is.pos),
								call=match.call(), S=structure(S,is.pos=is.pos))) 
		}
		
		# get acceptance ratio from multivariate normal distribution
		p <- mvtnorm::pmvnorm(lower=lb, upper=ub, mean=x, sigma=S)
		if (p < pmin) {
			msg <- .makeMessage(sprintf("Sampling local candidates is too inefficient (coverage prob = %.6f of bounding box).", p))
			message(msg)
			return(.qleError(message=msg,call=match.call(),error=p) )
		}
		# sampling
		Y <- try(.rejectSampling(n, p, x, S, lb, ub),silent=TRUE)
		if(inherits(Y,"try-error") || !attr(Y,"success")) {
			msg <- .makeMessage("Rejection sampling failed.")
			message(msg)
			return(.qleError(message=msg,call=match.call(),error=Y))
		}						
	} 
	if(!is.null(X)) {
		# get min distances
		dists <- .min.distXY(X,Y)
		# remove too close sample candidates							
		idx <- which(dists < eps)					
		if(length(idx) > 0L){
			Y <- Y[-idx,,drop=FALSE]
			if(nrow(Y) < floor(0.1*length(dists))) {											
				msg <- .makeMessage("Number of sampled candidates too smallis less than 10% of sample size:, ", nrow(Y))
				message(msg)
				return(.qleError(message=msg,call=match.call(),error=Y))
			}
			attr(Y,"idx") <- idx
		}			
	}
	return ( structure(Y, p = p) )
}

#' @name qscoring
#' 
#' @title Quasi-scoring iteration
#'
#' @description The function numerically solves the quasi-score equation by a root finding algorithm similar to Fisher's scoring method.
#'
#' @param qsd    	object of class \code{\link{QLmodel}}
#' @param x0		(named) numeric vector, the starting parameter
#' @param opts		quasi-scoring options, see details
#' @param Sigma	    a prespecified constant variance matrix estimate
#' @param ...	    further arguments passed to the function \code{\link{covarTx}}
#' @param check		logical, \code{TRUE} (default), whether to check input arguments
#' @param cvm		optional, either list of covariance models of the statistics for cross-validation based estimation of prediction variances of the statistics or
#' 					of class \code{cv} or list of cross-validation models of class \code{cvfull} of the QL model \code{qsd} for computation of the quasi-deviance
#' 					and error estimation of the quasi-score approximation w.r.t the kriging prediction models (see \code{\link{prefitCV}})  	
#' @param Iobs		optional, if \code{TRUE} (default), compute observed quasi-information by finite difference approximations of the quasi-score vector
#' @param verbose   logical, \code{FALSE} (default), otherwise print intermediate output
#'
#' @return List of results of quasi-scoring iteration with elements:
#'  \item{convergence}{ integer, why scoring iterations stopped}
#'  \item{message}{ string, corrsponding to `\code{convergence}`}
#'  \item{iter}{ number of iterations}
#'  \item{value}{ quasi-deviance value}
#'  \item{par}{ solution vector}
#'  \item{score}{ quasi-score vector}
#'  \item{I}{ quasi-information matrix}
#'  \item{start}{ starting point}
#'  \item{method}{ simply: "\code{qscoring}"}
#'  \item{criterion}{ equal to "\code{qle}"} 
#'
#' @details The function implements a step-length controlled quasi-scoring iteration with simple bound
#'  constraints (see also [1,3]) specifically tailored for quasi-likelihood parameter estimation. Due to the typical
#'  nonconvex nature of the (unknown and not further specified) quasi-likelihood function as an objective
#'  function one needs some kind of globalization strategy in order to stabilize the descent step and to avoid a premature
#'  termination. Therfore, we use the quasi-deviance function as a monitor function (see vignette) though it does not
#'  inherit all of the appreciable properties of a true objective function such as among others, for example,
#'  identifying appropriate descent directions. However, these are general numerical obsticles in case of pure root
#'  finding algorithms and need to be addressed elsewhere. 
#'  
#'  \subsection{Quasi-scoring under uncertainty}{ 
#'  The quasi-scoring iteration includes both kinds of prediction variances, kriging-based and those derived from a cross-validation (CV) approach,
#'  which account for the additinal uncertainty induced by the quasi-score approximation model. By default kriging variances
#'  are included in the computation during all iterations through the inverse quasi-information matrix as part of the search direction. If fitted covariance models `\code{cvm}` are supplied by the user
#'  in advance (see \code{\link{prefitCV}}), the variances of prediction errors of each statistic are separately evaluated by the proposed CV
#'  approach for each new point. For the price of relatively high computational costs those prediction variances
#'  are intended to increase the robustness against false roots due to simulation and approximation errors of the quasi-score function.#'   
#'  }
#' 
#'  The following algorithmic options, which can be set by `\code{opts}`, are available:
#'  \itemize{
#'   	\item{\code{ftol_stop}:}{ minimum value of the quasi-deviance for stopping the scoring iteration}
#' 		\item{\code{ftol_abs}:}{ minimum value of the quasi-deviance which is used as a reference value for a local minimizer or if the minimum steplength is reached}
#' 		\item{\code{xtol_rel}:}{ minimum relative stepsize between consecutive estimates }
#' 		\item{\code{step_tol}:}{ tolerance of two scaled consecutive steps }
#' 		\item{\code{grad_tol}:}{ upper bound on the quasi-score vector components,
#' 				 testing for a local minimum of the quasi-deviance in case of a line search failure}
#' 		\item{\code{score_tol}:}{ upper bound on the quasi-score vector components, testing for an approximate root}
#'      \item{\code{maxiter}:}{ maximum allowed number of iterations}
#' 		\item{\code{xscale}:}{ numeric, default is vector of 1, typical magnitudes of vector components of `\code{x0}`, e.g. the order of upper bounds of the parameter space}
#'      \item{\code{fscale}:}{ numeric, default is vector of 1, typical magnitudes of quasi-score components}
#' } 
#'
#' @examples 
#' data(normal)
#' QS <- qscoring(qsd,x0=c("mu"=3.5,"sigma"=0.5),
#'          opts=list("score_tol"=1e-4))
#' 
#' @seealso \code{\link{prefitCV}}
#' 
#' @author M. Baaske
#' @rdname qscoring
#' @export
qscoring <- function(qsd, x0, opts = list(), Sigma = NULL, ..., check = TRUE,
				cvm = NULL, Iobs = TRUE, verbose = FALSE)
{		
	if(check)
	 .checkArguments(qsd,x0,Sigma)
  	if(qsd$criterion != "qle")
	  stop("Quasi-scoring is only valid for criterion `qle`.")
    # just project to be sure
    xdim <- attr(qsd$qldata,"xdim")
	x0 <- .PROJMED(x0,qsd$lower,qsd$upper)    
	Sigma <- varMatrix(qsd,Sigma,...,cvm=cvm)
	if(.isError(Sigma)){
		msg <- paste0("Could not compute variance matrix for CV models.")
		message(msg)
		return(.qleError(message=msg,error=Sigma))
	} 
	qlopts <- list("varType"=qsd$var.type, "useCV"=!is.null(cvm)  && class(cvm)=="cv", "Iobs"=Iobs)
	args.qsd <- list("qlm"=qsd,	"VTX"=Sigma, "X"=as.matrix(qsd$qldata[seq(xdim)]))
	
	# LOO CV models patch
	if(!is.null(cvm) && class(cvm) == "cvfull") {
		cvm <- doInParallel(cvm,
				function(qlm, ...) {				
					Sigma <- varMatrix(qlm,...)
					if(.isError(Sigma))						  
					 return(.qleError(message=.makeMessage("Could not compute variance matrix for CV models."),error=Sigma))					  
					list("qlm"=qlm,"VTX"=Sigma,"X"=as.matrix(qlm$qldata[seq(attr(qlm$qldata,"xdim"))]))				
				}, Sigma=Sigma,...) 
		hasErr <- which(sapply(cvm,function(x) .isError(x)))
		if(length(hasErr)>0) {
			msg <- .makeMessage("Failed to compute variance matrix for cross-validation models.")
			message(msg)
			return(.qleError(message=msg,error=cvm))
		}		  
	}
	# set quasi-scoring options	
	opts <- .qsOpts(opts,xdim)		
	# call quasi-scoring
	structure(try(.Call(C_QSopt, x0, args.qsd, qlopts, cvm, opts), silent=TRUE),call=match.call())	
}

#' @rdname dsample
#' @title Random Samples Generation Through The Wang-Lee and Fu-Wang Algorithms
#' @description \code{dsample} generates a sample of specified size \code{n} from the target density funciton (up to a normalizing constant) based on the Wang-Lee algorithm
#' 
#' @param y 	  numeric values of the density to discretize, must be positive
#' @param X 	  must be a \code{data.frame}.  See \sQuote{Details}.
#' @param nk      positive integer, the number of contours.  See \sQuote{Details}.
#' @param n       non-negative integer, the desired sample size. 
#' @param wconst  scalar numeric between \code{0} and \code{1}.  See \sQuote{Details}.
#' @param info	  logical, default \code{FALSE} additional information about the sampling process
#' 
#' @details  \code{X} has the number of rows equals to the number of discrete base points. In each row, the first element contians the funcitonal value of the target density and the rest elements are the coordinates at which the density is evaluated.  
#' \code{wconst} is a constant for adjusting the volumn of the last contour.
#' 
#' @return The function gives a random sample as a \code{data.frame} with number of rows equals the specified size \code{n} and number of columns equals \code{ncol(x)-1}.
#' 
#' @references
#' Wang, L. and Lee, C.H. (2014). Discretization-based direct random sample generation. Computational Statistics and Data Analysis, 71, 1001-1010. 
#' Lee, C.H. (2009). Efficient Monte Carlo Random Sample Generation through Discretization, MSc thesis, Department of Satistics, University of Manitoba, Canada
#' Wang, L. and Fu, J. (2007). A practical sampling approach for a bayesian mixture model with unknown number of components. Statistical Papers, 48(4):631-653.
#' Fu, J. C. and Wang, L. (2002). A random-discretization based Monte Carlo sampling method and its application. Methodology and Computing in Applied Probability, 4, 5-25.
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca}, Liqun Wang \email{liqun.wang@@umanitoba.ca}
#' @keywords sampling, discretization
#' 
#' @examples 
#' ## The following example is taken from West (1993, page 414).
#' ## West, M. (1993). Approximating posterior distributions by mixture.
#' ##   Journal of the Royal Statistical Society - B, 55, 409-422.
#' 
#' expr <- expression((x1*(1-x2))^5 * (x2*(1-x1))^3 * (1-x1*(1-x2)-x2*(1-x1))^37)
#' sets <- list(x1=runif(1e5), x2=runif(1e5))
#' smp <- dsample(expr=expr, rpmat=sets, nk=1e4, n=1e3)
#' 
#' ##
#' ## More accurate results can be achieved by increasing the number 
#' ## of dicretization points and the number of contours.  
#' @export 
dsample <- function(y, X, n=1e3, nk=1e4, wconst=1, info=FALSE){
		
	stopifnot(!missing(X), nrow(X) > nk)
	stopifnot(is.numeric(nk), is.numeric(n), is.numeric(wconst))	
		
	X <- as.data.frame(X)
	yX <- cbind(y, X)
	yX <- yX[which(y>0),] # data.frame
	yX <- yX[order(yX$y, decreasing=TRUE), ]
	
 	cnt <- graphics::hist(yX$y, breaks=seq(from=min(yX$y), to=max(yX$y), length.out=nk+1), plot=FALSE)
	cnames <- paste("e", seq_len(nk), sep="") # contour name (in order)
	yX$cid <- cut(yX$y, cnt$breaks, include.lowest=TRUE)
	yX$cnt.name <- cnames[match(yX$cid, names(table(yX$cid)))]
	cnt.counts <- cnt$counts
	cnt.mids <- cnt$mids
	names(cnt.mids) <- names(cnt.counts) <- cnames
	
	gpdf <- rev(cnt.counts * cnt.mids)
	gpdf[nk] <- wconst*gpdf[nk]
	rev.cntc <- rev(cnt.counts)
	
	cdf <- cumsum(gpdf)/sum(gpdf)
	cumcdf <- cumsum(cdf)
	
	# counting how many samples are needed from each contour
	pptns <- graphics::hist(stats::runif(n), breaks=c(0,cdf), plot=FALSE)$counts
	names(pptns) <- rev(paste("e", seq_len(nk), sep=""))
	
	# sampling from each contour 
	scnt <- mapply(FUN=sample, MoreArgs=list(replace=TRUE), rev.cntc, pptns)
	idx <- unlist( mapply("+", as.list( c(0, cumsum(rev.cntc))[-(nk+1)] ), scnt) )
	yX <- yX[order(yX$y, decreasing=TRUE)[idx], ]
	robj <- list(yX=yX, X=yX[names(X)], cnt.counts=cnt.counts,
			cnt.mids=cnt.mids, gpdf=gpdf, cdf=cdf, cumcdf=cumcdf, pptns=pptns, scnt=scnt, idx=idx)
	class(robj) <- "dsample"
	return(robj)
}


#' @rdname summary.dsample
#' @title Generating Basic Summary Statistics of Marginal Distributions
#' @description  Producing basic summary statistics (the mean, the standard deviation and the first five modes) from the sample drawn via either the Fu-Wang algorithm or the Wang-Lee algorithm, for all marginal distributions of the target distribution.
#' @param object a \code{data.frame}, contains the sample drawn via either the Fu-Wang algorithm or the Wang-Lee algorithm 
#' @param n the first n samples
#' @param ... more arguments
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca}, Liqun Wang \email{liqun.wang@@umanitoba.ca}
#' @export 
summary.dsample <- function(object, n=5, ...) {
	stopifnot(inherits(object, "dsample"))
	stopifnot(is.numeric(n) && n>0)
	X <- object$X
	means <- colMeans(X)
	modes <- cbind(X)[1:n,]
	stdevs <- do.call(c, lapply(X, stats::sd, na.rm=TRUE))		
	
	robj <- list(means=means, stdevs=stdevs, modes=modes[!duplicated(modes),])
	class(robj) <- "dsample"
	return(robj)
}


