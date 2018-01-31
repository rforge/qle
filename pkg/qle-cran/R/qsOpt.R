# Copyright (C) 2017 Markus Baaske. All Rights Reserved.
# This code is published under the GPL (>=3).
#
# File: 	qsOpt.R
# Date:  	20/10/2017
# Author: 	Markus Baaske
#
# Functions for Quasi-likelihood based estimation as well as
# the quasi-scoring as a local root finding method 

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

# internal function to check 
# the arguments `args` with
# function `fun` 
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

.qsOpts <- function(options = list(), xdim = 1L, pl = 0L) {
	opts <- .addQscoreOptions(xdim)
	opts$pl <- pl
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
	list( "ftol_stop" = 1e-10,								# also used to select best roots
		  "xtol_rel"  = 1e-7,
		  "grad_tol"  = 1e-5,
		  "ftol_rel"  = 1e-8,
		  "ftol_abs"  = 1e-6,								# only for local minima if grad_tol reached as a more restrictive check
		  "score_tol" = 1e-4,								# also used to select best roots
		  "slope_tol" = 1e-4,
		  "maxiter"   = 100,
		  "xscale" = rep(1,xdim),							# scaling independent variables, e.i. parameter theta
		  "fscale" = rep(1,xdim),							# and function values, i.e. QS components 
		  "pl" = 0L)
}

.getDefaultGLoptions <- function(xdim) {
	list("stopval" = .Machine$double.eps,			 		# global stopping value
		 "C_max"   = 1e-3,
		 "xtol_rel" = .Machine$double.eps^0.25,
		 "maxiter" = 100,									# max number of global iterations
		 "maxeval" = 100,									# max number of global and local iterations
		 "sampleTol" = .Machine$double.eps^0.25,			# minimum (euclidean) distance between samples		 
	 	 "weights"=c(50,25,10,5,2,1),		 
		 "nsample" = (xdim+1)*2000,							# number of global random samples
		 "NmaxRel" = 5,		 
		 "NmaxCV" = 3,		 
		 "NmaxSample" = 3,
		 "NmaxLam" = 3,
		 "NmaxQI" = 3,		 		 
		 "Nmaxftol"= 3,
		 "nstart" = 25)										# number of starting points for multistart version at global phase
}

.getDefaultLOCoptions <- function(xdim) {
	list("ftol_rel" = .Machine$double.eps^(1/3),
		 "ftol_abs"	= .Machine$double.eps^0.5,			   # whether local minimizer is numerically zero
		 "lam_max" = 1e-2,								   # less restrictive
		 "pmin" = 0.05,									   # minimum accepted probability of coverage of sample points within search domain
		 "weights" = c(0.005,0.1,0.2,0.4,0.6,0.8,0.995),   # only for sampling with criterion `score`
		 "nsample" = (xdim+1)*1000,						   # number of local random samples
		 "perr_tol" = rep(0.05,xdim),					   # empirical error is 5% smaller than predicted error by inverse QI
		 "nobs"=100,									   # sampling size (root testing)
		 "alpha" = 0.05,							       # significance level testing a root		 
		 "eta" = c(0.025,0.05),							   # c("decrease"=0.05,"increase"=0.075) additive step size	
		 "nfail" = 3,									   # number of failed (not yet improved) iterations until next decrease of weights 
		 "nsucc" = 3,									   # number of successful iterations until next increase of weights 
		 "nextSample" = "score",						   # default selection criterion
		 "useWeights" = TRUE,							   # do not dynamically adjust weights and cycle through the weights								   
		 "test" = FALSE)								   # do not test approximate root		 
}

.setControls <- function(globals,locals) {
	defaults <- c("C_max","lam_max","xtol_rel","stopval","sampleTol",
				  "nfail","nsucc","ftol_rel","maxiter","maxeval")
	optlist <- c(globals,locals,"score_tol")
	namc <- match.arg(names(optlist), choices = defaults, several.ok = TRUE)
	ctls <- data.frame(cbind("cond" = unlist(optlist[namc]),
						     "val" = 0, "tmp"=0, "stop" = 0,
							 "count" = c(1,globals$NmaxCV,globals$NmaxRel,1,1,
									 	 globals$NmaxSample,globals$Nmaxftol,globals$NmaxLam,1,1)),
			row.names = namc, check.names = FALSE)
	
	# init some controls	
	ctls["sampleTol","val"] <- 1E100
	ctls[c("C_max","lam_max"),"val"] <- rep(1,2)		
	return (ctls)
}

#' @name getDefaultOptions
#' 
#' @title Print default options for optimization
#' 
#' @description Print default options for global and local optimization in function \code{\link{qle}}
#' 
#' @param xdim 		dimension of the unknown model parameter
#' 
#' @return List of options.
#' 
#' @details The function returns a lists of available options
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
		 "qle_local_opts" = .getDefaultLOCoptions(xdim),
		 "qle_global_opts" = .getDefaultGLoptions(xdim))
}


## Internal, not exported
## 'points' is best chosen as matrix
cverrorTx <- function(points, Xs, dataT, cvm, Y, type, cl = NULL) {	
	# extract full predictions	
	dfx <- as.data.frame(extract(Y,type="mean"))
	useMax <- (attr(cvm,"type") == "max")
	dfs2 <-
	 if(useMax) {
	   as.data.frame(extract(Y,type="sigma2"))			
	 } else NULL
	# number of fitted cov models equals 
	# number of blocks for jackknife variance
	n <- length(cvm)
	np <- nrow(dfx)			

    # prediction function for CV
	statsCV <- function(covT,points,Xs,dataT) {				
		 id <- attr(covT,"id")
		.COL2LIST(predictKM(covT,points,Xs[-id,],dataT[-id,]))	   			
	} 	
	
	L <- tryCatch(
		  doInParallel(cvm, statsCV, points=points, Xs=Xs, dataT=dataT, cl=cl)
			,error = function(e) {
				msg <- .makeMessage("Cross-validation prediction failed: ",
						conditionMessage(e))
				message(msg)
				.qleError(message=msg,call=match.call(),error=e)
			}
	)
	# on error
	if(.isError(L))
	  return(L)
  
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
								 mod <- cvm[[k]]
								 id <- attr(mod,"id")
								 varKM(mod[[i]],points[k,], Xs[-id,],dataT[-id,i])								 
							 }
				 	)	 		
					(dfx[,i]-diag(y))^2/sigK 		  	 
			 },			 
			 "msd" = { rowMeans((y - dfx[,i])^2) },				# CV based mean squared deviation (prediction uncertainty)
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
				 if(n!=np)										    # at sample points (model validity)
				   stop("Standardized cross-validation MSE can be computed only if number of samples equals number of CV models.")			     	 
				   sigK <- sapply(1:length(cvm),
							 function(k) {
								mod <- cvm[[k]]
								id <- attr(mod,"id")
								varKM(mod[[i]],points[k,], Xs[-id,],dataT[-id,i])										
		         			 })	 
			     mean((dfx[,i]-diag(y))^2/sigK) 			  	 
			 },
			 "sigK" = { 										# standardized (by kriging variance) cross-validation error				 
				if(n!=np)
				 stop("Leave-one-out kriging variance is only available if the number of sample points equals the number of CV models.")				   	
				sapply(1:length(cvm),
						 function(k) {
							 mod <- cvm[[k]]
							 id <- attr(mod,"id")
							 varKM(mod[[i]],points[k,], Xs[-id,],dataT[-id,i])								 
						 }
				 )				 		  	 
			 }			 
		)		
	}))	
}

#' @name crossValTx 
#'
#' @title Prediction variances by cross-validation
#'
#' @description The function estimates the prediction variances by a cross-validation approach (see vignette) applied to
#'  each sample means of the involved statistics. 
#'
#' @param qsd   	object of class \code{\link{QLmodel}}
#' @param cvm		list of prefitted covariance models from function \code{\link{prefitCV}}
#' @param theta		optional, default \code{NULL}, list or matrix of points where to estimate prediction variances
#' @param type		name of prediction variance measure 
#' @param cl	    cluster object, \code{NULL} (default), of class "\code{MPIcluster}", "\code{SOCKcluster}", "\code{cluster}"
#' 						 
#' 	
#' @return A matrix of estimated prediction variances for each point given by the argument \code{theta} (rows)
#'  and for each statistic (columns).  
#'
#' @details	Other than the kriging prediction variance, which solely depends on interdistances of sample points
#'  and estimated covariance parameters of some assumed to be known spatial covariance structure, the cross-validation
#'  based approach (see [4] and the vignette) even takes into account the predicted values at `\code{theta}` and thus can be seen as a more robust
#'  measure of variability between different spatial locations. By default, `\code{theta}` equals the current sampling set 
#'  stored in the object `\code{qsd}`.
#' 
#'  If we set the error `\code{type}` equal to "\code{cve}", the impact on the level of accuracy (predicting at unsampled
#'  points) is measured by the \emph{delete-k jackknifed variance} of prediction errors. This approach does not require further
#'  simulations as a measure of uncertainty for predicting the sample means of statistics at new candidate points accross the parameter space.
#'  Note that if the attribute \code{attr(cvm,"type")} equals "\code{max}", then the maximum of kriging and CV-based prediction
#'  variances is returned. 
#' 
#'  In addition, other measures of prediction uncertainty are available, such as the \emph{root mean square deviation}
#'  (\code{rmsd}) and \emph{mean square deviation} (\code{msd}) or the \emph{standardized cross-validation error}
#'  (\code{scve}). The details are explained in the vignette. In order to assess the predictive quality of possibly
#'  different covariance structures (also depending on the initial sample size), including the comparison of different
#'  sizes of initial sampling designs, the following measures [8] are
#'  also available for covariance model validation and adapted to the cross-validation approach here by using an
#'  \emph{average cross-validation error} (\code{acve}), the \emph{mean square error} (\code{mse}) or the
#'  \emph{average standardized cross-validation error} (\code{ascve}). These last measures can only be computed in case the total number
#'  of sample points equals the number of leave-one-out covariance models. This requires to fit each cross-validation
#'  covariance model by \code{\link{prefitCV}} using the option `\code{reduce}`=\code{FALSE} which is then based on exactly
#'  one left out point. Also, we can calculate the kriging variance at the left-out sample points if we set the option `\code{type}`
#'  equal to "\code{sigK}". 
#'
#' @examples
#' data(normal)
#' 
#' # design matrix and statistics
#' X <- as.matrix(qsd$qldata[,1:2])
#' Tstat <- qsd$qldata[grep("^mean.",names(qsd$qldata))]
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
		          type = c("rmsd","msd","cve","scve","acve","mse","ascve","sigK"),
				    cl = NULL)
{		
 	stopifnot(!is.null(cvm))
    type <- match.arg(type)
	
	dx <- attr(qsd$qldata,"xdim")
	Xs <- as.matrix(qsd$qldata[seq(dx)])
	# set sample points as default
	# points to predict the CV error
	if(is.null(theta))
	 theta <- Xs
	# dataT has to be list (of class data.frame)
	dataT <- qsd$qldata[(dx+1):(dx+length(qsd$covT))]
		
	tryCatch({
			Y <- estim(qsd$covT,theta,Xs,dataT,krig.type="var")
			# cross-validation variance/RMSE of statistics
			cv <- cverrorTx(theta,Xs,dataT,cvm,Y,type,cl)
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
updateCV <- function(i, qsd, fit, ...) {
	covT <- qsd$covT
	qsd$qldata <- qsd$qldata[-i,]	# keep ith observations (k-fold CV)

	xdim <- attr(qsd$qldata,"xdim")
	Xs <- as.matrix(qsd$qldata[seq(xdim)])
	fitit <- (fit && !(nrow(Xs) %% qsd$nfit))

	cvm <- lapply(1:length(covT),
			function(j) {
				xm <- covT[[j]]					
				xm$start <- xm$param[xm$free]				
			    if(!is.null(xm$fix.nugget))
				  xm$fix.nugget <- xm$fix.nugget[-i]
				if(fitit) {
				  xm$dataT <- qsd$qldata[[xdim+j]]
				}			
				xm
			}
	)	
	
	if(fitit) {
	  res <- lapply(cvm, doREMLfit, Xs=Xs, ...)
	  if(!inherits(res,"error")) {
		 return(structure(.extractCovModels(res),"id"=i,"class"="krige"))	    
	  } else {		
		 msg <- message("Could not update covariance parameters because `REML` failed.")
	     message(msg)
		 return(.qleError(message=msg,error=res,"id"=i))
   	  } 
    } else {
	  return(structure(cvm,"id"=i,"class"="krige"))
	}
}

#' @name prefitCV 
#'
#' @title Covariance parameter estimation for cross-validation 
#'
#' @description The function constructs a list of covariance models of statistics in order to estimate the prediction error
#'  variances by a cross-validation (CV) approach at unsampled points. 
#'
#' @param qsd   	  object of class \code{\link{QLmodel}}
#' @param reduce	  if \code{TRUE} (default), reduce the number of covariance models to refit
#' @param type		  type of prediction variances, "\code{cv}" (default), see \code{\link{qle}}
#' @param control	  control arguments for REML estimation passed to \code{\link[nloptr]{nloptr}}  	
#' @param cl	      cluster object, \code{NULL} (default), of class "\code{MPIcluster}", "\code{SOCKcluster}", "\code{cluster}"
#' @param verbose	  if \code{TRUE}, print intermediate output
#'
#' @return A list of certain length depending on the current sample size (number of evaluated points).
#'  Each list element corresponds to a (reduced) number of sample points with at most \eqn{k} points
#'  (see details) left out for fitting the covariance models. 
#'
#' @details Using the CV-based approach (see vignette) for estimating the prediction variances 
#' 	might require a refit of covariance parameters of each statistic based on leaving out a certain number of sample points.
#'  The covariance models can be refitted if `\code{fit}` equals \code{TRUE} and otherwise are simply updated without fitting which
#'  saves some computational resources. The number of points left out is dynamically adjusted depending on the number
#'  of sample points in order to prevent the main estimation algorithm to fit as many models as there are points already evaluated.  
#' 
#'  For CV the number \eqn{n_c} of covariance models still to fit, that is, the number of partitioning sets of sample points, is limited by
#'  \eqn{n_c\leq n}, with maximum \eqn{k} sampling points deleted from the full sample set with overall \eqn{n} sample points such that
#'  \eqn{n=n_c k} (see vignette for further details). 
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
prefitCV <- function(qsd, reduce = TRUE, type = c("cv","max"),
		              control = list(),	cl = NULL, verbose = FALSE)
{	
	N <- nrow(qsd$qldata)
	p <- if(reduce) {
			ifelse(N>20,
			 ifelse(N>30,
			   ifelse(N>40,
				ifelse(N>50,
				 ifelse(N>200,0.1,0.3),0.4),0.6),0.8),1.0)
		 } else 1
	nb <- floor(p*N)
	k <- ceiling(N/nb) # block size
	S <-
	 if((N-k) >= qsd$minN){
        Ni <- seq_len(N) 
		split(Ni, sort(Ni%%nb))
	 } else stop(paste0("Total number of points must be at least ",qsd$minN," for cross-validation."))
	
	fit <- isTRUE(qsd$cv.fit)
    type <- match.arg(type)
	# Leave-k-Out CV
	tryCatch({			 
		 if(length(control) > 0L) {		
			opts <- nloptr::nl.opts()
			opts[names(control)] <- control
		 } else {
			opts <- attr(qsd,"opts")		
		 }			
		 return(
		   structure(doInParallel(S, updateCV, qsd=qsd, fit=fit, 
						opts=opts, cl=cl, verbose=verbose),
	        type=type)
		 )

	  },error = function(e) {
		 msg <- paste0("Prefitting covariance models failed.\n")
		 if(verbose)
		   message(msg)
		 stop(e)
	  }
	)	
}

# internal, alloc C structure
# Also for QL, a pre-set (inverse) variance matrix can be supplied by VTX
# No predictions variances here (done at C level), theta is only needed
# for the weighted version of avergage variance approximation
.qdAlloc <- function(qsd, Sigma = NULL, ..., inverted = FALSE, cvm = NULL) {	
	X <- as.matrix(qsd$qldata[seq(attr(qsd$qldata,"xdim"))])
	useSigma <- (!is.null(Sigma) && qsd$var.type == "const")
		
	if(qsd$var.type != "kriging" && is.null(Sigma)){
		if(qsd$var.type %in% c("wcholMean","wlogMean")){
			nms <- names(list(...))
			if(!all( c("W","theta") %in% nms))
			 message(paste0("Found `var.type`=\"",qsd$var.type, "\" but no weighting matrix `W` or estimate `theta` was supplied!."))		
		}
		Sigma <- covarTx(qsd,...,cvm=cvm)[[1]]$VTX	
	} else if(useSigma && !inverted){
		# Only for constant Sigma, which is used as is!
		Sigma <- try(gsiInv(Sigma),silent=TRUE)
		if(inherits(Sigma,"try-error")) {
			msg <- paste0("Inversion of constant variance matrix failed.")
			message(msg)
			return(.qleError(message=msg,error=Sigma))
		}
	}		
	# init QL data and kriging models	
	qlopts <- list("varType"=qsd$var.type,
				   "useCV"=!is.null(cvm),
				   "useSigma"=useSigma)
	
	 # return TRUE for success othewise signal error
	try(.Call(C_initQL,qsd,qlopts,X,Sigma,cvm))	
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
		
		# even for `qle` we can use a kind of constant Sigma but do not need to
		# invert it. In this case prediction variances are always used at C level
		# Sigma is inverted after adding these as diagonal terms 	  	 
		if(qsd$var.type == "kriging"){
			stop("`Sigma` must be `NULL` if using kriging approximation of variance matrix.")	    
		} else if(qsd$var.type == "const" && qsd$criterion == "qle")
			stop("`Sigma` cannot be used as a constant variance matrix for criterion `qle`.")			
				
	} else if(qsd$var.type == "kriging" && is.null(qsd$covL))
		stop("Covariance models for kriging variance matrix must be given, see function `setQLdata`.")	
	  else if(qsd$var.type == "const") 
		stop("`Sigma` must not be NULL for `const` variance matrix approximation.")
	  	
}


#' @name searchMinimizer
#'
#' @title Minimize a criterion function 
#'
#' @description The function searches for a root of the quasi-score vector or minimizes one of the criterion functions.
#'
#' @param x0		  (named) numeric vector, the starting point
#' @param qsd   	  object of class \code{\link{QLmodel}}
#' @param method	  names of possible minimization routines (see details) 
#' @param opts		  list of control arguments for quasi-scoring iteration, see \code{\link{qscoring}}
#' @param control 	  list of control arguments passed to the auxiliary routines
#' @param ...		  further arguments passed to \code{\link{covarTx}}
#' @param obs		  numeric vector of observed statistics, overwrites `\code{qsd$obs}`
#' @param info		  additional information at found minimizer
#' @param check		  logical, \code{TRUE} (default), whether to check input arguments
#' @param pl		  numeric value (>=0), the print level 
#' @param verbose	  if \code{TRUE} (default), print intermediate output
#'
#' @details The function provides an interface to local and global numerical minimization routines
#'  using the approximate quasi-deviance (QD) or Mahalanobis distance (MD) as an objective function.
#'  
#'  The function does not require additional simulations to find an approximate minimizer or root. The
#'  numerical iterations always take place on the fast to evaluate criterion function approximations.
#'  The main purpose is to provide an entry point for minimization without the
#'  need of sampling new candidate points for evaluation. This is particularly useful if we search
#'  for a "first-shot" minimizer. 
#' 
#'  The criterion function is treated as a deterministic (non-random) function during minimization
#'  (or root finding) whose surface depends on the sample points. Because of the typical nonconvex nature of the
#'  criterion functions one cannot expect a global minimizer by applying any local search method like,
#'  for example, the scoring iteration \code{\link{qscoring}}.
#'  Therfore, if the scoring iteration or some other available method gets stuck in a possibly local
#'  minimum of the criterion function showing at least some kind of numerical convergence we use such
#'  minimizer as it is and finish the search, possibly being unlucky, having not found an approximate root
#'  of the quasi-score vector (or minimum of the Mahalanobis distance). If there is no convergence practically,
#'  the search is restarted by switching to the next user supplied minimization routine defined in `\code{method}`. 
#' 
#'  \subsection{Choice of auxiliary minimization methods}{  
#'  Besides the local quasi-scoring (QS) iteration, `\code{method}` equal to "\code{qscoring}", the following
#'  (derivative-free) auxiliary methods from the \code{\link[nloptr]{nloptr}} package are available for minimizing
#'  both criterion functions:
#'  
#' 	\itemize{
#' 	  \item{}{ \code{\link[nloptr]{bobyqa}}, \code{\link[nloptr]{cobyla}} and \code{\link[nloptr]{neldermead}}}
#'    \item{}{ \code{\link[nloptr]{direct}}, global search with a locally biased version named \code{directL}}
#' 	  \item{}{ \code{\link[nloptr]{lbfgs}},  for minimizing the MD with constant `\code{Sigma}` only}
#' 	  \item{}{ \code{\link[nloptr]{nloptr}}, as the general optimizer, which allows to use further methods}
#'  }
#'    
#'  Using quasi-scoring first, which is only valid for minimizing the QD function, is always a good idea since we might have done
#'  a good guess already being close to an approximate root. If it fails we switch to any of the above alternative methods
#'  (e.g. \code{\link[nloptr]{bobyqa}} as the default method) or eventually - in some real hard situations - to the
#'  method `\code{direct}` or its locally biased version `\code{directL}`. The order of processing is determined
#'  by the order of appearance of the names in the argument `\code{method}`. Any method available from package `\code{nloptr}` can be
#'  chosen. In particular, setting \code{method="nloptr"} and `\code{control}` allows to choose a multistart algorithm such
#'  as \code{\link[nloptr]{mlsl}}.
#' 
#'  Only if there are reasonable arguments against quasi-scoring, such as expecting a local
#'  minimum rather than a root first or an available limited computational budget, we can always apply
#'  the direct search method `\code{direct}` leading to a globally exhaustive search. Note that we must always supply a starting
#'  point `\code{x0}`, which could be any vector valued parameter of the parameter space unless method `\code{direct}` is
#'  chosen. Then `\code{x0}` is still required but ignored as a starting point since it uses the "center point" of
#'  the (hyper)box constraints internally. In addition, if CV models `\code{cvm}` are given, the CV based prediction variances
#'  are inherently used during consecutive iterations of all methods. This results in additional computational efforts
#'  due to the repeated evaluations of the statistics to calculate these variances during each new iteration.  
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
#' @seealso \code{\link[nloptr]{nloptr}}, \code{\link{qscoring}}
#' 			
#' @rdname searchMinimizer
#' @author M. Baaske
#' @export
#' @importFrom nloptr direct directL cobyla bobyqa lbfgs neldermead
searchMinimizer <- function(x0, qsd, method = c("qscoring","bobyqa","direct"),
					 opts = list(), control = list(), ...,  
					   obs = NULL, info = TRUE, check = TRUE, 
					     pl = 0L, verbose = FALSE)
{
	if(check)
	 .checkArguments(qsd,x0,...)
    stopifnot(is.numeric(pl) && pl >= 0L )
	
	x0 <- 
	 if(is.matrix(x0))
		structure(as.numeric(x0),names=colnames(x0))	
	 else unlist(x0)
 	
    fun.name <- ""
	nms <- names(x0)	
	# current sample points
	xdim <- attr(qsd$qldata,"xdim")
	if(xdim != length(x0))
	 stop("Dimension of `x0` does not match.")
	
	# may overwrite (observed) statistics	
	if(!is.null(obs)) {
		obs <- unlist(obs)
		if(anyNA(obs) | any(!is.finite(obs)))
			warning("`NA`, `NaN` or `Inf` values detected in argument `obs`.")
		if(!is.numeric(obs) || length(obs)!=length(qsd$covT))
		  stop("Object `obs` must be a (named) numeric vector or list of length equal to the number of given statistics in `qsd`.")
		qsd$obs <- obs
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
		    qscoring(qsd,x0,opts,...,check=FALSE,pl=pl,verbose=verbose)
		   }, error = function(e) {	e }
  		 )
		} else NULL
	}
	
    if(!is.null(S0) && (.isError(S0) || S0$convergence < 0L)){
	   if(pl > 0L) { 
		 msg <- .makeMessage("Minimization by `",fun.name,"` did not converge: ")
		 if(!is.null(S0$convergence))
		  msg <- c(msg, paste0(" (status=",S0$convergence,")") )
	  	 if(inherits(S0,"error"))
			msg <- c(msg, conditionMessage(S0)) 
		 message(msg)
	   }
	   if(pl >= 10L){
	   	   message("Failed minimization: \n\n")
		   print(S0)
		   cat("\n\n")
	   }
	   method <- method[-1]
	   if(is.na(method[1])){
			message("No convergence and only one method supplied.")
			return(S0)	
	   }		
    }
	
	if(is.null(S0) || S0$convergence < 0L) {	  	
	  S0 <- 
		tryCatch({			
			if(length(control) == 0L){
			  control <- list("stopval"=0,"maxeval"=1000,
							  "ftol_rel"=1e-7,"xtol_rel"=1e-6)		  	  	
	  		}			
			# alloc C level
			if(!.qdAlloc(qsd,...))
			 stop("Could not allocate C memory and construct QL model.")
			
			fn <-
			 switch(qsd$criterion,
				"qle" = { function(x) .Call(C_qDValue,x) },
				"mahal" = { function(x) .Call(C_mahalValue,x) }
			)			
		 	repeat {
				if(!is.na(method[1])) {
					if(pl > 0L)
					  cat(paste0("Using method: ",method[1],"...\n"))
					fun.name <- method[1]					
				} else {
					return(.qleError(message="No convergence and only one method supplied: ",
							call = sys.call(),
							   error = if(inherits(S0,"error")) conditionMessage(S0) else NULL,
							   	S0=S0, method = method[1]))	
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
									if(qsd$criterion != "mahal" || qsd$var.type == "kriging")
									  stop("`lbfgs` only for criterion `mahal` using a constant `Sigma` or an average variance approximation.")
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
									if(is.null(control$algorithm)){
										control["algorithm"] <- "NLOPT_LN_BOBYQA"
										message(paste0("Using default derivative-free method: ",control$algorithm))									
									}
									ret <- do.call(nloptr::nloptr, list(x0, eval_f=fn, lb=qsd$lower,
													ub=qsd$upper, opts=control))
									structure(list("par"=ret$solution,
												   "value"=ret$objective,
												   "iter"=ret$iterations,
												   "convergence"=ret$status,
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
					  }, error = function(e) {e})
			    
				if(!inherits(S0,"error") && S0$convergence >= 0L) {				
					break
				} else {
					msg <- .makeMessage("Minimization failed by: ",fun.name,".")
					message(msg, if(inherits(S0,"error")) conditionMessage(S0) else "",sep=" ")
				  	method <- method[-1]
				}
			}
			S0
		}, error = function(e) {
			 msg <- .makeMessage("Surrogate minimization failed: ",
					  conditionMessage(e))
			 message(msg)
			 return(.qleError(message=msg,call=sys.call(),error=e,method=fun.name))			
		}, finally = { 
			 if(!.qdDealloc())
			   stop("Could not release C memory.")
		})	
	}
	if(.isError(S0))
	  return(S0)		
	if(!is.null(nms))
 	  names(S0$par) <- nms     
 	
    if(class(S0) != "QSResult") {	 
	  S0 <- structure(
	    	    c(S0,list("method"=fun.name,				   	  
						  "criterion"=qsd$criterion,						 
				 		  "start"=x0)),
	  		   class="QSResult")
		 
	  if(info){
		qd <-
		  tryCatch({				
				if(qsd$criterion == "mahal")
					mahalDist(S0$par,qsd,...,check=FALSE,verbose=verbose)
				else
					quasiDeviance(S0$par,qsd,...,check=FALSE,verbose=verbose)
			}, error = function(e) {
				 msg <- .makeMessage("Error in criterion function: ",
						   conditionMessage(e))				 
				.qleError(message=msg,call=sys.call(),error=e)		
		  })
		if(!.isError(qd)){			
	 		S0 <- structure(
					  c(S0,qd[[1]][which(!(names(qd[[1]]) %in% names(S0)))]),
					 Sigma = attr(qd,"Sigma"),
				   class = "QSResult")				 	
	 	} else { 
			message(qd$message)
			return(structure(S0, error = qd))
		}
	  }
    }	
	if(verbose){
	  cat(paste0("Successful minimization by: ",fun.name," (status=",S0$convergence,")","\n\n"))
	  if(pl >= 10L){
		  print(S0)
		  cat("\n\n")
	  }
    }  		
    return(S0)   
}

#' @name multiSearch
#'
#' @title A multistart version of local searches for parameter estimation
#'
#' @description  The function is multistart version of \code{\link{searchMinimizer}} which selects the best
#' 	root of the quasi-score (if there is any) or a local minimum from all found minima according to the criteria described in the vignette.
#' 
#' @param xstart 	numeric, \code{NULL} default, list, vector or matrix of starting parameters
#' @param qsd		object of class \code{\link{QLmodel}}
#' @param ...    	arguments passed to \code{\link{searchMinimizer}} 
#' @param nstart 	number of random samples from which to start local searches (if `\code{xstart}`=\code{NULL}, then ignored)
#' @param optInfo 	logical, \code{FALSE} (default), whether to store original local search results
#' @param cl 	 	cluster object, \code{NULL} (default), of class "\code{MPIcluster}", "\code{SOCKcluster}", "\code{cluster}"
#' @param verbose	if \code{TRUE} (default), print intermediate output
#' 
#' @details The function performs a number of local searches depending which local method `\code{method}` was passed to
#'  \code{\link{searchMinimizer}}. Either the starting points are given by `\code{xstart}` or are generated as an augmented 
#'  design based on the sample set stored in `\code{qsd}`. The function evaluates all found solutions and selects the one which 
#'  is best according to the criteria defined in the vignette.
#' 
#' @return Object of class \code{QSResult} and attribute `\code{roots}`, e.t. the matrix of estimated parameters for which any of
#'  the available minimization methods has been successfully applied. If `code{optInfo}` is \code{TRUE}, then the originally estimtation reuslts
#'  are also returned. The best solution is stored as an attribute named `\code{par}` if any could have been found.
#' 
#' @seealso \code{\link{checkMultRoot}}
#'  
#' @examples 
#'  data(normal)
#'  x0 <- c("mu"=3.5,"sigma"=1.5)
#'  S0 <- multiSearch(xstart=x0,qsd,method=c("qscoring","bobyqa"),
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
multiSearch <- function(xstart=NULL, qsd, ..., nstart=10, optInfo=FALSE, cl=NULL, verbose=FALSE){
	if(nstart>0L){		
		X <- as.matrix(qsd$qldata[seq(attr(qsd$qldata,"xdim"))])
		Xs <- try(multiDimLHS(N=nstart,qsd$lower,qsd$upper,X=X,
				   method="augmentLHS",type="matrix"),silent=TRUE)
   		if(inherits(Xs,"try-error")) 
		  message("Could not generate random starting points in function `multiDimLHS`.")
	    else if(!is.null(xstart))
		 xstart <- rbind(xstart,Xs)
	    else xstart <- Xs
	} 
	if(is.null(xstart))
	 stop("No starting points given for local searches.")
	if(!is.list(xstart))
	 xstart <- .ROW2LIST(xstart)
	
	opt.args <- list(...)
	RES <- do.call(doInParallel,
			c(list(X=xstart,
				FUN=function(xstart,...){
					searchMinimizer(xstart,...)
				},
				cl=cl, qsd=qsd), opt.args)) 
	
	if(.isError(RES))
	 return(RES)
 	# do not evaluate solution for just a single parameter
    if(length(RES) == 1L){
		if(verbose)
		  message("We do not compare or evaluate the root criteria when only a single solution is available.")
		return (RES[[1]])
	}
	# check results again
	ok <- which(sapply(RES,function(x) !.isError(x) & x$convergence >= 0L))
	if(length(ok) == 0L){
		msg <- .makeMessage("All local searches have errors or did not converge.")
		message(msg)
		return(.qleError(message=msg,call=match.call(),error=RES))							
	} else if(length(ok) < length(RES)){
		message(paste0("A total of ",length(RES)-length(ok)," local searches have errors or did not converge."))							
	}
	
	hasError <- which(!(1:length(RES) %in% ok))
	if(length(hasError) > 0L)	
	 message(paste0("A total of ",length(hasError)," local searches failed."))
		
	roots <- .evalRoots(RES[ok])
	if(.isError(roots)) {
		msg <- .makeMessage("Could not evaluate best results of local searches")
		message(msg)
		attr(roots,"optInfo") <- RES
		return(.qleError(message=msg,call=match.call(),error=roots))	   
	}
	id <- attr(roots,"id")
	if(anyNA(id)){
		msg <- .makeMessage("Could not find any root.")
		message(msg)
		attr(roots,"optInfo") <- RES
		return(.qleError(message=msg,call=match.call(),error=roots))
 	}
	stopifnot(length(id)==1L)
	structure(RES[[id]],
		"roots"=if(optInfo) roots else NULL,
		"optRes"=if(optInfo) RES else NULL,				
		"hasError"=hasError) 	
}

#' @name qle
#'
#' @title Simulated quasi-likelihood parameter estimation
#'
#' @description  This is the main function of the simulated quasi-likelihood estimation (QLE) approach. 
#' 
#' @param qsd			object of class \code{\link{QLmodel}}
#' @param sim		    simulation function, see details
#' @param ...			further arguments passed to the simulation function `\code{sim}` 
#' @param nsim			optional, number of simulation replications at each new sample point,
#'  					`\code{qsd$nsim}` (default)
#' @param x0 			optional, numeric vector of starting parameters
#' @param obs			optional, numeric vector of observed statistics, overwrites `\code{qsd$obs}`
#' @param Sigma			optional, constant variance matrix estimate of statistics (see details) 
#' @param global.opts	options for global search phase
#' @param local.opts	options for local search phase
#' @param method		vector of names of local search methods	
#' @param qscore.opts   list of control arguments passed to \code{\link{qscoring}}
#' @param control		list of control arguments passed to any of the routines defined in `\code{method}` 
#' @param errType		type of prediction variances, choose one of "\code{kv,cv,max}" (see details)
#' @param multistart	logical, \code{FALSE} (default), whether to search for local minimia or roots from multiple starting points at global phase  
#' @param pl			print level, use \code{pl}>0 to print intermediate results
#' @param cl			cluster object, \code{NULL} (default), of class "\code{MPIcluster}", "\code{SOCKcluster}", "\code{cluster}" 
#' @param iseed			integer seed, \code{NULL} (default) for default seeding of the random number generator (RNG) stream for each worker in the cluster
#' @param plot 			if \code{TRUE}, plot newly sampled points (for 2D-parameter estimation problems only)
#'
#' @return List of the following objects:
#' 	  \item{par}{ final parameter estimate}
#' 	  \item{value}{ value of criterion function}
#'    \item{ctls}{ a data frame with values of stopping conditions}
#'    \item{qsd}{ final \code{\link{QLmodel}} object, including all sample points
#' 				  and covariance models}
#' 	  \item{cvm}{ CV fitted covariance models}
#'    \item{why}{ names of stopping conditions matched}
#'	  \item{final}{ final local minimization results of the criterion function, see \code{\link{searchMinimizer}} }
#'	  \item{score}{ quasi-score vector or gradient of the Mahalanobis distance}
#' 	  \item{convergence}{ logical, whether the iterates converged, see details} 	  
#' 
#'  Attributes: 	 
#'  
#'  \item{tracklist}{ an object (list) of class \code{QDtrack} containing the local minimization results,
#'     evaluated sample points and the status of the corresponding iteration}    
#'  \item{optInfo}{ a list of arguments related to the estimation procedure:}
#'  \itemize{
#'    \item{x0:}{ starting parameter vector}
#' 	  \item{W:}{ final weighting matrix (equal to quasi-information matrix at \code{theta}) used for both variance
#' 			 average approximation, if applicable, and as the predicted variance for (local) sampling of new candidate points
#' 			 according to a multivariate normal distribution with this variance and the current root as the mean parameter.}
#'    \item{theta:}{ the parameter corresponding to \code{W}, typically an approximate root or local minimzer} 
#' 	  \item{last.global:}{ logical, whether last iteration sampled a point globally}
#' 	  \item{minimized:}{ whether last local minimization was successful}
#' 	  \item{useCV:}{ logical, whether the CV approach was applied}
#' 	  \item{method:}{ name of final search method applied}
#'    \item{nsim:}{ number of simulation replications at each evaluation point}
#' 	  \item{iseed}{ the seed to initialize the RNG}
#'  }
#'     
#' @details
#'  The function sequentially estimates the unknown model parameter. Basically, the user supplies a simulation function `\code{sim}`
#'  which must return a vector of summary statistics (as the outcome of model simulations) and expects a vector of parameters
#'  as its first argument. Further arguments can be passed by the `\code{\ldots}` argument. The object
#'  `\code{qsd}` aggregates the type of variance matrix approximation, the data frame of observed and simulated data, the
#'  initial sample points and the covariance models of the involved statistics (see \code{\link{QLmodel}}). In addition, it defines
#'  the criterion function by `\code{qsd$criterion}`, which is either used to monitor the sampling process or minimized itself. The user
#'  also has the option to choose among different types of prediction variances: either "\code{kv}" (kriging variances), "\code{cv}"
#'  (CV variances) or the maximum of both, by "\code{max}", are available.
#' 
#'  \subsection{Criterion functions}{The QD criterion function follows the quasi-likelihood estimation principle (see vignette)
#'  and seeks a solution to the quasi-score equation. Besides, the Mahalanobis distance (MD) as an alternative (simulation-based)
#'  criterion function has a more direct interpretation. It can be seen as a (weighted or generalized) least squares criterion
#'  depending on the employed type of variance matrix approximation. For this reason, we support several types of variance matrix
#'  approximations. In particular, given `\code{Sigma}` and setting `\code{qsd$var.type}` equal to "\code{const}" treats `\code{Sigma}`
#'  as a constant estimate throughout the whole estimation procedure. Secondly, if `\code{Sigma}` is supplied and used as
#'  an average variance approximation (see \code{\link{covarTx}}), it is considered an initial variance matrix approximation and
#'  recomputed each time an approximate (local) minimizer of the criterion function is found. This is commonly known as an iterative update
#'  strategy of variance matrices in the context of GMM estimation. Opposed to this, setting `\code{qsd$var.type}` equal to
#'  "\code{kriging}" corresponds to continuously updating the variance matrix each time a new criterion function value is
#'  required at any point of the parameter space. In this way the algorithm can also be seen as a simulated version of a least squares
#'  method or even as a special case of a \emph{simulated method of moments} (see, e.g. [3]). Note that some input combinations
#'  concerning the variance approximation types are not applicable since the criterion "\code{qle}", which uses the
#'  QD criterion function, does not support a constant variance matrix at all.
#'  }
#'       
#'  \subsection{Monte Carlo (MC) hypothesis testing}{ The algorithm sequentially evaluates promising local minimizers of the criterion function during
#'  the local phase in order to assess the plausibility of being an approximate root of the corresponding quasi-score vector. We use essentially
#'  the same MC test procedure as in \code{\link{qleTest}}. First, having found a local minimum of the test statistic, i.e. the criterion
#'  function, given the data, new observations are simulated w.r.t. to the local minimizer and the algorithm re-estimates the approximate roots for each
#'  observation independently. If the current minimizer is accepted as an approximate root at the significance level `\code{local.opts$alpha}`, then the algorithm stays
#'  in its local phase and continues sampling around the current minimizer accoring to its asymptotic variance (measured by the inverse of the
#'  predicted quasi-information) and uses the additional simulations to improve the current kriging approximations. Otherwise we switch to the global phase and
#'  do not consider the current minimizer as an approximate root.
#' 
#'  This procedure also allows for a stopping condition derived from the reults of the MC test. We can compare the estimated mean squared error (MSE) with the
#'  predicted error of the approximate root by its relative difference and terminate in case this value drops below a user-defined bound `\code{perr_tol}`
#'  (either as a scalar value or numeric vector of length equal to the dimension of the unknown parameter). A value close to zero suggests a good match of both
#'  error measures. The testing procedure is disabled by default. Use `\code{local.opts$test=TRUE}` for testing approximate roots. A value of the criterion function smaller
#'  than `\code{local.opts$ftol_abs}` indicates that the corresponding minimizer could be an approximate root. Otherwise the last evaluation point is used as
#'  a starting point for next local searches which mimics a random multistart type minimization over the next iterations of the algorithm. This behaviour is
#'  also implemented for results of the above MC test when the local minimizer is not accepted as an approximate root. Note that this approach has the
#'  potential to escape regions where the criterion function value is quite low but, however, is not considered trustworthy in relation to the upper bound
#'  `\code{local.opts$ftol_abs}` or the results of the MC test procedure.
#'  
#'  If one of the other termination criteria is met in conjunction with a neglectable value of the criterion function, we
#'  say that the algorithm successfully terminated and converged to a local minimizer of the criterion function which could be an approximate root of the quasi-score
#'  vector. We then can perform a goodness-of-fit test in order to assess its plausibility (see \code{\link{qleTest}}) and quantify the empirical and predicted
#'  estimation error. If we wish to improve the final estimate the algorithm allows for a simple warm start strategy though not yet as an fully automated
#'  procedure. The algorithm can be easily restarted based on the final result of the preceeding run. We only need to extract the object
#'  `\code{OPT$qsd}` as an input argument to function \code{\link{qle}} again. 
#'  }
#' 
#'  \subsection{Sampling new points}{Our QLE approach dynamically switches from a \emph{local} to a \emph{global search phase} and vise versa for
#'  sampling new promising candidates for evaluation, that is, performing new simulations of the statistical model. Depending on the current value of the criterion
#'  function three different sampling criteria are used to select next evaluation points which aim on potentially improving the quasi-score
#'  or criterion function approximation. If a local minimizer of the criterion function has been accepted as an approximate root, then a local search
#'  tries to improve its accuracy. The next evaluation point is either selected according to a weighted minimum-distance criterion (see [2] and vignette),
#'  for the choice `\code{nextSample}` equal to "\code{score}", or by maximizing the weighted variance of the quasi-score vector in
#'  case `\code{nextSample}` is equal to "\code{var}". In all other cases, for example, if identifiable roots of the QS could not be found
#'  or the (numerical) convergence of the local solvers failed, the global phase of the algorithm is invoked and selects new potential
#'  candidates accross the whole search space based on a weighted selection criterion. This assigns large weights to candidates
#'  with low criterion function values and vise versa. During both phases the cycling between local and global candidates is
#'  controlled by the weights `\code{global.opts$weights}` and `\code{locals.opts$weights}`, respectively. Besides this, the smaller
#'  the weights the more the candidates tend to be globally selected and vice versa during the global phase. Within the local phase,
#'  weights approaching one result in selecting candidates close to the current minimizer of the criterion
#'  function. Weights approximately zero maximize the minimum distance between candidates and previously sampled points and
#'  thus densely sample the search space almost everywhere if the algorithm is allowed to run infinitely. The choice of weights
#'  is somewhat ad hoc but may reflect the users preference on guiding the whole estimation more biased towards either a local
#'  or global search. In addition the local weights can be dynamically adjusted if `\code{useWeights}` is \code{FALSE}
#'  depending on the current progress of estimation. In this case the first weight given by `\code{locals.opts$weights}` is 
#'  initially used for this kind of adjustment.   
#'  }
#' 
#'  Some notes: For a 2D parameter estimation problem the function can visualize the sampling and selection process, which
#'  requires an active 2D graphical device in advance. The function can also be run in an cluster environment
#'  using the `\code{parallel}` package. Make sure to export all functions to the cluster environment `\code{cl}` beforehand,
#'  loading required packages on each cluster node, which are used in the simulation function
#'  (see \code{\link{clusterExport}} and \code{\link{clusterApply}}).
#'  If no cluster object is supplied, a local cluster is set up based on forking (under Linux) or as a socket connection
#'  for other OSs. One can also set an integer seed value `\code{iseed}` to initialize each worker, see \code{\link{clusterSetRNGStream}},
#'  for reproducible results of estimation in case a local cluster is used, i.e. \code{cl=NULL} and option \code{mc.cores>1}. If
#'  using a prespecified cluster object \code{cl}, then the user is responsible for seeding whereas the seed can be stored
#'  in the return value as well, see attribute `\code{optInfo}$iseed`.  
#' 
#'  The following controls `\code{local.opts}` for the local search phase are available:
#'   \itemize{
#'   \item{\code{ftol_rel}:}{ upper bound on relative change in criterion function values}
#'   \item{\code{lam_max}:}{ upper bound on the maximum eigenvalue of the generalized eigenvalue decomposition of
#' 		the quasi-information matrix and estimated interpolation error (variance) of quasi-score.
#'  	This stops the main iteration sampling new locations following the idea that in this case
#' 		the quasi-score interpolation error has dropped below the estimated precision at best measured by
#' 		quasi-information matrix for `\code{global.opts$NmaxLam}` consecutive iterates.}
#' 	 \item{\code{pmin}:}{ minimum required probability that a new random candidate sample falls inside the parameter
#'                space. Dropping below this value triggers a global phase sampling step. This might indicate
#' 				  that the inverse of the quasi-information matrix does not sufficiently reflect the variance
#' 				  of the current parameter estimate due to a sparse sample or the (hyper)box constraints of the
#' 				  parameter space could be too restrictive.}
#' 	 \item{\code{nsample}:}{ sampling size of candidate locations at the local phase}
#' 	 \item{\code{weights}:}{ vector of weights, \eqn{0\leq\code{weights}\leq 1}, for local sampling}
#' 	 \item{\code{useWeights}:} {logical, if \code{FALSE} (default), dynamically adjust the weights, see vignette}
#'	 \item{\code{ftol_abs}:}{ upper bound on the function criterion: values smaller trigger the local phase
#'    treating the current minimzer as an approximate root otherwise forces the algorithm to switch to the global phase and vice versa.}
#'   \item{\code{eta}:}{ values for decrease and increase of the local weights, which is intended to faciliate convergence
#' 		 while sampling new points more and more around the current best parameter estimate.} 
#' 	 \item{\code{alpha}:}{ significance level for computation of empirical quantiles of one of the test statistics, that is,
#'          testing a parameter to be a	root of the quasi-score vector in probability.}
#'   \item{perr_tol}{ upper bound on the relative difference of the empirical and predicted error of an approximate root}
#'   \item{\code{nfail}:}{ maximum number of consecutive failed iterations}
#'   \item{\code{nsucc}:}{ maximum number of consecutive successful iterations}
#'   \item{\code{nextSample}:}{ either "\code{score}" (default) or "\code{var}" (see details)} 
#'   }
#' 
#'  The following controls `\code{global.opts}` for the global search phase are available:   
#' 	\itemize{
#'   \item{\code{stopval}:}{ stopping value related to the criterion function value, the main iteration terminates
#' 				     as soon as the criterion function value drops below this value. This might be preferable to a time consuming
#' 					 sampling procedure if one whishes to simply minimize the criterion function or find a first
#' 					 approximation to the unknown model parameter.}
#'   \item{\code{C_max}:}{ upper bound on the relative maximum quasi-score interpolation error. The algorithm terminates
#' 					its value drops below after a number of `\code{global.opts$NmaxCV}` consecutive iterations.}
#' 	 \item{\code{xtol_rel}:}{ relative change of found minimizer of the criterion function or root of quasi-score.}
#' 	 \item{\code{maxiter}:}{ maximum allowed global phase iterations }
#' 	 \item{\code{maxeval}:}{ maximum allowed global and local iterations }
#' 	 \item{\code{sampleTol}:}{ minimum allowed distance between sampled locations at global phase}	
#' 	 \item{\code{weights}:}{ vector of \eqn{\code{weights}>0} for global sampling}
#'   \item{\code{nsample}:}{ sampling size of candidate locations at the global phase}
#'   \item{\code{NmaxRel}:}{ maximum number of consecutive iterates until stopping according to `\code{xtol_rel}`}
#'   \item{\code{NmaxCV}:}{ maximum number of consecutive iterates until stopping according to `\code{C_max}`}
#'   \item{\code{NmaxSample}:}{ maximum number of consecutive iterations until stopping according to `\code{sampleTol}`}
#'   \item{\code{NmaxLam}:}{ maximum number of consecutive iterations until stopping for which the generalized eigenvalue of the variance
#' 		 of the quasi-score vector within the kriging approximation model and its total variance measured by the quasi-information matrix
#'       at some estimated parameter drops below the upper bound `\code{local.opts$lam_max}` }
#'   \item{\code{NmaxQI}:}{ maximum number of consecutive iterations until stopping for which the relative difference of the empirical error
#'      and predicted error of an approximate root drops below `\code{perr_tol}`}
#'   \item{\code{Nmaxftol}:}{ maximum number of consecutive iterations until stopping for which the relative change in the values
#'    of the criterion function drops below `\code{local.opts$ftol_rel}`}
#'  }
#'  
#' 
#' @seealso \code{\link{mahalDist}}, \code{\link{quasiDeviance}}, \code{\link{qleTest}} 
#' 
#' @examples
#' data(normal)
#'  
#' # main estimation with new evaluations
#' # (simulations of the statistical model)
#' OPT <- qle(qsd,qsd$simfn,nsim=10,
#' 		    global.opts=list("maxeval"=1),
#'  		local.opts=list("test"=FALSE))
#' 
#' 
#' @author M. Baaske
#' @rdname qle
#' 
#' @useDynLib qle, .registration = TRUE, .fixes = "C_"
#' @export 
#' 
#' @import parallel stats
#' @importFrom graphics points
qle <- function(qsd, sim, ... , nsim, x0 = NULL, obs = NULL,
		        Sigma = NULL, global.opts = list(), local.opts = list(),
				  method = c("qscoring","bobyqa","direct"),
				   qscore.opts = list(), control = list(),
				    errType = "kv", multistart=FALSE, pl = 0, 
					 cl = NULL, iseed = NULL, plot=FALSE)
{		
	# print information 	
	.printInfo = function(){		
		if(pl > 0L) {
			cat("\n")
			nmax <- nglobal+nlocal
			if(nmax >= maxEval || nglobal >= maxIter)
			  cat("Final results: \n\n")
   		    cat("Total evaluations...",nmax,"\n")						
			cat("Criterion value.....",formatC(ft, digits=4, format="e", big.mark=","),"\n\n")
			cat("Current iterate: \n\n")
			print.default(formatC(signif(as.numeric(xt), digits=8), digits=8, format="fg", flag="#"),
					print.gap = 4, quote = FALSE)
			cat("\n")						
			if(!is.null(Stest) && !.isError(Stest)){
			 qt <- attr(Stest$test,"qt")
			 cat("Criterion value < ",names(qt),"quantile: ",formatC(ft, digits=4, format="e", big.mark=","))
			 cat(paste0(" <",formatC(qt, digits=4, format="e", big.mark=",")),"\n")			 
			}			 
			cat("\n")
		}		
	}

	.showConditions = function() {
		if(pl > 1L) {
		    cat("Iterations......",paste0("global=",nglobal,", local=",nlocal,"\n"))
			cat("Sampling:.......",paste0(if(status[["global"]]>1L) "global" else "local", " (status=",status[["global"]],")\n"))
			cat("Local search:...",paste0(ifelse(status[["minimized"]],"success","failed"),"\n"))			
			if(locals$nextSample=="score")
				cat("weight factor:..",w,"\n")
			cat("\n")
			df <- as.data.frame(
					cbind(
					 c(formatC(signif(as.numeric(x),digits=6),digits=6,format="fg", flag="#"),formatC(signif(f,digits=4),digits=4,format="e")),
					 c(formatC(signif(as.numeric(xt),digits=6),digits=6,format="fg", flag="#"),formatC(signif(ft,digits=4),digits=4,format="e")),
					 c(formatC(signif(as.numeric(Snext$par),digits=6),digits=6,format="fg", flag="#"),formatC(signif(Snext$value,digits=4),digits=4,format="e"))))
			dimnames(df) <- list(c(names(x0),"value"),c("Start","Estimate", "Sample"))
			print(format(df, digits=6),
				print.gap = 2, right=FALSE, quote = FALSE)	
			cat("\n\n")		
			# other conditions
			# max of quasi-score depends on whether criterion was minimized (local) or not
			cat("Current stopping conditions: \n\n")			
			cond <- 
			 if(status[["minimized"]]) {
			   	c("|score_max|" = max(abs(S0$score)))				
			 } else c("|score_max|" = max(abs(Snext$score)))
	 
			cond <-
			 c(cond,
			   "lam_max"=unlist(ctls["lam_max","val"]),
			   "varTol"=unlist(ctls["C_max","val"]),
			   "ftol_rel"=unlist(ctls["ftol_rel","val"]),
			   "xtol_rel"=unlist(ctls["xtol_rel","val"]),
			   "sampleTol"=unlist(ctls["sampleTol","val"]))	
			
		 	print.default(format(cond, digits = 4, justify="left"),
							print.gap = 2, quote = FALSE)
								
			if(pl > 2L) {
				if(!is.null(Stest) && !.isError(Stest)){
				   cat("\n\n")
				   cat("MC testing: \n\n")
				   print(Stest)
				}
				cat("\n")
			}
			cat("----\n\n")
		}	
	}	
	
	if(!is.numeric(pl) || pl < 0L)
	  stop("Print level `pl` must be some positive numeric value.")	
	if(missing(nsim))
	  nsim <- attr(qsd$qldata,"nsim")  	
	if(is.null(nsim) || !is.numeric(nsim))
	  stop("Number of simulations must be given.")
    
	 # may overwrite (observed) statistics	
	 if(!is.null(obs)) {
		  obs <- unlist(obs)
		  if(anyNA(obs) | any(!is.finite(obs)))
			  warning("`NA`, `NaN` or `Inf` values detected in argument `obs`.")
		  if(!is.numeric(obs) || length(obs)!=length(qsd$covT))
			  stop("Object `obs` must be a (named) numeric vector or list of length equal to the number of given statistics in `qsd`.")
		  qsd$obs <- obs
	} 
  
    # wrapping simulator function
    sim <- match.fun(sim)	
	# silently remove not required
	args <- list(...)
	.checkfun(sim,args,remove=TRUE)
	simFun <-
	  structure(
		 function(x) {
			 try(do.call(sim,c(list(x),args)))			
		 },
	     class=c("simObj","function")
	  )
		
    # set default starting point
	x0 <-
	 if(is.null(x0))
	   (qsd$lower + qsd$upper)/2 
     else if(is.list(x0) || is.vector(x0)) {		
		unlist(x0)
	 } else { stop("Starting vector 'x0' must be list or a numeric vector.") }
		
    # check here 
    .checkArguments(qsd,x0,Sigma)	
	
	# clean or invert Sigma if supplied
	if(!is.null(Sigma)){
		if(qsd$var.type == "kriging"){
			Sigma <- NULL
			message("Ignoring `Sigma` because kriging approximation of variance matrix is set.")
		} else if(qsd$var.type == "const") {
			Sigma <- try(gsiInv(Sigma),silent=TRUE)
			if(inherits(Sigma,"try-error"))
				stop("Failed to invert initial estimate `Sigma` as a constant variance matrix.")		
		}
	}
	xdim <- attr(qsd$qldata,"xdim")
	# available local optimization method(s) to choose
	nlopt.fun <- c("cobyla","bobyqa","neldermead",
			       "direct","directL","lbfgs","nloptr")		   
	all.local <- 
	 if(qsd$criterion == "qle") {		
		c("qscoring",nlopt.fun)
	 } else nlopt.fun	
	# quasi-score options
	qscore.opts <- .qsOpts(qscore.opts,xdim)

	loc <- pmatch(method,all.local)
	if(length(loc) == 0L || anyNA(loc)) {
		stop(paste0("Invalid local method(s) for criterion `mahal`. Choose one or more of: ",
						paste(all.local, collapse = ", ")))
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
	
	# set 'globals'
	globals <- .getDefaultGLoptions(xdim)
	if(length(global.opts) > 0L) {
		.checkOptions(globals,global.opts)
		globals[names(global.opts)] <- global.opts				
	}
	
	## set 'locals'
	locals <- .getDefaultLOCoptions(xdim)
	if(length(local.opts) > 0L) {
		.checkOptions(locals,local.opts)
		locals[names(local.opts)] <- local.opts	
	}
	# setting controls as data frame
	ctls <- .setControls(globals,locals)
	# data frame for storing relative estimation error deviation
	# add to `ctls` if testing is enabled
	perr <- 
	 if(locals$test) {
	  if(length(locals$perr_tol)!=xdim)
		locals$perr_tol <- rep(locals$perr_tol,length.out=xdim)  
	  data.frame(cbind("cond" = locals$perr_tol, "val" = rep(Inf,xdim),
			"tmp" = rep(Inf,xdim), "stop" = rep(0,xdim), "count" = rep(globals$NmaxQI,xdim)),
		  row.names = xnames, check.names = FALSE)
	 } else NULL

	# local weights
	if(any(locals$weights > 1L) || any(locals$weights < 0L))
		stop("Weights for local sampling must be
			in the interval [0,1] for sampling criterion `score`!")	
	mWeights <- length(locals$weights)
	# global weights	
	if(any(globals$weights < 0L))
	  stop("Weights for global sampling must be positive!")
  	mWeightsGL <- length(globals$weights)
	# parallel options: seeding
	# the seed is stored if given
	noCluster <- is.null(cl)
	tryCatch({
		if(noCluster){
			type <- if(Sys.info()["sysname"]=="Linux")
						"FORK" else "PSOCK"
			cores <- getOption("mc.cores",1L)
			if(cores > 1L) 
			  try(cl <- parallel::makeCluster(cores,type=type),silent=FALSE)
		    # re-initialize in any case (see `set.seed`)		    
			if(!is.null(cl)){
				clusterSetRNGStream(cl,iseed)			  	 
			} else noCluster <- FALSE
		}				
	},error = function(e)  {
		noCluster <- FALSE
		message(.makeMessage("Could not initialize cluster."))
	})	
	   	
 	# select criterion function	
	criterionFun <- 
		switch(qsd$criterion,
			"mahal" = {				  		  
				  function(x,...) {				  
					mahalDist(x,qsd,Sigma,W=W,theta=theta,
					 cvm=cvm,inverted=TRUE,check=FALSE,...,cl=cl)
				  }  
			 },
			 "qle" = {				  
				 function(x,...)
					quasiDeviance(x,qsd,NULL,W=W,theta=theta,
						cvm=cvm,check=FALSE,...,cl=cl)					
			 }, { stop("Unknown criterion function!") }
		) 
 
	## loop init	
	nglobal <- nlocal <- 0L
	EPS <- .Machine$double.eps^(2/3)
	
	# miniumum distance of sampling candidates
	# to previously sampled points
	eps <- 0.001*prod(abs(qsd$upper-qsd$lower))^(1/xdim)
	# set starting values global/local weights
	maxIter <- globals$maxiter 	# number of global iterations
	maxEval <- globals$maxeval	# amount of evaluation of sample locations
		
	## record iteration status
	status <- list("global"=0L, "minimized"=FALSE) 
			
	# first time CV:
	# this is also for testing 
	# before iterating many times 
	cvm <- NULL
	errId <- pmatch(errType,c("kv","cv","max"))
	if(anyNA(errId) || length(errId)!=1L)
	  stop("Invalid argument `errType`. Please choose one of `kv`, `cv`, `max`")
		   	
    # initialize	
	info <- reset <- TRUE	
	W <- theta <- Stest <- NULL
	QD <- criterionFun(x0)	
	if(.isError(QD)){					
	  return(.qleError(message="Could not compute criterion function.",
		call=match.call(), error=QD))
	}		
	xt <- x <- xold <- x0 									# xt: current, x: starting point, xold: old, x0: initial point
	ft <- f <- fold <- QD[[1]]$value						# see above!
	Snext <- c(QD[[1]],"fval"=f)
	# but then reset so it can be computed again
	if(qsd$var.type != "const")
	 Sigma <- NULL
 
	dummy <- 
	  tryCatch({						
		repeat{		
				# refit
				if(useCV <- (errId > 1)) {
					if(pl > 0L)
					 cat("Update cross-validation covariance models...\n")
					cvm <- try(prefitCV(qsd,type=errType,cl=cl),silent=TRUE) 
					if(.isError(cvm)) {						
						cvm <- NULL
						message("Prefit of CV models failed during final surrogate minimization.")			
					} 
				}					
				# either start a (multistart) local search
				
				if(multistart && status[["global"]] > 1L){
					# always include last sample point `x` as a starting point
					S0 <- multiSearch(x, qsd, method, qscore.opts, control,
							Sigma=Sigma, W=W, theta=theta, inverted=TRUE, cvm=cvm,
							 check=FALSE, pl=0L, nstart=globals$nstart,
							  cl=cl, verbose=pl>0L)
					if(.isError(S0))		
					 message("Could not complete multistart local search.")										
				} else {
				    S0 <- searchMinimizer(x, qsd, method, qscore.opts, control,
					    	Sigma=Sigma, W=W, theta=theta, inverted=TRUE, cvm=cvm,
								check=FALSE, pl=pl, verbose=pl>0L)
				}
				# store local minimization results
				tmplist <- list("S0"=S0)				
				# Set current iterate to last sampled point in case of no convergence
				# during the global phase, eventually a sampled point 'Snext' also becomes
				# a minimizer after a number of steps cycling through the vector of global weights						
				if(!.isError(S0) && S0$convergence >= 0L){					
					 I <- S0$I
					xt <- S0$par
					ft <- S0$value					 
					varS <- S0$varS					
					status[["minimized"]] <- TRUE					
				} else {														
					I  <- Snext$I
					xt <- Snext$par
					ft <- Snext$value
					varS <- Snext$varS					
					status[["minimized"]] <- FALSE					
				}			
																
			    # current iteration did not stop 
			    # choose next type of phase: pure local or global				
				if(status[["minimized"]]){											# found a local minimizer
					if(S0$convergence == 1L){										# found root of quasi-score
						status[["global"]] <- 0L									# start local phase
					} else if(ft < locals$ftol_abs){								# found approximate root
					   if(locals$test){ 
						 if(pl > 0L)
						   cat("Testing local minimizer...\n")
					     Stest <-
							  tryCatch({					    
								  newObs <- simQLdata(simFun, nsim=locals$nobs, X=rbind(xt), cl=cl, verbose=pl>0)
								  if(.isError(newObs))
									  stop(paste(c("Cannot generate data at approximate root: \n\t ",
															  format(xt, digits=6, justify="right")),collapse = " "))				  
								  # test for an approximate root (based on criterion function)
								  .rootTest(xt, ft, I, newObs[[1]], locals$alpha, qsd$criterion,
										  qsd, method, qscore.opts, control, Sigma=Sigma, W=W,
										  theta=theta, cvm=cvm, cl=cl)	
								  
							  }, error = function(e){
								  msg <- .makeMessage("Testing approximate root failed: ",
										  conditionMessage(e))
								  message(msg)
								  .qleError(message=msg,call=match.call(),error=e)
							  }		
						  )
						  # store results in temporary list for failure analysis
						  tmplist <- c(tmplist,list("Stest"=Stest))						  
						  # set status
						  status[["global"]] <-
							  if(.isError(Stest)){
								  msg <- paste(c("Cannot test approximate root at: \n\t ",
												  format(xt, digits=6, justify="right")),collapse = " ")
								  message(msg)
								  2L															# switch to global in case of error
							  } else if(attr(Stest$test,"passed")) { 1L }						# found approximate root  
							    else 2L 														# did not pass the test	
					    } else { status[["global"]] <- 0L } 											
					} else { status[["global"]] <- 2L }
				} else {
					status[["global"]] <- 2L										
					# Though we might have found a local minimizer, we do not
					# sample there because we trust the selection criteria
					# which additionally considers point-interdistances. 
					# We start the next local search from the current sample point
					# and thus imitate some kind of random multistart local search.
				}				
																				
				# find new sample point
				Snext <- 
					tryCatch({						
						if(status[["global"]] < 2L){
							# weighting matrix for variance
							# average approximation and local sampling
							W <- I 	 		
							# and its inverse is used as the variance
							# of theta for sampling from MVN
							theta <- xt	 	
							# found approximate root
						    							
							# stopping conditions for
							# relative estimation error deviation (see qleTest)							
							if(locals$test && !is.null(Stest) && !.isError(Stest)) {
								perr["val"] <- attr(Stest,"relED")	
								if(anyNA(perr["val"]))
								 message("Cannot test stopping conditions while testing local minimizer. `NAs` detected.")
								else if( any(perr["val"] < perr["cond"]) ){									 # can be a vector for each component of the parameter					        
									perr["stop"] <- perr["stop"] + as.integer(perr["val"] < perr["cond"])	 # and count separately for each one
									if(any(perr["stop"] >= globals$NmaxQI))											
										break																			
								} else { perr["stop"] <- rep(0L,xdim) }										 # reset	
							}																	
																	  
							# generate local candidates							
							Y <- nextLOCsample(W,theta,
									locals$nsample,lb=qsd$lower,
									  ub=qsd$upper,pmin=locals$pmin,invert=TRUE)
							
							if(.isError(Y)) {
								status[["global"]] <- 2L
								message(.makeMessage("Sampling local candidates failed. Try global sampling."))								
							} else {							
								# get min distances
		    					dists <- .min.distXY(X,Y)
								# remove too close sample candidates							
								idx <- which(dists < eps)					
								if(length(idx) > 0L){
									Y <- Y[-idx,,drop=FALSE]
									if(nrow(Y) < floor(0.1*length(dists))) {											
										msg <- paste("Number of local candidates is less than 10% of sample size:, "
												,nrow(Y)," try global sampling now.")
										message(msg)
										status[["global"]] <- 2L
									}
									dists <- dists[-idx]
								}
								# sampling might cause switch to global phase
								if(status[["global"]] < 2L){		 
									 nlocal <- nlocal + 1L
									 dmin <- min(dists)
									 dmax <- max(dists)									 
									 id <- 
									   switch(
									    locals$nextSample,
									     "score" = {			
											 # use user defined weights
											   if(locals$useWeights) {										
												  k <- nlocal %% mWeights
												  w <- ifelse(k != 0L,
														  locals$weights[k],
														  locals$weights[mWeights] )
											   } else {
												 
												if(ft < 0.9999*fold || ctls["lam_max","val"] < ctls["lam_max","tmp"])													 
												 {
													 ctls["nfail","val"] <- 0L
													 ctls["nsucc","val"] <- ctls["nsucc","val"] + 1L												
												 } else {
													 ctls["nsucc","val"] <- 0L
													 ctls["nfail","val"] <- ctls["nfail","val"] + 1L												
												 }
												 if(reset){
													 reset <- FALSE
													 w <- locals$weights[1]													 
												 }
												 # update weights									
												 if(ctls["nfail","val"] > 0L && 
												  !(ctls["nfail","val"] %% ctls["nfail","cond"])){											 											 
													 w <- max(w-locals$eta[1],0)
												 } else if(ctls["nsucc","val"] > 0L && 
														 !(ctls["nsucc","val"] %% ctls["nsucc","cond"])){										 
													 w <- min(w+locals$eta[2],1)
												 }									 										 
											  }						
											  # minimize ballanced criterion
										  	  # criterion funtion values at candidates
										  	  # Sigma is re-calculated here at theta)
											  fd <- criterionFun(Y,value.only=2L)
											  if(.isError(fd) || !is.numeric(fd)){
												  stop(paste("Criterion function evaluation failed: ",fd))
											  }
											  smin <- min(fd)
											  smax <- max(fd)
											  sw <- if(abs(smax-smin) < EPS) 1 
												 	 else (fd-smin)/(smax-smin)	
											  dw <- if(abs(dmax-dmin) < EPS) 1		
													 else (dmax-dists)/(dmax-dmin)
											  which.min( w*sw + (1-w)*dw )								
										  },
										  "var" = {						
											  # maximize trace criterion
											  # (same for quasi-deviance and mahalanobis distance)
											  # Sigma is re-calculated here at theta
											  fd <- criterionFun(Y,value.only=3L)								  
											  dw <- if(abs(dmax-dmin) < EPS) 1		
													   else (dists-dmin)/(dmax-dmin)
											  which.max( fd*dw )
										  })									
								} 								
							}							
						} # end local sampling
						
						# start global sampling
						if(status[["global"]] > 1L){							
							reset <- TRUE							
							W <- theta <- NULL   										# no weighting in global phase		
							nglobal <- nglobal + 1L
													
							# sample new candidates
							Y <- sapply(seq_len(ncol(X)),
							  		function(i) {
									  runif(globals$nsample,qsd$lower[i],qsd$upper[i])
									})							
							colnames(Y) <- xnames					
							dists <- .min.distXY(X,Y)									# check for minimum distance between sample points
							idx <- which(dists < eps)
							if(length(idx) > 0L)	{
								Y <- Y[-idx,,drop=FALSE]
								if(nrow(Y) < floor(0.1*length(dists))) {										
							 	   message(.makeMessage("Number of candidates ",nrow(Y)," is not enough to proceed sampling."))									 											
								   break;
								}
								dists <- dists[-idx]
							}							
							# quasi-deviance or Mahalanobis distance as
							# a criterion for global sampling
							fval <- criterionFun(Y,value.only=TRUE)									  																		
							if(.isError(fval) || !is.numeric(fval)){
							  stop(paste("Criterion function evaluation failed: ",fval))
							}
							# next sampling location							
							fmin <- min(fval)
							fmax <- max(fval)
							k <- nglobal %% mWeightsGL
							w <- ifelse(k != 0L, globals$weights[k], globals$weights[mWeightsGL] )
							fd <- if(abs(fmax-fmin) < EPS) 1 
							      	else (fval-fmin)/(fmax-fmin)							
							# next sampling point
							id <- which.max( exp(-w*fd) * dists )
											
						} # end global					
						if(!is.numeric(id) || length(id) == 0L)
						  stop("Could not find index of selection candidate.")
						
					    # compute criterion function at new sample point						
						c( criterionFun(Y[id,])[[1]],"fval"=fd[id] )						
						
					}, error = function(e) {						
						msg <- .makeMessage("Sampling new candidates failed: ",
								  conditionMessage(e))
						message(msg)
						.qleError(message=msg,call=match.call(),error=e)
					}
				)
				
				# criterion at selected point `Snext`
				tracklist <- c(tracklist,
				 list(c(tmplist,"Snext"=list(Snext),"status"=list(status))))				
				
		  		if(.isError(Snext)){					
				  stop(attr(Snext,"error"))
				}								  
			  
				# next sampling location				
				# optional: plot iterates (2D) 
				if(plot && xdim < 3L) {
					p <- rbind(Snext$par,xt)
					if(xdim == 1L) 
					 p <- cbind(p,c(0,0))					 
					cols <- if(status[["global"]]>1L) "blue" else "green"
					try(points(p[1,,drop=FALSE],pch=8,cex=1,col=cols,bg=cols))
					if(status[["global"]] == 1L){
						try(points(p[2,,drop=FALSE],pch=8,cex=1,col="green",bg="green"))	
					} else try(points(p[2,,drop=FALSE],pch=21,cex=0.5,col="magenta",bg="magenta"))
					
				}												
								
				# ---------------- check stopval only for global phase -------------------------
				
				if(status[["global"]] && ft < ctls["stopval","cond"]) 
					ctls["stopval","stop"] <- 1L 													
				
				# Minimum sampling distance reached ?
				dm <- attr(.check.distAll(X,xTol=ctls["sampleTol","cond"]),"min")
				if(dm < ctls["sampleTol","val"])
				  ctls["sampleTol","val"] <- dm
				if(dm < ctls["sampleTol","cond"]) {
					ctls["sampleTol","stop"] <- ctls["sampleTol","stop"] + 1L				
				} else {  ctls["sampleTol","stop"] <- 0L }					  
				
				# -------------------- check xtol and ftol ----------------------------
				
				ctls["xtol_rel","val"] <- max(abs(xt-xold)/pmax(abs(xold),1.0))
				if( ctls["xtol_rel","val"] < ctls["xtol_rel","cond"]) {
					ctls["xtol_rel","stop"] <- ctls["xtol_rel","stop"] + 1L					
				} else { ctls["xtol_rel","stop"] <- 0L }
				
				# ftol_rel (global and local) (low priority)
				ctls["ftol_rel","val"] <- abs(ft-fold)/max(abs(fold),EPS)
				if( ctls["ftol_rel","val"] < ctls["ftol_rel","cond"]) {
					ctls["ftol_rel","stop"] <- ctls["ftol_rel","stop"] + 1L					
				} else { ctls["ftol_rel","stop"] <- 0L }				
				
				# generalized EVD (either based on CV or Kriging variances)
				ctls["lam_max","tmp"] <- ctls["lam_max","val"]
				ctls["lam_max","val"] <- max(geneigen(varS,I,only.values=TRUE), na.rm=TRUE)
				
				# generalized eigenvalues
				if( ctls["lam_max","val"] < ctls["lam_max","cond"]) {
					ctls["lam_max","stop"] <- ctls["lam_max","stop"] + 1L									
				} else { ctls["lam_max","stop"] <- 0L }
				
				# Maximum prediction variance of the quasi-score vector:
				# either CV based or evaluated by kriging variances.				
				if(qsd$krig.type == "var") {
					ctls["C_max","tmp"] <- ctls["C_max","val"]				
					ctls["C_max","val"] <- max(diag(varS))
					test <- abs(ctls["C_max","val"]-ctls["C_max","tmp"])/ctls["C_max","tmp"]
					if(test < ctls["C_max","cond"]) {
						ctls["C_max","stop"] <- ctls["C_max","stop"] + 1L
						if(ctls["C_max","stop"] >= globals$NmaxCV) {
							ctls["C_max","val"] <- test							
						}
					} else { ctls["C_max","stop"] <- 0L }
				}
				
				# show info
				.printInfo()
				# print stopping conditions				
				.showConditions()
								
				# update current iterate
				if(status[["global"]] > 1L){
					xold <- xt
					fold <- ft	
					x <- Snext$par
					f <- Snext$value					
				} else {
					x <- xold <- xt
					f <- fold <- ft											
				}
			
				# ----------------- maximum iterations/evaluations -----------------------------
				
                # If stopping at global phase, then the current sample point
				# at maximum weight corresponds to a sampled minimum of the
                # criterion function if not locally minimized.
				
				if(nglobal >= maxIter){
					if(status[["global"]] > 1L){
						if(w == max(globals$weights)){  		# stop if minimum of criterion function is sampled
							ctls["maxiter","stop"] <- 1L
						    xt <- x; ft <- f					# set to globally sampled point
							status[["minimized"]] <- FALSE
						}
					} else ctls["maxiter","stop"] <- 1L 
				} else if((nglobal+nlocal) >= maxEval){
					if(status[["global"]] > 1L){
						if(w == max(globals$weights)) { 		# stop if minimum of criterion function is sampled
							ctls["maxeval","stop"] <- 1L
							xt <- x; ft <- f					# set to globally sampled point
							status[["minimized"]] <- FALSE
						}
					} else ctls["maxeval","stop"] <- 1L 
				}
				# stop main loop
				if(any(ctls[,"stop"] >= ctls[,"count"])) break;
				
				# simulate at new locations
				# new simulations, qsd$nsim is default
				newSim <-
					tryCatch(
						simQLdata(simFun, nsim=nsim, X=Snext$par, cl=cl, verbose=pl>0L),
						error = function(e) {
							msg <- .makeMessage("Simulating the model failed: ",conditionMessage(e))
					 		.qleError(message=msg,call=match.call(),error=e)
						})		 
				if(.isError(newSim)){
					tracklist[[length(tracklist)]]$newSim <- newSim
					msg <- paste(c("Cannot simulate data at candidate point: \n\t ",
					  format(Snext$par, digits=6, justify="right")),collapse = " ")
					stop(msg)					
				}

				# update QL
				qsd <-
				 if(status[["minimized"]] && locals$test &&
					!is.null(Stest) && !.isError(newObs)) {				
					  # `d`= sample point, `x` is a local minimum
					  updateQLmodel(qsd, rbind("d"=Snext$par,"x"=xt), 
							 structure(c(newSim,newObs),nsim=c(nsim,locals$nobs),class="simQL"),						 
							 fit=TRUE, cl=cl, verbose=pl>0L)					 
				 } else {
					  updateQLmodel(qsd, Snext$par, newSim,						 
							 fit=TRUE, cl=cl, verbose=pl>0L)
				 }
									
				# check results of kriging
				if(.isError(qsd)){
					msg <- paste(c("Cannot update QL model at candidate point: \n\t ",
						     format(as.numeric(Snext$par), digits=6, justify="right")),collapse = " ")				 			    
					e <- attr(qsd,"error")
					if(inherits(e,"error"))
					  msg <- c(msg, conditionMessage(e))
					stop(msg)
				}					
				Stest <- NULL
				# update X sample			
				X <- as.matrix(qsd$qldata[seq(xdim)])				
			}																				# end main loop
		}, error = function(e) {
			msg <- .makeMessage("Current iteration stopped unexpectedly: ",
					  conditionMessage(e))
			message(msg)
			structure(
				list("par"=xt,
					 "value"=ft,
					 "ctls"=rbind(ctls,perr),					 
					 "qsd"=qsd,
					 "cvm"=cvm,
					 "why"=NULL,					 
					 "final" = S0,															# local results
					 "convergence"=FALSE),				
				tracklist = tracklist,				
				optInfo = list("x0"=x0,
							   "W"=W,
							   "theta"=theta,						   
							   "last.global"=(status[["global"]] == 2L),
							   "minimized"=status[["minimized"]],
							   "useCV"=useCV,
							   "method"=S0$method,
							   "nsim"=nsim,
							   "iseed"=iseed),
				class = c("qle","error"), call = sys.call(), error=e)	
							
		}, finally = {
		  if(noCluster) {
			if(inherits(try(stopCluster(cl),silent=TRUE),"try-error")){
			  	rm(cl)				
				message("Error in stopping cluster.")
		  	} else {
				cl <- NULL				
			}
			invisible(gc())
		  }
		}
	) # end outer tryCatch	
	
	# stop on error 
	if(.isError(dummy))
	 return(dummy)
	
 	# Last iteration was done at global phase, so try to minimize again
	# either from the last sample point since it this most local supposed to
	# small for high (global) weights or by a multistart approach.
	
 	if(status[["global"]] == 2L){
		if(multistart){
			# always include last sample point `x` as a starting point
			S0 <- multiSearch(x, qsd, method, qscore.opts, control,
					Sigma=Sigma, W=W, theta=theta, inverted=TRUE, cvm=cvm,
					 check=FALSE, pl=0L, nstart=globals$nstart,
					cl=cl, verbose=pl>0L)
			if(.isError(S0))		
				message("Could not complete multistart local search.")										
		} else {
			S0 <- searchMinimizer(x, qsd, method, qscore.opts, control,
					Sigma=Sigma, W=W, theta=theta, inverted=TRUE, cvm=cvm,
					 check=FALSE, pl=pl, verbose=pl>0L)
		}
	   	# overwrite last sample point if local minimization was successful
		if(!.isError(S0) && S0$convergence >= 0L){					
			xt <- S0$par
			ft <- S0$value					 
			status[["minimized"]] <- TRUE					
		}
		# store local minimization results as these are overwritten now
		# where `Snext` does not change anymore, however, add it
		tracklist <- c(tracklist,
		 list("S0"=S0,"Snext"=Snext,"status"=list(status)))
	    
        # show results (update) 
		.printInfo()
		.showConditions()	
	}
 
	## only for estimte theta=(xt,ft)	
	ctls["stopval",c(2,4)] <- c(ft,ft < ctls["stopval","cond"])	
	ctls["maxiter",c(2,4)] <- c(nglobal,nglobal >= maxIter)
	ctls["maxeval",c(2,4)] <- c(nglobal+nlocal, nglobal+nlocal >= maxEval)	

	# remove `nfail`, `nsucc`
	ctls <-
	 if(.isError(S0)) {
	   message("Last search results have errors. Please see list element `\"final\")`.")	   
	   ctls[1:8,-3]
	 } else {	  	
		val <- max(abs(S0$score))
		ctls <- rbind(ctls[1:8,-3],   												# remove `tmp` column
		    	as.data.frame(cbind("cond" = qscore.opts$score_tol,
			 				        "val" = val ,
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
			 "final" = S0,			 												# local results (NULL if `last.global` is TRUE)
			 "convergence"=(status[["minimized"]] && length(arg.names) > 0L)),	 	
		tracklist = tracklist,		
		optInfo = list("x0"=x0,
				 	   "W"=W,
					   "theta"=theta,						   
			    	   "last.global"=(status[["global"]]==2L),
				  	   "minimized"=status[["minimized"]],
				  	   "useCV"=useCV,
				  	   "method"=S0$method,
				  	   "nsim"=nsim,
					   "iseed"=iseed),
		 class = "qle", call = sys.call())  	
 	 
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
#' @param ... 	 ignored, additional arguments
#' 
#' @rdname print.qle
#' @method print qle
#' @export 
print.qle <- function(x, pl = 2, digits = 5,...){	
	if(.isError(x)) 
	  stop(.makeMessage("Estimation result has errors.","\n"))	
	if(!inherits(x,"qle"))
	  stop("Method is only for objects of class `qle`.")
    if(!is.numeric(pl) || pl < 0L )
	  stop("Print level must be a positive numeric value.")
  
	if(pl > 0L) {
	  	if(x$qsd$criterion == "mahal")
		 cat("Mahalanobis distance: \n\n",x$value,"\n\n")
		else			
		 cat("Quasi-deviance: \n\n",x$value,"\n\n")
	  	
		cat("Estimate:\n\n")
		print.default(formatC(signif(x$par, digits = digits), digits = digits, format="fg", flag="#"),
				print.gap = 4, quote = FALSE)
		
		cat("\n")
		if(pl > 1L) {			
			if(x$convergence >= 0L) {
				by <- x$ctls[x$why,"val"]
				names(by) <- x$why
				cat("Convergence: ")
				if("maxeval" %in% x$why)
				  cat("maximum evaluations reached.\n\n")
			    else cat("\n\n")
				print.default(formatC(by, format="e", digits=digits), right=FALSE, print.gap=2,
						quote=FALSE)
			} else cat("(None of convergence criteria matched.)\n\n")
		} 
	}
  	if(pl > 1L) {
		cat("\n")	
		nsim <- unlist(x$ctls["maxeval","val"])*attr(x,"optInfo")$nsim	
		cat("Evaluations: ",unlist(x$ctls["maxeval","val"])," (simulations: ",nsim,")\n\n")
		cat("Variance approximation type: ",x$qsd$var.type,"\n")	
	}
	if(pl > 2L) {
		if(!.isError(x$final)) {
			cat("\n\n ***  Final results *** \n\n\n")			
			print(x$final)
	 	}
		if(x$qsd$var.type != "kriging"){
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
#' @param ... 	    ignored, additional arguments
#' 
#' @rdname print.QSResult
#' @method print QSResult
#' @export 
print.QSResult <- function(x, pl = 1, digits = 5,...) {	
	if(.isError(x)) 
	  stop(.makeMessage("Quasi-scoring iteration had errors.","\n"))
	if(!inherits(x,"QSResult"))
	  stop("Method is only for objects of class `QSResult`.")	
    if(!is.numeric(pl) || pl < 0L )
	 stop("Print level must be a positive numeric value.")
 
	cat(paste0("Local method:\n\n `",x$method,"`\n\n"))		
	cat("Solution: \n\n")
	print.default(formatC(signif(x$par, digits = digits), digits = digits, format="fg", flag="#"),
			print.gap = 4, quote = FALSE)
	cat("\n")
	cat("Quasi-deviance:\n\n",x$value,"\n\n")
	cat("Iterations....",x$iter,"\n")		
	cat("Convergence...",x$convergence >= 0L,"\n")				
	if(!is.null(x$score)){
		cat("\nStopped by: \n\n",unlist(strsplit(paste(x$message, "\n" ),':')), fill=TRUE )
		if(x$criterion == "qle") {
			cat("Quasi-score:\n\n")
		}  else cat("Gradient:\n\n")
		print.default(formatC(signif(x$score, digits=8), digits=8, format="fg", flag="#"),
				print.gap = 4, quote = FALSE)				
	}
	cat("\n")
	
	if(pl > 1L) {
		Sigma <- attr(x,"Sigma")		
		if(!is.null(Sigma)) {	
			cat("\nApproximation of variance matrix: \n\n Sigma = \n\n")
			print(Sigma)
			cat("\n\n")
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
#' @description Generate a random sample of points as a set of candidates for evaluation 
#'
#' @param S			    variance matrix of sample points (usually chosen as the information matrix)
#' @param x			    an approximate root as the mean value if the MVN distribution
#' @param n				number of points to sample
#' @param lb			vector of lower bounds of the hypercube
#' @param ub			vector of upper bounds of the hypercube
#' @param pmin			minimum required probability to cover the hypercube (parameter space)
#' @param invert	    optional, \code{invert=FALSE} (default) for no inversion of `\code{S}`
#'
#' @return Matrix of sampled locations.
#'
#' @details The function generates a random sample of points with mean and variance given by `\code{x}`
#' 	 and `\code{S}`, respectively, according to a (truncated) multivariate normal distribution (using
#'   rejection sampling) to match the parameter space given by the lower and upper bound vectors. 	 
#'
#' @examples
#'  X <- nextLOCsample(matrix(c(1,0,0,1),nr=2), c(0,0), 10, c(-0.5,-0.5), c(0.5,0.5))
#'  
#' @author M. Baaske
#' @rdname nextLOCsample
#' 
#' @importFrom mvtnorm pmvnorm rmvnorm
#' @export
nextLOCsample <- function(S, x, n, lb, ub, pmin = 0.05, invert = FALSE) {
	if(invert) 
	  S <- try(gsiInv(S),silent=TRUE)
	if(.isError(S)) {
		msg <- paste0("Failed to invert matrix.")
		message(msg)
		return(.qleError(message=msg,call=match.call(),error=S))		
	}
	if(anyNA(S) || !is.matrix(S))
		warning("Variance matrix has `NaN` values. A matrix?")
	if( rcond(S) < sqrt(.Machine$double.eps)) {
	  warning(" Variance matrix is nearly ill conditioned.")
	  if( (is.pos = .isPosDef(X)) != 0L )
		 return(.qleError(
				  message=paste0("Variance matrix not positive definite: ",is.pos),
				  call=match.call(),
				  S=structure(S,is.pos=is.pos))) 
 	}
	
	# Akzeptanzrate p aus der Multivariaten Normalverteilung bestimmen	
	p <- mvtnorm::pmvnorm(lower=lb, upper=ub, mean=x, sigma=S)
	if (p < pmin) {
	   msg <- paste0(sprintf("Sampling local candidates is too inefficient (p=%.6f).", p))
	   message(msg)
	   return( .qleError(message=msg,call=match.call(),error=p) )
 	}
	# sampling
	X <- try(.rejectSampling(n, p, x, S, lb, ub),silent=TRUE)
	if(inherits(X,"try-error") || !attr(X,"success")) {
		msg <- paste0("Rejection sampling failed. ")
		message(msg)
		return(.qleError(message=msg,call=match.call(),error=X))
	}	
	return ( structure(X, p = p) )			
}

#' @name qscoring
#' 
#' @title Quasi-scoring iteration
#'
#' @description The function solves the quasi-score equation by a root finding algorithm similar to Fisher's scoring iteration.
#'
#' @param qsd    	object of class \code{\link{QLmodel}}
#' @param x0		(named) numeric vector, the starting parameter
#' @param opts		quasi-scoring options, see details
#' @param Sigma	    a pre-specified variance matrix estimate
#' @param ...	    further arguments passed to function \code{\link{covarTx}}
#' @param inverted  currently ignored
#' @param check		logical, \code{TRUE} (default), whether to check input arguments
#' @param cvm		list of covariance models for cross-validation (see \code{\link{prefitCV}})
#' @param Iobs	    logical, \code{FALSE} default, whether to compute observed quasi-information matrix
#' @param pl	    numeric, print level, use \code{pl}>0 to give intermediate output  	
#' @param verbose   \code{FALSE} (default), otherwise print intermediate output
#'
#' @return List of results of quasi-scoring iteration.
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
#'  constraints (see, e.g. [1,3]) specifically tailored for quasi-likelihood parameter estimation. Due to the typical
#'  nonconvex nature of the (unknown and not further specified) quasi-likelihood function as an objective
#'  function one needs some kind of globalization strategy in order to stabilize the descent step and to avoid a premature
#'  termination. Therfore, we use the quasi-deviance function as a monitor function (see vignette) though it does not
#'  inherit all of the appreciable properties of a true objective function such as among others, for example,
#'  identifying appropriate descent directions. However, these are general numerical obsticles in case of pure root
#'  finding algorithms and need to be addressed elsewhere. 
#'  
#'  \subsection{Quasi-scoring under uncertainty}{ 
#'  The quasi-scoring iteration covers both kinds of prediction variances, kriging-based and those by a CV approach, which account for
#'  the uncertainty induced by the quasi-score approximation model. By default kriging variances
#'  are included in the computation during all iterations. If fitted covariance models `\code{cvm}` are supplied by the user
#'  in advance (see \code{\link{prefitCV}}), the variances of prediction errors of each statistic are separately evaluated by the proposed CV
#'  approach for each new point. For the price of relatively high computational costs those prediction variances
#'  are intended to increase the robustness against false roots due to simulation and approximation errors of the quasi-score function.
#' 
#'  Opposed to this, the user also has the option to carry out a "pure version" of quasi-scoring without accounting for
#'  these errors. This can be set earlier as an option in \code{\link{QLmodel}}. See also \code{\link{covarTx}} and
#'  \code{\link{mahalDist}} for details on how to choose the variance matrix approximation of the statistics.
#'  }
#' 
#'  The following algorithmic options, which can be set by `\code{opts}`, are available:
#'  \itemize{
#'   	\item{\code{ftol_stop}:}{ minimum value of the quasi-deviance for stopping the scoring iteration}
#' 		\item{\code{ftol_abs}:}{ minimum value of the quasi-deviance which is used as a reference value for a local minimizer}
#' 		\item{\code{xtol_rel}, \code{ftol_rel}:}{ see \code{\link{qle}} }
#' 		\item{\code{grad_tol}:}{ upper bound on the quasi-score vector components,
#' 				 testing for a local minimum of the quasi-deviance in case of a line search failure}
#' 		\item{\code{score_tol}:}{ upper bound on the quasi-score vector components, testing for an approximate root}
#'      \item{\code{maxiter}:}{ maximum allowed number of iterations}
#' 		\item{\code{xscale}:}{ numeric, default is vector of 1, typical magnitudes of vector components of `\code{x0}`, e.g. order of upper bounds of the parameter space}
#'      \item{\code{fscale}:}{ numeric, default is vector of 1, typical magnitudes of quasi-score components}
#' 	    \item{\code{pl}:}{ print level (>=0), use \code{pl}=10 to print individual
#' 							 iterates and further values}
#'  } 
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
qscoring <- function(qsd, x0, opts = list(), Sigma = NULL, ...,
			    	  inverted = FALSE, check = TRUE, cvm = NULL, 
				 	   Iobs = TRUE, pl = 0L, verbose = FALSE)
{	
	if(check)
	 .checkArguments(qsd,x0,Sigma)
 	stopifnot(is.numeric(pl) && pl >= 0L )
 
  	if(qsd$criterion!="qle")
	  stop("Quasi-scoring is only valid for criterion `qle`.")
  
    xdim <- attr(qsd$qldata,"xdim")
  	X <- as.matrix(qsd$qldata[seq(xdim)])  		
	if(qsd$var.type != "kriging" && is.null(Sigma)){
		# Only mean covariance matrix is estimated here. 
		# Adding prediction variances (kriging/CV) at C level		
		if(qsd$var.type %in% c("wcholMean","wlogMean")){
			nms <- names(list(...))
			if(!all( c("W","theta") %in% nms))
			 message(paste0("Found `var.type`=\"",qsd$var.type, "\" but no weighting matrix `W` or estimate `theta` was supplied!."))		
		}
		Sigma <- covarTx(qsd,...,cvm=cvm)[[1]]$VTX		
	} 
	
	# set quasi-scoring options
	opts <- .qsOpts(opts,xdim,pl)
	qlopts <- list("varType"=qsd$var.type,					  
				   "useCV"=!is.null(cvm),
				   "useSigma"=FALSE,
		   		   "Iobs"=Iobs)
		   
	try(.Call(C_QSopt, x0, qsd, qlopts, X, Sigma, cvm, opts), silent=TRUE)	
}


## intern only
## Conduct next simulations,
## and update covariance models
updateQLmodel <- function(qsd, Xnew, newSim, fit = TRUE, cl = NULL, verbose = FALSE ){	
    if(verbose)
	 cat("Setup next data...\n")
    stopifnot(nrow(Xnew)==length(newSim))
	
	nextData <-
		tryCatch(
		  setQLdata(newSim,
				    Xnew,
					qsd$var.type,
				    attr(qsd$qldata,"Nb"),
					verbose=verbose),
			error = function(e) {
				msg <- .makeMessage("Construction of quasi-likelihood data failed: ",
						   conditionMessage(e))				
				.qleError(message=msg,call=match.call(),error=e)
			}
	)
    if(.isError(nextData))
	  return(nextData)
	# combine to new data and update
	if(verbose)
	  cat("Update covariance models...\n")	
    updateCovModels(qsd,nextData,fit,cl=cl,verbose=verbose)
}
