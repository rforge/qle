# Copyright (C) 2017 Markus Baaske. All Rights Reserved.
# This code is published under the L-GPL.
#
# File: 	qsOpt.R
# Date:  	12.04.2017
# Author: 	Markus Baaske
#
# Quasi-likelihood estimation

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
	( class(x)=="error" ||
	  !is.null(attr(x,"error")) ||
	  inherits(x, "error") ||
	  inherits(x, "try-error") 
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

.qsOpts <- function(options = list(),pl = 0L) {
	opts <- .addQscoreOptions()
	opts$pl <- pl
	if(length(options) > 0L) {
	 .checkOptions(opts,options)
	 namc <- match.arg(names(options), choices = names(opts), several.ok = TRUE)
	 if (!is.null(namc))
	 opts[namc] <- options[namc]
	}	
	return(opts)
}

.addQscoreOptions <- function() {
	list( "ftol_stop" = .Machine$double.eps^(1/2),
		  "xtol_rel"  = .Machine$double.eps^(1/2),
		  "grad_tol"  = .Machine$double.eps^0.25,
		  "ftol_rel"  = .Machine$double.eps^(1/3),
		  "score_tol" = 1e-3,
		  "slope_tol" = 1e-2,
		  "maxiter" = 100,
		  "pl"=0L)
}

.getDefaultGLoptions <- function(xdim) {
	list("stopval" = 0,								 		# global stopping value
		 "C_max" = 1e-3,
		 "xtol_rel" = .Machine$double.eps^(1/3),
		 "maxiter" = 50,									# max number of global iterations
		 "maxeval" = 25,									# max number of global and local iterations
		 "sampleTol" = .Machine$double.eps^0.25,			# minimum (euclidean) distance between samples		 
	 	 "weights"=c(0.01,0.1,2,4,6,8,10,20,50),		 
		 "nsample" = (xdim+1)*2000,							# number of global random samples
		 "NmaxRel" = 5,
		 "NmaxCV" = 3,		 
		 "NmaxSample" = 3,
		 "NmaxLam" = 3,
		 "Nmaxftol"= 3)
}

.getDefaultLOCoptions <- function(xdim) {
	list("ftol_rel"  = .Machine$double.eps^(1/3),		 	
		 "lam_max" = 1e-2,
		 "pmin" = 0.05,
		 "weights"=c(0.005,0.1,0.2,0.4,0.6,0.8,0.995),	   # only for sampling with criterion `score`
		 "nsample" = (xdim+1)*1000,						   # number of local random samples
		 "ftol_abs" = 0.01,								   # local minimum value tolerance		 
		 "eta" =c(0.05,0.075),							   # c("decrease"=0.05,"increase"=0.075) additive step size
		 "nfail" = 3,
		 "nsucc" = 2,
		 "nextSample"="score",
		 "useWeights"=FALSE)							   # dynamically adjust weights									   
}

.setControls <- function(globals,locals) {
	defaults <- c("C_max","lam_max","xtol_rel","stopval","sampleTol",
				  "nfail","nsucc","ftol_rel","maxiter","maxeval")
	optlist <- c(globals,locals,"score_tol")
	namc <- match.arg(names(optlist), choices = defaults, several.ok = TRUE)
	ctls <- data.frame(cbind("cond" = unlist(optlist[namc]), "val" = 0, "tmp"=0, "stop" = 0),
			row.names = namc, check.names = FALSE)
	
	# init some controls	
	ctls["sampleTol","val"] <- 1E100
	ctls[c("C_max","lam_max"),"val"] <- rep(1,2)		
	return (ctls)
}

#' @name getDefaultOptions
#' 
#' @title Show default options
#' 
#' @description Print all default options for optimization
#' 
#' @param xdim 		size of problem dimension (number of model parameters to estimate)
#' 
#' @return List of options.
#' 
#' @details The function returns a lists of all available options
#'  for functions \code{\link{qscoring}} and \code{\link{qle}}.
#' 
#' @examples
#' \dontrun{
#' # options for a 3D parameter estimation
#' getDefaultOptions(3)
#' }
#' 
#' @author M. Baaske
#' @rdname getDefaultOptions
#' @export
getDefaultOptions <- function(xdim) {
	if(!is.numeric(xdim))
	  stop("`xdim` mus be a numeric value.")
  
	list("qscoring" = .addQscoreOptions(),
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
#' @description Estimate prediction variances based on cross-validation for
#'  each sample means of statistics 
#'
#' @param qsd   	object of class \code{\link{QLmodel}}
#' @param cvm		list of prefitted covariance models (see \code{\link{prefitCV}})
#' @param theta		optional, default \code{NULL}, list or matrix of points where to estimate prediction variances
#' @param type		name of prediction variance measure 
#' @param cl	    cluster object, \code{NULL} (default), see \code{\link[parallel]{makeCluster}} 
#'
#' @return Matrix of estimated prediction variances for each point given by \code{theta} (as a row)
#'  and statistics (as columns).  
#'
#' @details	Other than the kriging prediction variance, which solely depends on interdistances of sample points
#'  and estimated covariance parameters of some assumed to be known spatial covariance structure, the cross-validation
#'  based approach even takes into account the predicted values at `\code{theta}` and thus can be seen as a more robust
#'  measure of variability between different spatial locations. By default, `\code{theta}` equals the current sampling set 
#'  stored in the object `\code{qsd}`.
#' 
#'  If we set the error `\code{type}` equal to "\code{cve}", the impact on the level of accuracy (predicting at unsampled
#'  points) is measured by the \emph{delete-k jackknifed variance} of prediction errors, which then are calculated by a
#'  particular cross-validation approach (see [4] and the vignette). This approach does not require further simulations as a
#'  measure of uncertainty for predicting the sample means of statistics at new candidate points accross the parameter space.
#'  Note that if attribute \code{attr(cvm,"type")} equals "\code{max}", then the maximum of kriging and CV-based prediction
#'  variances is returned (see \code{\link{prefitCV}}) for option `\code{type}`=="\code{cve}". 
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
#'  covariance model (see \code{\link{prefitCV}}) using the option `\code{reduce}` equal to \code{FALSE} based on exactly
#'  one left out point. Also, we can calculate the kriging variance at the left-out sample points if we set the option `\code{type}`
#'  equal to "\code{sigK}". 
#'
#' @examples
#' \dontrun{
#' data(normal)
#' 
#' # design matrix and statistics
#' X <- as.matrix(qsd$qldata[,1:2])
#' Tstat <- qsd$qldata[grep("^mean.",names(qsd$qldata))]
#' 
#' # fit CV models
#' cvm <- prefitCV(qsd, fit = TRUE)
#' theta0 <- c("mu"=2,"sd"=1)
#' 
#' # get RMSD by cross-validation 
#' crossValTx(qsd, cvm, theta0, type = "rmsd")
#' 
#' # and kriging standard error  
#' sqrt(varKM(qsd$covT,theta0,X,Tstat)) 	 
#' }
#' 
#' @seealso \code{\link[parallel]{makeCluster}}
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
				xm$start <- xm$param[xm$fixed]				
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
#' @description The function constructs a list of covariance models
#'  of statistics in order to estimate the prediction error
#'  variances by a cross-validation (CV) approach at unsampled points. 
#'
#' @param qsd   	  object of class \code{\link{QLmodel}}
#' @param fit		  if \code{TRUE} (default), re-estimate covariance parameters
#' @param reduce	  if \code{TRUE} (default), reduce the number of covariance models to refit
#' @param type		  type of prediction variances, "\code{cv}" (default), see \code{\link{qle}}
#' @param control	  control arguments for REML estimation passed to \code{\link[nloptr]{nloptr}}  	
#' @param cl	      cluster object, \code{NULL} (default), see \code{\link[parallel]{makeCluster}}
#' @param verbose	  if \code{TRUE}, print intermediate output
#'
#' @return 
#' 	A list of certain length depending on the current sample size (number of evaluated points).
#'  Each list element corresponds to a (reduced) number of sample points with at most \eqn{k} points
#'  (see details) left out for fitting the covariance models. 
#'
#' @details
#'  Using the CV-based approach (see vignette) for estimating the prediction variances 
#' 	requires a refit of covariance parameters of each statistic based on leaving out a certain number of sample points.
#'  The covariance models are refitted if `\code{fit}` equals \code{TRUE} and otherwise simply updated without fitting
#'  saving some computational resources. The number of points left out is dynamically adjusted depending on the number
#'  of sample points in order to prevent the main estimation algorithm to fit as many models as there are points already
#'  sampled, which would surely be prohibitive in most practical situations. 
#' 
#'  For CV the number \eqn{n_c} of covariance models still to fit, that is, the number of
#'  partitioning sets of sample points, is limited by \eqn{n_c\leq n}, with maximum \eqn{k} sampling
#'  points deleted from the full sample set with overall \eqn{n} sample points such that \eqn{n=n_c k} (see vignette). 
#' 
#' @examples 
#'  \dontrun{
#'   data(normal)
#'   cvm <- prefitCV(qsd)
#'  }
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
	
	fit <- isTRUE(qsd$cvfit)
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

.qdFunx <- function(x) { .Call(C_qDValue,x) } 

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
			 warning("Found `var.type`=\"",qsd$var.type, "\" but no weighting matrix `W` or estimate `theta` was supplied!.")		
		}
		Sigma <- covarTx(qsd,...,cvm=cvm)[[1]]$VTX	
	} else if(useSigma && !inverted){
		# Only for constant Sigma, which is used as is!
		Sigma <- try(solve(Sigma),silent=TRUE)
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
	if(class(qsd) != "QLmodel")
	  stop("`qsd` object must be of class `QLmodel`.")
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
			stop("`Sigma` must be `Null` if using kriging approximation of variance matrix.")	    
		} else if(qsd$var.type == "const" && qsd$criterion == "qle")
			stop("`Sigma` cannot be used as a constant variance matrix for criterion `qle`.")			
		  else if(!isTRUE(attr(qsd$qldata,"chol"))) {
			stop("Cholesky decomposed terms must be given. See function `setQLdata`.")
		}		 	    
		
	} else if(qsd$var.type == "kriging" && is.null(qsd$covL))
		stop("Covariance models for kriging variance matrix must be given, see function `setQLdata`.")	
	  else if(qsd$var.type == "const") 
		stop("`Sigma` must not be NULL for `const` variance matrix approximation.")
	  else if(!isTRUE(attr(qsd$qldata,"chol")))
		stop("Cholesky decomposed terms must be given. See function `setQLdata`.")		
}


##
## TODO: 
##	SQP algorithm
##	(linear inequality constraints)
##
## 	Algorithms: 
##		With bounds and ineq. constr.: cobyla, slsqp (as qscoring variant)
##		Only bounds: direct, bobyqa, neldermead, lbfgs (and qscoring)
##
##	Comment: Sampling candidates must be adapted,
##		(rejection sampling: compute lb and ub at theta! )
##
#' @name searchMinimizer
#'
#' @title Minimize a criterion function 
#'
#' @description The function finds a root of the
#'   quasi-score vector or minimizes one of the criterion functions.
#'
#' @param x0		  (named) numeric vector, the starting point
#' @param qsd   	  object of class \code{\link{QLmodel}}
#' @param method	  names of minimization routines to be used (see details) 
#' @param opts		  list of control arguments for quasi-scoring passed to \code{\link{qscoring}}
#' @param control 	  list of control arguments passed to the auxiliary routines
#' @param ...		  further arguments passed to \code{\link{covarTx}}
#' @param obs		  numeric vector of observed statistics, overwrites `\code{qsd$obs}`
#' @param info		  additional information at found minimizer
#' @param check		  logical, \code{TRUE} (default), whether to check input arguments
#' @param pl		  numeric value (>=0), the print level 
#' @param verbose	  if \code{TRUE} (default), print intermediate output
#'
#' @details The function provides an interface to local and global numerical minimization
#'  using the approximate quasi-deviance (QD) or Mahalanobis distance (MD) as a criterion.
#'  
#'  The function does not require additional simulations to find an approximate minimizer. The
#'  numerical iterations always take place on the fast to evaluate criterion function approximations.
#'  The main purpose is to provide an entry point for minimization without the
#'  need of sampling new candidate points for evaluation. This is particularly useful if we search
#'  for a "first-shot" minimizer. 
#' 
#'  The criterion function is treated as a deterministic (non-random) function during minimization
#'  (or root finding) and depends on the sample points. Because of the typical nonconvex nature of the
#'  criterion functions one cannot expect a global minimizer by applying any local search method like,
#'  for example, the scoring iteration \code{\link{qscoring}}.
#'  Therfore, if the scoring iteration or some other available method gets stuck in a possibly local
#'  minimum of the criterion function showing at least some kind of numerical convergence we use such
#'  minimizer as it is and finish the search possibly being unlucky having not found an appropriate root
#'  of the quasi-score vector (or global minimum of the Mahalanobis distance). If there is no convergence
#'  the search is restarted by switching to the next user supplied minimization routine defined by
#'  `\code{method}`. 
#' 
#'  \subsection{Choice of auxiliary minimization methods}{  
#'  Besides the local quasi-scoring (QS) iteration, `\code{method}` equal to "\code{qscoring}", the following
#'  (derivative-free) auxiliary methods from the \code{\link[nloptr]{nloptr}} package are available for minimizing
#'  both criterion functions:
#'  
#' 	\itemize{
#' 	  \item{}{ \code{\link[nloptr]{bobyqa}}, \code{\link[nloptr]{cobyla}} and \code{\link[nloptr]{neldermead}}}
#'    \item{}{ \code{\link[nloptr]{direct}}, locally biased version: \code{directL}}
#' 	  \item{}{ \code{\link[nloptr]{lbfgs}} for minimizing the MD with constant `\code{Sigma}` only}
#' 	  \item{}{ \code{\link[nloptr]{nloptr}} as the general optimizer, which allows to use further methods}
#'  }
#'    
#'  Using quasi-scoring first, which is only valid for minimizing the QD function, is always a good idea since we might have done
#'  a good guess being in the vicinity of an approximate root. If it fails we switch to any of the above alternative methods
#'  (e.g. \code{\link[nloptr]{bobyqa}} as the default method) or eventually - in some real hard situations - to the
#'  method `\code{direct}` or its locally biased version `\code{directL}`. The order of processing is determined
#'  by the order of appearance of the names in `\code{method}`. Any method available from package `\code{nloptr}` can be
#'  chosen. In particular, setting \code{method="nloptr"} and `\code{control}` allows to choose a multistart algorithm such
#'  as \code{\link[nloptr]{mlsl}}. Only if there are reasonable arguments against quasi-scoring, such as expecting a local
#'  minimum rather than a root first or an available limited computational budget, we can always apply
#'  the direct search method `\code{direct}` leading to a globally exhaustive search. Note that we must always supply a starting
#'  point `\code{x0}`, which could be any vector valued parameter of the parameter space unless method `\code{direct}` is
#'  chosen. Then `\code{x0}` is still required but ignored as a starting point since it uses the "center point" of
#'  the (hyper)box constraints internally. In addition, if CV models `\code{cvm}` are given, the CV based prediction variances
#'  are inherently used during consecutive iterations of all methods. This results in additional computational efforts
#'  due to the repeated evaluations of the statistics to calculate these variances for each new iterate.  
#' }
#' 
#' @return A list as follows
#' 	  \item{par}{solution vector}
#' 	  \item{value}{objective value}
#' 	  \item{method}{applied method}
#' 	  \item{status}{termination code}
#' 	  \item{score}{if applicable, quasi-score vector (or gradient of MD)}
#' 	  \item{convergence}{logical, indicates numerical convergence}	 
#' 
#' @examples
#' data(normal)
#' # many iterations needed with `cobyla` 
#' searchMinimizer(c("mu"=2.5,"sd"=0.2),qsd,method="cobyla",verbose=TRUE)
#' # if quasi-scoring fails, use DFO method 
#' searchMinimizer(c("mu"=2.5,"sd"=0.2),qsd,
#'  method=c("qscoring","bobyqa"),verbose=TRUE) 
#' 
#' @seealso \code{\link[nloptr]{nloptr}}
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
	
 	fun.name <- ""
	x0 <- unlist(x0)
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
		if(!is.na(m1)) {
		 if(m1!=1)
		  method <- c("qscoring",method[-m1])		
		 tryCatch({			
		    qscoring(qsd,x0,opts,...,check=FALSE,pl=pl,verbose=verbose)
		   }, error = function(e) {		   
			   e			   
		   }
  		 )
		} else NULL
	}
	
    if(!is.null(S0) && (inherits(S0,"error") || S0$convergence <= 0L)){
	   if(pl > 0L) { 
		 cat("Minimization by `",fun.name,"` did not converge")
		 if(!is.null(S0$status))
			cat(" (",S0$status,")")
		 cat(".","\n")
	   }
		method <- method[-1]
		S0 <- NULL
    }
	
	if(is.null(S0)){	  	
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
				if(!inherits(S0,"error") && S0$convergence > 0L) {				
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
					   	  "status"=S0$convergence,
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
	 		S0 <- structure(c(S0,qd[[1]]),
					 Sigma = attr(qd,"Sigma"),
				   class = "QSResult")				 	
	 	} else { 
			message(qd$message)
			return(structure(S0, error = qd))
		}
	  }
    }
	
	if(verbose)
	  cat("Successful minimization by: ",fun.name,"\n\n")  
  		
    return(S0)   
}

#' @name qle
#'
#' @title Simulated quasi-likelihood parameter estimation
#'
#' @description  This is the main function using a simulated quasi-likelihood estimation (QLE) approach. 
#' 
#' @param qsd			object of class \code{\link{QLmodel}}
#' @param sim		    simulation function (see details)
#' @param ...			arguments passed to function `\code{sim}`
#' @param nsim			optional, number of simulation replications at each new sample point,
#'  					`\code{qsd$nsim}` (default)
#' @param x0 			optional, numeric vector of starting parameters
#' @param Sigma			optional, constant variance matrix estimate of statistics (see details) 
#' @param global.opts	options for global search phase
#' @param local.opts	options for local search phase
#' @param method		vector of names of local search methods	
#' @param qscore.opts   list of control arguments passed to \code{\link{qscoring}}
#' @param control		list of control arguments passed to any of the routines defined in `\code{method}` 
#' @param errType		type of prediction variances (see details)
#' @param pl			print level, use \code{pl}>0 to print intermediate output
#' @param cl			cluster object, \code{NULL} (default), see \code{\link[parallel]{makeCluster}} 
#' @param iseed			integer seed, \code{NULL} (default) for no seeding of the RNG stream for each worker
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
#' 	  \item{minima}{ matrix of all local minima found so far}
#' 	  \item{score}{ quasi-score or gradient of MD}
#' 	  \item{convergence}{ logical, whether the iterates converged, see details} 	  
#' 
#'  Attributes: 	 
#'  \item{final}{ final minimization results of criterion function (see \code{\link{searchMinimizer}}) }
#'  
#'  \item{optInfo}{ a list:}
#'    \itemize{
#'    \item{x0:}{ starting vector (parameter)}
#' 	  \item{W:}{ final weighting matrix for average approximation of variance matrix and
#' 				 quasi-information matrix at \code{theta} for local sampling of new candidate points}
#'    \item{theta:}{ parameter corresponding to \code{W}} 
#' 	  \item{last.global:}{ logical, whether last iteration sampled a point globally}
#' 	  \item{minimized:}{ whether last minimization was successful}
#' 	  \item{useCV:}{ logical, whether the CV approach was applied}
#' 	  \item{method:}{ final search method applied}
#'    \item{nsim:}{ number of simulation replications at each point}
#' 	  \item{iseed}{ the seed given by the user, see arguments}
#'  }
#'    
#' 
#' @details
#'  The function sequentially estimates the unknown model parameter. Basically, the user supplies a simulation function `\code{sim}`
#'  which must return a vector of summary statistics (as the outcome of model simulations) and expects a vector of parameters
#'  as its first argument. Further arguments can be passed by the `\code{\ldots}` ellipsis. The object
#'  `\code{qsd}` aggregates the type of variance matrix approximation, the data frame of observed and simulated data, the
#'  initial sample points and the covariance models of the involved statistics (see \code{\link{QLmodel}}). In addition, it defines
#'  the criterion function by `\code{qsd$criterion}`, which is either used to monitor the sampling process or minimized itself. The user
#'  also has the option to choose among different types of prediction variances: either "\code{kv}" (kriging variances), "\code{cv}"
#'  (cross-validation variances) or the maximum of both, by "\code{max}", are available.
#' 
#'  \subsection{Criterion functions}{  
#'  The QD criterion function follows the quasi-likelihood estimation principle (see vignette) and seeks a solution to the quasi-score
#'  equation. Besides, the Mahalanobis distance (MD) as an alternative (simulation-based) criterion function has a more direct
#'  interpretation. It can be seen as a (weighted or generalized) least squares criterion
#'  depending on the employed type of variance matrix approximation. For this reason, we support several types of variance matrix
#'  approximations. In particular, given `\code{Sigma}` and setting `\code{qsd$var.type}` equal to "\code{const}" treats `\code{Sigma}`
#'  as a constant estimate throughout the whole estimation procedure. Secondly, if `\code{Sigma}` is supplied and used as
#'  an average variance approximation (see \code{\link{covarTx}}), it is considered an initial variance matrix approximation and
#'  recomputed each time an approximate (local) minimizer of the criterion function is found. This is commonly known as an iterative update
#'  strategy of variance matrices in the context of GMM estimation. Opposed to this, setting `\code{qsd$var.type}` equal to
#'  "\code{kriging}" corresponds to continuously updating the variance matrix each time a new criterion function value is
#'  required at any point of the parameter space. In this way the algorithm can also be seen as a simulated version of a least squares
#'  method or even as a special case of a \emph{simulated method of moments} (see, e.g. [3]). Note that some input combinations
#'  concerning the variance approximation types are not applicable since the criterion "\code{qle}", which exploits the
#'  QD criterion function, does not use a constant variance at all. 
#'  }
#'       
#'  \subsection{Sampling new points}{  
#'  Our QLE approach dynamically switches from a \emph{local} to a \emph{global search phase} and vise versa for sampling new promising
#'  candidates for evaluation, that is, performing new simulations of the statistical model. Depending on the current value of the criterion
#'  function three different sampling criteria are used to select next evaluation points which aim on potentially improving the quasi-score
#'  or criterion function approximation. Values less than a user defined tolerance `\code{local.opts$ftol_abs}` trigger the local
#'  phase for sampling and searching. In this case it is assumed to have found a local minimizer such that it seems
#'  to be worth sampling in its local vicinity. Then local searches mostly try to improve on the current best parameter estimate.
#'  Either the next evaluation point is selected according to a weighted minimum-distance criterion (see [2] and the vignette),
#'  for the choice `\code{nextSample}` equal to "\code{score}", or by maximizing the weighted variance of the quasi-score vector in
#'  case `\code{nextSample}` is equal to "\code{var}". In all other cases, for example, if identifiable roots of the QS could not be found
#'  or the (numerical) convergence of the local solvers failed, the global phase of the algorithm is invoked and selects new potential
#'  candidates accross the whole search space based on a weighted selection criterion. This assigns large weights to candidates
#'  with low criterion function values and vise versa. During both phases the cycling between local and global candidates is
#'  controlled by the weights `\code{global.opts$weights}` and `\code{locals.opts$weights}`, respectively. (If the algorithm is within the 
#'  global phase at last iteration, then finally a global exhaustive search (if "\code{direct}" is part of the argument `\code{method}`)
#'  for the best minimizer is done.) Besides this, the smaller the weights the more the candidates tend to be globally selected and vice versa.
#'  Within the local phase, weights approaching one result in selecting candidates close to the current minimizer of the criterion
#'  function. Weights approximately zero maximize the minimum distance between candidates and previously sampled points and
#'  thus densely sample the search space almost everywhere if the algorithm is allowed to run infinitely. The choice of weights
#'  is somewhat ad hoc but may reflect the users preference on guiding the whole estimation more biased towards either a local
#'  or global search. In addition the local weights can be dynamically adjusted if `\code{useWeights}` is \code{FALSE}
#'  depending on the current progress of estimation. In this case the first weight given by `\code{locals.opts$weights}` is 
#'  used as a starting point for this kind of adjustment.   
#'  
#'  If any of the termination criteria is met, in conjunction with a neglectable value of the criterion function, we
#'  say that the algorithm successfully terminated and converged to a reasonable parameter estimate. In addition, we 
#'  can access the precision of the found estimator by a goodness-of-fit test (see \code{\link{qleTest}}). Also
#'  in case of multiple roots, stored in the matrix of found roots `\code{minima}`, the function
#'  \code{\link{checkMultRoot}} may be helpful to select the \emph{best} of these roots (see vignette). If we wish to improve the
#'  final estimate the algorithm allows for a simple warm start strategy though not yet as an fully automated procedure. The
#'  algorithm can be easily restarted based on the final result of the preceeding run. We only need to extract the object
#'  `\code{OPT$qsd}` as an input argument to function \code{\link{qle}} again.  
#'  }
#' 
#'  \subsection{Some notes}{   
#' 	For a 2D parameter estimation problem the function can visualize the sampling and selection process, which
#'  requires an active 2D graphical device in advance. The function can also be run in an cluster environment
#'  using the `parallel` package. Make sure to export all functions to the cluster environment `\code{cl}` beforehand,
#'  loading required packages on each cluster node, which are used in the simulation function
#'  (see \code{\link[parallel]{clusterExport}} and \code{\link[parallel]{clusterApply}} in \code{parallel} package).
#'  If no cluster object is supplied, a local cluster is set up based on forking (under Linux) or as a socket connection
#'  for other OSs. One can also set an integer seed value `\code{iseed}` to initialize each worker, see \code{\link[parallel]{clusterSetRNGStream}},
#'  for reproducible results of estimation in case a local cluster is used, i.e. \code{cl=NULL} and option \code{mc.cores>1}. If
#'  using a prespecified cluster object \code{cl}, then the user is responsible for seeding whereas the seed can be stored
#'  in the return value as well, see attribute `\code{optInfo}$iseed`.     
#'  }
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
#'	 \item{\code{ftol_abs}:}{ upper bound on function criterion, values smaller than this trigger the local phase
#' 					  following the reasoning that the algorithm has found an approximate root of the quasi-score
#' 					  (or the gradient of the Mahalanobis distance).}
#'   \item{\code{eta}:}{ values for decrease and increase of the local weights, which is intended to faciliate convergence
#' 		 while sampling new points more and more around the current best parameter estimate.} 				
#'   \item{\code{nfail}:}{ maximum number of consecutive failed iterations}
#'   \item{\code{nsucc}:}{ maximum number of consecutive successful iterations}
#'   \item{\code{nextSample}:}{ either "\code{score}" (default) or "\code{var}" (see details)} 
#'   }
#' 
#'  The following controls `\code{global.opts}` for the global search phase are available:   
#' 	\itemize{
#'   \item{\code{stopval}:}{ stopping value for criterion function, if `\code{stopval>0}` then the main iteration terminates
#' 				     as soon as the function value drops below. This might be preferable to a time consuming
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
#'   \item{\code{NmaxRel}:}{ maximum number of consecutive iterates for `\code{xtol_rel}`}
#'   \item{\code{NmaxCV}:}{ maximum number of consecutive iterates for `\code{C_max}`}
#'   \item{\code{NmaxSample}:}{ maximum number of consecutive iterates for `\code{sampleTol}`}
#'   \item{\code{NmaxLam}:}{ maximum number of consecutive iterates for `\code{local.opts$lam_max}` }
#'   \item{\code{Nmaxftol}:}{ maximum number of consecutive iterates for `\code{local.opts$ftol_rel}`}
#'  }
#'  
#' 
#' @seealso \code{\link{mahalDist}}, \code{\link{quasiDeviance}}
#' 
#' @example /inst/examples/qle.R
#' 
#' @author M. Baaske
#' @rdname qle
#' 
#' @useDynLib qle, .registration = TRUE, .fixes = "C_"
#' @export 
#' 
#' @import parallel stats
#' @importFrom graphics points
qle <- function(qsd, sim, ... , nsim, x0 = NULL, Sigma = NULL,
				 global.opts = list(),local.opts = list(),
				  method = c("qscoring","bobyqa","direct"),
				   qscore.opts = list(), control = list(),
				    errType = "kv", pl = 0L, 
					 cl = NULL, iseed = NULL, plot=FALSE)
{		
	# print information 	
	.printInfo = function(k,n){		
		if(pl > 0L) {
			cat("\n\n")			
			cat("Evaluation: ",k+n,"\n")
			cat(paste0("[global=",k,", local=",n,"]"),"\n\n")			
			cat("Criterion: ",qsd$criterion,"\n\n")			
			cat("Current solution: \n\n")
			print.default(formatC(signif(xt, digits=8), digits=8, format="fg", flag="#"),
					print.gap = 4, quote = FALSE)
			cat("\n")
			cat("Objective: ",formatC(ft, digits=4, format="e", big.mark=","),"\n\n")			
		}		
	}

	.showConditions = function() {
		if(pl > 1) {
			cat("Update parameter: \n\n")				
			cat("global:..........",ifelse(status[["global"]],"yes","no"),"\n")
			cat("minimized:.......",ifelse(status[["minimized"]],"yes","no"),"\n")
			cat("weight factor:...",w,"\n")
			cat("start at:........","[",formatC(signif(x[1:length(x)],digits=6),digits=6,format="fg", flag="#"),"]","\n")
			cat("solution:........","[",formatC(signif(xt[1:length(xt)],digits=6),digits=6,format="fg", flag="#"),"]","\n")
			cat("new sample:......","[",formatC(signif(Xnext[1:length(Xnext)],digits=6),digits=6,format="fg", flag="#"),"]","\n")			
			cat("\n")
			
			# other conditions
			cat("Current stopping conditions: \n\n")			
			cond <- 
			 if(status[["minimized"]] && !is.null(QD$score)) {
			  c("|score_max|" = max(abs(QD$score)))				
			 } else NULL
			cond <-
			 c(cond,
			   "lam_max"=unlist(ctls["lam_max","val"]),
			   "varTol"=unlist(ctls["C_max","val"]),
			   "ftol_rel"=unlist(ctls["ftol_rel","val"]),
			   "xtol_rel"=unlist(ctls["xtol_rel","val"]),
			   "sampleTol"=unlist(ctls["sampleTol","val"]))	
            
			print.default(format(cond, digits = 4, justify="left"),
					print.gap = 2, quote = FALSE)						
			cat("\n")
			if(pl > 2) {				
				cat("\nVariance matrix: \n\n Sigma = \n\n")
				print(attr(QD,"Sigma"))				
				if(!is.null(W)) {					
					cat("\nWeighting matrix: \n\n W = \n\n")
					print(W)					
				}
				cat("\n")
			}							
		}
	}	
	
	if(!is.numeric(pl) || pl < 0L)
	  stop("Print level `pl` must be some positive numeric value.")	
	if(missing(nsim))
	  nsim <- attr(qsd$qldata,"nsim")  	
	if(is.null(nsim) || !is.numeric(nsim))
	  stop("Number of simulations must be given.")
    
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
	xnames <- NULL
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
			Sigma <- try(solve(Sigma),silent=TRUE)
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
	qscore.opts <- .qsOpts(qscore.opts)

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
		
	# matrix of approximate roots
	loc_min <- matrix(,nrow=0L,ncol=xdim)	
	# get initial sample
	X <- as.matrix(qsd$qldata[seq(xdim)])
	xnames <- colnames(X)
	# initial sample size
	nstart <- nrow(X)	
	
	## set 'globals'
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
	
	if(any(locals$weights > 1L) || any(locals$weights < 0L))
		stop("Weights for local sampling must be
			in the interval [0,1] for sampling criterion `score`!")	
	w <- locals$weights[1]	
	mWeights <- length(locals$weights)	
	
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
				if(length(iseed) == 0L)
				   message("You did not set a seed!")
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
	status <- list("global"=FALSE, "minimized"=TRUE) 
			
	# first time CV:
	# this is also for testing 
	# before iterating many times 
	cvm <- NULL
	errId <- pmatch(errType,c("kv","cv","max"))
	if(anyNA(errId) || length(errId)!=1L)
	  stop("Invalid argument `errType`. Please choose one of `kv`, `cv`, `max`")
	useCV <- (errId > 1)	   	
		
	if(useCV) {
		if(pl > 0L)
		  cat("Update cross-valiadation covariance models...\n")
		cvm <- try(prefitCV(qsd,type=errType,cl=cl),silent=TRUE) 
		if(.isError(cvm)) {						
			cvm <- NULL
			message(.makeMessage("Prefit of CV models failed during final surrogate minimization."))			
		} 
	}
	
	# initialize
	W <- theta <- NULL
	Xnext <- xt <- x <- x0
	fs <- ft <- f <- 1E+100
	
	# first search might be with given starting Sigma 
	S0 <- searchMinimizer(x0, qsd, method, qscore.opts, control,
			Sigma=Sigma, W=NULL, theta=NULL, inverted=TRUE, cvm=cvm,
			 info=TRUE, check=FALSE, pl=pl, verbose=pl>0L)
	
	# but then reset so it can be computed again
	if(qsd$var.type != "const")
	 Sigma <- NULL
	
	# was it successful ?
	if(inherits(S0,"error") || !(S0$convergence > 0L))
	   warning("First local search failed or did not converge.")	 
				
	dummy <- 
	  tryCatch({						
		repeat{								
				# This is just for the worst cases:
				# Set current iterate to last sampled point in case of no convergence,
				# eventually a sampled point 'Xnext' also becomes a minimizer after a number
				# of steps cycling through the weights
				if(!inherits(S0,"error") && S0$convergence) {
					xt <- S0$par
					ft <- S0$value
					status[["minimized"]] <- TRUE
				} else {
					# old sampling location
					ft <- fs
					xt <- Xnext					
					status[["global"]] <- TRUE	 
					status[["minimized"]] <- FALSE
				}	
											
				# Test stop maximum iterations
				if(!((nglobal+nlocal) < maxEval && nglobal < maxIter))
				  break					
				# check stopval only for global phase
				if(status[["global"]] && ft < ctls["stopval","cond"]) 
				  break
			  	# Minimum sampling distance reached ?
			  	dm <- attr(.check.distAll(X,xTol=ctls["sampleTol","cond"]),"min")
			  	if(dm < ctls["sampleTol","val"])
				  ctls["sampleTol","val"] <- dm
				if(dm < ctls["sampleTol","cond"]) {
				  ctls["sampleTol","stop"] <- ctls["sampleTol","stop"] + 1L
				  if(ctls["sampleTol","stop"] >= globals$NmaxSample)
					break;
			    } else {  ctls["sampleTol","stop"] <- 0L }					  
			  
				# check xtol and ftol_abs
				ctls["xtol_rel","val"] <- max(abs(xt-x)/pmax(abs(x),EPS))
				if( ctls["xtol_rel","val"] < ctls["xtol_rel","cond"]) {
					ctls["xtol_rel","stop"] <- ctls["xtol_rel","stop"] + 1L
					if(ctls["xtol_rel","stop"] >= globals$NmaxRel) {				
						break
					}
				} else { ctls["xtol_rel","stop"] <- 0L }
				
				# ftol_rel (global and local) (low priority)
				ctls["ftol_rel","val"] <- abs(ft-f)/max(abs(f),EPS)
				if( ctls["ftol_rel","val"] < ctls["ftol_rel","cond"]) {
					ctls["ftol_rel","stop"] <- ctls["ftol_rel","stop"] + 1L
					if(ctls["ftol_rel","stop"] >= globals$Nmaxftol) {					
						break
					} else { ctls["ftol_rel","stop"] <- 0L }
				}				
				## choose criterion function
				QD <- criterionFun(xt)[[1]]
				 
				# generalized EVD (either based on CV or Kriging variances)
				ctls["lam_max","tmp"] <- ctls["lam_max","val"]
				ctls["lam_max","val"] <- max(geneigen(QD$varS,QD$I,
								           only.values=TRUE), na.rm=TRUE)
				
				# generalized eigenvalues
				if( ctls["lam_max","val"] < ctls["lam_max","cond"]) {
					ctls["lam_max","stop"] <- ctls["lam_max","stop"] + 1L
					if(ctls["lam_max","stop"] >= globals$NmaxLam) 
						break					
				} else { ctls["lam_max","stop"] <- 0L }
				
				# Maximum prediction variance of the quasi-score vector:
				# either CV based or evaluated by kriging variances.				
				if(qsd$krig.type == "var") {
					ctls["C_max","tmp"] <- ctls["C_max","val"]				
					ctls["C_max","val"] <- max(diag(QD$varS))
					
					test <- abs(ctls["C_max","val"]-ctls["C_max","tmp"])/ctls["C_max","tmp"]
					if(test < ctls["C_max","cond"]) {
						ctls["C_max","stop"] <- ctls["C_max","stop"] + 1L
						if(ctls["C_max","stop"] >= globals$NmaxCV) {
							ctls["C_max","val"] <- test		
							break;
						}
					} else { ctls["C_max","stop"] <- 0L }
				}
								
				# current iteration did not stop 
				# -> next is local/global ?				
				if(ft < locals$ftol_abs) {					    
					# W is weighting matrix for
					# variance average approximation
					# and its inverse is used as the 
					# variance of theta for sampling from MVN 
					W <- QD$I		
					theta <- xt
					# found minimizer/root?
					loc_min <- rbind(loc_min,xt)
					if(status[["global"]]) {
						w <- locals$weights[1]
						status[["global"]] <- FALSE
					}												
				} else status[["global"]] <- TRUE		
				
				Snext <- 
					tryCatch({						
						if(!status[["global"]])
						{
							# generate local candidates							
							Y <- nextLOCsample(W,theta,
									locals$nsample,lb=qsd$lower,
									ub=qsd$upper,pmin=locals$pmin,invert=TRUE)
							
							if(.isError(Y)) {
								message(paste0("Sampling local candidates failed. Try global sampling."))
								status[["global"]] <- TRUE
							} else {							
								# get min distances
								nlocal <- nlocal + 1L
								dists <- .min.distXY(X,Y)
								# remove too close sample candidates							
								idx <- which(dists < eps)					
								if(length(idx) > 0L)	{
									Y <- Y[-idx,,drop=FALSE]
									if(nrow(Y) < floor(0.1*length(dists))) {											
										msg <- paste("Number of local candidates is less than 10% of sample size:, "
												,nrow(Y)," try global sampling now.")
										message(msg)
										status[["global"]] <- TRUE
									}
									dists <- dists[-idx]
								}																
								
								if(!status[["global"]])
								{		 
									 dmin <- min(dists)
									 dmax <- max(dists)									 
									 id <- 
									   switch(
									    locals$nextSample,
									     "score" = {			
											 # use user defined weights
											   if(locals$useWeights) {										
												  k <- (nlocal+1L) %% mWeights
												  w <- ifelse(k != 0L,
														  locals$weights[k],
														  locals$weights[mWeights] )
											   } else {
												 if( ft < (f-1e-4) ||
														 (ctls["lam_max","val"] < ctls["lam_max","tmp"] &&
														  ctls["C_max","val"] <= ctls["C_max","tmp"]) )
												 {
													 ctls["nfail","val"] <- 0L
													 ctls["nsucc","val"] <- ctls["nsucc","val"] + 1L												
												 } else {
													 ctls["nsucc","val"] <- 0L
													 ctls["nfail","val"] <- ctls["nfail","val"] + 1L												
												 }
												 
												 # update weights									
												 if(ctls["nfail","val"] > 0L && 
												  !(ctls["nfail","val"] %% ctls["nfail","cond"])){											 											 
													 w <- max(w-locals$eta[1],0.005)
												 } else if(ctls["nsucc","val"] > 0L && 
														 !(ctls["nsucc","val"] %% ctls["nsucc","cond"])){										 
													 w <- min(w+locals$eta[2],0.995)
												 }									 										 
											  }						
											  # minimize ballanced criterion
										  	  # criterion funtion values at candidates
										  	  # Sigma is re-calculated here at theta)
											  fd <- criterionFun(Y,value.only=2)
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
											  fd <- criterionFun(Y,value.only=3)								  
											  dw <- if(abs(dmax-dmin) < EPS) 1		
													   else (dists-dmin)/(dmax-dmin)
											  which.max( fd*dw )
										  })									
								} 								
							}							
						} # end local sampling
						
						if(status[["global"]])
						{			
							nglobal <- nglobal + 1L
							# sample new candidates
							Y <- sapply(seq_len(ncol(X)),
							  		function(i) {
									  runif(globals$nsample,qsd$lower[i],qsd$upper[i])
									}
							)							
							colnames(Y) <- xnames							
							# reset: no weighting in global phase
							W <- theta <- NULL				
							# distances check
							dists <- .min.distXY(X,Y)
							idx <- which(dists < eps)
							if(length(idx) > 0L)	{
								Y <- Y[-idx,,drop=FALSE]
								if(nrow(Y) < floor(0.1*length(dists))) {										
							 	   message("No further sampling locations could be found. Stop sampling now.")									 											
								   break;
								}
								dists <- dists[-idx]
							}							
							# quasi-deviance or Mahalanobis distance as
							# a criterion for global sampling
							fval <- criterionFun(Y,value.only=TRUE)									  																		
							
							# next sampling location							
							fmin <- min(fval)
							fmax <- max(fval)
							k <- (nglobal+1L) %% mWeightsGL
							w <- ifelse(k != 0L, globals$weights[k], globals$weights[mWeightsGL] )
							fd <- if(abs(fmax-fmin) < EPS) 1 
							      	else (fval-fmin)/(fmax-fmin)									
							
							id <- which.max( exp(-w*fd) * dists )				
						
						} # end global					
						
						# new sample point and value
						list("par"=Y[id,],"value"=fd[id])
						
					}, error = function(e) {
						msg <- .makeMessage("Sampling new candidates failed: ",
								  conditionMessage(e))
						stop(msg)
					}
				) # end sampling tryCatch			
					
				# next sampling location
				fs <- Snext$value
				Xnext <- rbind(Snext$par)
				# plot iterates (2D) 
				if(plot && xdim <= 2L) {
					p <- rbind(Xnext,xt)
					if(xdim == 1L) 
					 p <- cbind(p,c(0,0))					 
					cols <- if(status[["global"]]) "blue" else "green"
					try(points(p[1,,drop=FALSE],pch=8,cex=1,col=cols,bg=cols))
					try(points(p[2,,drop=FALSE],pch=21,cex=0.15,col="magenta",bg="magenta"))
				}		
						
				# show info
				.printInfo(nglobal,nlocal)				
				# print stopping conditions 
				.showConditions()
				
				if(pl > 0L) {
				 cat("---\n\n")	
			 	}
				
				# simulate at new locations and update (fit): qsd, data,				
				qsd <- updateQLmodel(qsd, Xnext, simFun, nsim=nsim,						 
						 fit=TRUE, cl=cl, verbose=pl>0L)
						
				# check results of kriging
				if(.isError(qsd)){
					msg <- paste(c("Cannot simulate at candidate point: \n\t ",
								format(Xnext[1:ncol(Xnext)], digits=6, justify="right"),"\n"),
							collapse = " ")
					e <- attr(qsd,"error")
					if(inherits(e,"error"))
					  msg <-  c(msg, conditionMessage(e))
					message(msg)
					return(.qleError(message=msg,error=qsd))
				}	
								
				# update iterates				
				# if global then start minimization at random sample points
				# and eventually start very near a (local) minimum of QD in case
				# the weights w are relatively large;
				# (QS iteration starts at previously found approx. root)
				f <- ft				
				x <- if(status[["global"]])	Snext$par else xt				
				# update X sample			
				X <- as.matrix(qsd$qldata[seq(xdim)])	
				
				# refit
				if(useCV) {
					if(pl > 0L)
					  cat("Update cross-validation covariance models...\n")
					cvm <- try(prefitCV(qsd,type=errType,cl=cl),silent=TRUE) 
					if(.isError(cvm)) {						
						cvm <- NULL
						message("Prefit of CV models failed during final surrogate minimization.")			
					} 
				}	
				# minimizing surrogate, root finding or local/global search
				# if maxeval reached, then also return criterion function results;
				# Also if "direct" is supplied, then use it in case of last global sampling
				info <- (nglobal==maxIter || ((nglobal+nlocal)>=maxEval))				
				if(info && status[["global"]]) {
					if("direct" %in% method)
					  method <- "direct" 
				}
				S0 <- searchMinimizer(x, qsd, method, qscore.opts, control,
						Sigma=Sigma, W=W, theta=theta, inverted=TRUE, cvm=cvm,
						info=info, check=FALSE, pl=pl, verbose=pl>0L)
			}
			# end main loop
			NULL
		}, error = function(e) {
			msg <- .makeMessage("Current iteration stopped unexpectedly: ",
					conditionMessage(e))
			message(msg)
			# return last state of loop
			.qleError("qle",msg,sys.call(),
					  error=list("par"=xt,
							     "value"=ft,
								 "ctls"=ctls,	
								 "qsd"=qsd,	
								 "cvm"=cvm,
								 "convergence"=FALSE,
								 "why"=NULL,
								 "minima"=loc_min,
						    	 "final" = S0,
						   		 "x0"=x0,
								 "W"=W,
								 "theta"=theta,
								 "last.global"=status[["global"]],
								 "minimized"=status[["minimized"]],
								 "useCV"=useCV,
								 "nsim"=nsim,
								 "iseed"=iseed))					
		}, finally = {
		  if(noCluster) {
			if(inherits(try(stopCluster(cl),silent=TRUE),"try-error"))
			  message("Error in stopping cluster.")
		  }
		}
	) # end outer tryCatch	
	
	# stop on error 
	if(.isError(dummy))
	 return(dummy)
	 
	## only for estimte theta=(xt,ft)	
	ctls["stopval",c(2,4)] <- c(ft,ft < ctls["stopval","cond"])	
	ctls["maxiter",c(2,4)] <- c(nglobal,nglobal >= maxIter)
	ctls["maxeval",c(2,4)] <- c(nglobal+nlocal, nglobal+nlocal >= maxEval)
	
	# remove `nfail`, `nsucc`
	ctls <- ctls[1:8,c(1,2,4)]	
	
	# last fitted qsd model
	if(.isError(S0)) {
	 message("Last search results have errors. Please see `attr(,\"final\")`.")
	} else {	  	
		ctls <- rbind(ctls,
		    	as.data.frame(cbind("cond"=qscore.opts$score_tol,
			 			"val"=max(abs(S0$score)) ,"stop" = 0L),
		   		     row.names = ifelse(qsd$criterion == "qle","score","grad"),
	   			check.names = FALSE))							  	
	}
	# why stopped		
	arg.names <- row.names(ctls[which(ctls[,"stop"] > 0L),])
				
	return (
		structure(
		    list("par"=xt,
				 "value"=ft,
				 "ctls"=ctls,
			 	 "qsd"=qsd,
				 "cvm"=cvm,
				 "why"=arg.names,
				 "minima"=loc_min,				 
				 "convergence"=(status[["minimized"]] && length(arg.names) > 0L)),
		 	"final" = S0,
		 	"optInfo" = list("x0"=x0,
					 	     "W"=W,
						     "theta"=theta,						   
				    	     "last.global"=status[["global"]],
					  	     "minimized"=status[["minimized"]],
					  	     "useCV"=useCV,
					  	     "method"=S0$method,
					  	     "nsim"=nsim,
							 "iseed"=iseed),
		 class = "qle", call = sys.call()))	
}

#' @name print.qle
#' 
#' @title print results of class \code{qle}
#' 
#' @description S3 method to print the results of \code{\link{qle}}.
#' 
#' @param x      object of class \code{qle} from call to function \code{\link{qle}}
#' @param pl 	 numeric (positive) value, the print level (higher values give more information)
#' @param digits number of digits to display
#' @param ... 	 ignored, additional arguments
#' 
#' @rdname print.qle
#' @method print qle
#' @export 
print.qle <- function(x, pl = 1, digits = 5,...){	
	if(.isError(x)) 
	  stop(.makeMessage("Estimation result has errors.","\n"))	
	if(!inherits(x,"qle"))
	  stop("Method is only for objects of class `qle`.")
    if(!is.numeric(pl) || pl < 0L )
	  stop("Print level must be a positive numeric value.")
  
	if(pl > 0L) {
	  	if(x$qsd$criterion == "mahal")
		 cat("Mahalanobis distance value: \n\n",x$value,"\n\n")
		else			
		 cat("Quasi-deviance value: \n\n",x$value,"\n\n")
	  	
		cat("Estimated parameter:\n\n")
		print.default(formatC(signif(x$par, digits = digits), digits = digits, format="fg", flag="#"),
				print.gap = 4, quote = FALSE)
		
		cat("\n")
		if(pl > 1) {
			cat("Convergence: ", x$convergence, collapse =" ","\n\n")
			if(x$convergence) {
				by <- x$ctls[x$why,"val"]
				names(by) <- x$why		
				print.default(formatC(by, format="e", digits=digits), right=FALSE, print.gap=2,
						quote=FALSE)
			} else cat("(No convergence criteria matched.)\n\n")
		} 
	}
  	if(pl > 1) {
		cat("\n")	
		nsim <- unlist(x$ctls["maxeval","val"])*attr(x,"optInfo")$nsim	
		cat("Evaluations: ",unlist(x$ctls["maxeval","val"])," (simulations: ",nsim,")\n\n")
		cat("Variance approximation by: ",x$qsd$var.type,"\n")	
	}
	if(pl > 2) {
		if(!.isError(attr(x,"final"))) {
			cat("\n\n ***  Final search results *** \n\n\n")			
			print(attr(x,"final"))
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
#' @description S3 method to print the results of \code{\link{qscoring}}.
#' 
#' @param x  		object of type \code{QSResult}, obtained from \code{\link{qscoring}}
#' @param pl		numeric positive value, the print level (higher values give more information) 
#' @param digits 	number of digits to display
#' @param ... 	    ignored, additional arguments
#' 
#' @rdname print.QSResult
#' @method print QSResult
#' @export 
print.QSResult <- function(x, pl = 1, digits = 5,...) {	
	if(.isError(x)) 
	  stop(.makeMessage("Quasi-scoring result has errors.","\n"))
	if(!inherits(x,"QSResult"))
	  stop("Method is only for objects of class `QSResult`.")	
    if(!is.numeric(pl) || pl < 0L )
	 stop("Print level must be a positive numeric value.")
 
	if(pl > 0L) {
		cat("Local method: ",x$method,"\n\n")		
		cat("Solution: \n\n")
		print.default(formatC(signif(x$par, digits = digits), digits = digits, format="fg", flag="#"),
				print.gap = 4, quote = FALSE)
		cat("\n")
		cat("Objective:\n\n",x$value,"\n\n")
		cat("Iterations....",x$iter,"\n")		
		cat("Convergence...",x$convergence,"\n")				
		if(!is.null(x$score)){
			cat("\nStopped by: \n\n",unlist(strsplit(paste(x$message, "\n" ),':')), fill=TRUE )
			if(x$criterion == "qle") {
				cat("Quasi-score:\n\n")
			}  else cat("Gradient:\n\n")
			print.default(formatC(signif(x$score, digits=8), digits=8, format="fg", flag="#"),
					print.gap = 4, quote = FALSE)				
		}
		cat("\n")
	}
	if(pl > 2) {
		Sigma <- attr(x,"Sigma")		
		if(!is.null(Sigma)) {	
			cat("\nVariance matrix approximation: \n\n Sigma = \n\n")
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
#' @description Generate a random sample of points. 
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
#' 	 For an example see file "normal.R".		
#'  
#' @author M. Baaske
#' @rdname nextLOCsample
#' 
#' @importFrom mvtnorm pmvnorm rmvnorm
#' @export
nextLOCsample <- function(S, x, n, lb, ub, pmin = 0.05, invert = FALSE) {
	if(invert) 
	  S <- try(solve(S),silent=TRUE)
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
#' @description The function solves the quasi-score equation by a 
#'  a root finding algorithm similar to Fisher's scoring iteration.
#'
#' @param qsd    	object of class \code{\link{QLmodel}}
#' @param x0		(named) numeric vector, the starting parameter
#' @param opts		quasi-scoring options (see details)
#' @param Sigma	    a pre-specified variance matrix estimate
#' @param ...	    further arguments passed to function \code{\link{covarTx}}
#' @param inverted  currently ignored
#' @param check		logical, \code{TRUE} (default), whether to check input arguments
#' @param cvm		list of covariance models for cross-validation (see \code{\link{prefitCV}})
#' @param pl	    numeric, print level, use \code{pl}>0 to give intermediate output  	
#' @param verbose   \code{FALSE} (default), otherwise print intermediate output
#'
#' @return List of results of quasi-scoring iteration.
#'  \item{status}{ integer, why scoring iterations stopped}
#'  \item{message}{ string, corrsponding to `\code{status}`}
#'  \item{iter}{ number of iterations}
#'  \item{value}{ quasi-deviance value}
#'  \item{par}{ solution vector}
#'  \item{score}{ quasi-score vector}
#'  \item{I}{ quasi-information matrix}
#'  \item{start}{ starting point}
#'  \item{convergence}{ integer, for convergence=1 and local convergence=10, no convergence=-1}
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
#'  The quasi-scoring iteration covers both kinds of prediction variances, kriging-based and those by a cross-validation
#'  approach, which account for the uncertainty induced by the quasi-score approximation model. By default kriging variances
#'  at each iterate are included in the computation. If fitted covariance models `\code{cvm}` are supplied by the user
#'  in advance, the variances of prediction errors of each statistic are separately evaluated by the proposed cross-validation
#'  approach (see vignette) at each new iterate. For the price of relatively high computational costs those prediction variances
#'  are intended to increase the robustness against false roots due to simulation and approximation errors of the quasi-score function.
#' 
#'  Opposed to this, the user also has the option to carry out a "pure version" of quasi-scoring without accounting for
#'  these errors. This can be set earlier as an option in \code{\link{QLmodel}}. See also \code{\link{covarTx}} and
#'  \code{\link{mahalDist}} for details on how to choose the variance matrix approximation of the statistics.
#'  }
#' 
#'  The following algorithmic options, which can be set by `\code{opts}`, are available:
#'  \itemize{
#'   	\item{\code{ftol_stop}:}{ minimum value of quasi-deviance for stopping the scoring iteration}
#' 		\item{\code{xtol_rel}, \code{ftol_rel}:}{ see \code{\link{qle}} }
#' 		\item{\code{grad_tol}:}{ upper bound on the quasi-score vector components,
#' 				 testing for a local minimum of the quasi-deviance in case of a line search failure}
#' 		\item{\code{score_tol}:}{ upper bound on the quasi-score vector components, testing for an approximate root}
#'		\item{\code{slope_tol}:}{ upper bound on the 2-norm of the quasi-score vector, testing for an approximate descent step}
#'      \item{\code{maxiter}:}{ maximum allowed number of iterations}
#' 	    \item{\code{pl}:}{ print level (>=0), use \code{pl}=10 to print individual
#' 							 iterates and further values}
#'  } 
#'
#' @example  inst/examples/scoring.R
#' 
#' @author M. Baaske
#' @rdname qscoring
#' @export
qscoring <- function(qsd, x0, opts = list(), Sigma = NULL, ...,
			    	  inverted = FALSE, check = TRUE, cvm = NULL, 
				 	   pl = 0L, verbose = FALSE)
{	
	if(check)
	 .checkArguments(qsd,x0,Sigma)
 	stopifnot(is.numeric(pl) && pl >= 0L )
 
  	if(qsd$criterion!="qle")
	  stop("Quasi-scoring is only valid for criterion `qle`.")
  	X <- as.matrix(qsd$qldata[seq(attr(qsd$qldata,"xdim"))])
  		
	if(qsd$var.type != "kriging" && is.null(Sigma)){
		# Only mean covariance matrix is estimated here. 
		# Adding prediction variances (kriging/CV) at C level		
		if(qsd$var.type %in% c("wcholMean","wlogMean")){
			nms <- names(list(...))
			if(!all( c("W","theta") %in% nms))
				warning("Found `var.type`=\"",qsd$var.type, "\" but no weighting matrix `W` or estimate `theta` was supplied!.")		
		}
		Sigma <- covarTx(qsd,...,cvm=cvm)[[1]]$VTX		
	} 
	
	# set quasi-scoring options
	opts <- .qsOpts(opts,pl)
	qlopts <- list("varType"=qsd$var.type,					  
				   "useCV"=!is.null(cvm),
				   "useSigma"=FALSE)
		   
	try(.Call(C_QSopt, x0, qsd, qlopts, X, Sigma, cvm, opts), silent=TRUE)	
}


## intern
## Conduct next simulations,
## and update covariance models
updateQLmodel <- function(qsd, Xnew, sim, nsim, fit = TRUE,
							 cl = NULL, verbose = FALSE ){
	# new simulations, qsd$nsim is default
	if(verbose)
	 cat("Simulating the model...\n")
	newSim <-
	   tryCatch(simQLdata(sim, nsim=nsim, X=Xnew, cl=cl, verbose=verbose),
			error = function(e) {
				msg <- .makeMessage("All simulations failed: ",
						conditionMessage(e))
				message(msg)
			   .qleError(message=msg,call=match.call(),error=e)
			}
	)		 
	if(.isError(newSim))
	  return(newSim)  
	# get new data
	if(verbose)
	 cat("Initialize current iteration...\n")
	nextData <-
		tryCatch(
		  setQLdata(newSim,
				    X=Xnew,
				    chol=(qsd$var.type!="const"),
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
