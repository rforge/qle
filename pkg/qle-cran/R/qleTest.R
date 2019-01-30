# Copyright (C) 2018 Markus Baaske. All Rights Reserved.
# This code is published under the GPL (>=3).
#
# File: 	krige.R
# Date:  	14/03/2018
# Author: 	Markus Baaske
#
# Monte Carlo based testing procedure using an approximate
# efficient score test statistic (here: quasi-deviance) or a generalized
# least squares criterion as the test statistic for testing the goodness-of-fit

# collect test results
.qleTest <- function(obj,alpha = 0.05) {	
	sb <- attr(obj,"sb")
	Sb <- attr(obj,"Sb")	
	pb <- (1+sum(Sb>=sb))/(length(Sb)+1)  # similar to : 1-ecdf(Sb)(sb)
	qt <- try(quantile(Sb,1-alpha,na.rm=TRUE),silent=TRUE)
	if(inherits(qt,"try-error"))
	 stop("Could not compute quantiles of test statistic.")	
	
	ans <- list()
	ans$param <- cbind(obj[c("par","se","rmse","bias","mean")])
	dimnames(ans$param) <- list(row.names(obj),
			c("Estimate", "Std. Error", "RMSE", "Bias", "Mean"))
	
	tname <- attr(obj,"test")	
	ans$test <- noquote(cbind(sb, pb))# , ifelse(degf>0,pchisq(sb, degf, lower.tail = FALSE),"****")))	
	dimnames(ans$test) <- list("Test E(Q)=0:  ", c("s value", "Pr(>s)")) #, "Pr(>Chisq)"))
	ans$Stest <- noquote(paste0("Bootstrap ",ifelse(tname=="mahal", "LS criterion test:", "Score-test:"))) 	
	ans$tname <- tname
	attr(ans$test,"Sb") <- Sb
	attr(ans$test,"qt") <- qt
	attr(ans$test,"passed") <- (sb < qt)
	
	return(
	  structure(ans,
	    class=c("qleTest"))
	)	
}

# choose the best root, if any, according to the  
# criteria (see vignette) and the smallest value
# of the maximum of the quasi-score vector
.evalBestRoot <- function(dm, opts)
{
	stopifnot(is.data.frame(dm))	
	# remove all rows of criteria wher any NA can be found
	notNa <- !(apply(dm,1,anyNA))	
	if(sum(notNa) < nrow(dm)){		
	  message("`NA`s removed from quasi-score solutions.")
	} else if(sum(notNa) == 0L) {
	  warning("Cannot select a best solution. All rows in data frame of estimated parameters contain `NA`s.")
	  return( structure(row.names(dm),"id"=NA,"valid"=FALSE))
    }	
	# matrix to check the parameter at least being a plausible root
	A <- matrix(FALSE,nrow(dm),3)
	A[notNa,] <-
	 cbind(as.numeric(dm[notNa,"minor"]) == 0, 
		   as.numeric(dm[notNa,"value"]) < opts$ftol_abs,
		   as.numeric(dm[notNa,"|score_max|"]) < opts$score_tol)
	
	ki <- 1L 
	ok <- which(A[,1]==TRUE) 										# all where 'Iobs' pos.def. 
  	if(!is.numeric(ok) || length(ok) == 0L){
	  ki <- c(ki,2L)												# exclude 'minor' and 'det' 
	  ok <- which(notNa==TRUE)
	  message(.makeMessage("'Observed quasi-information is not positive definite so cannot use it."))	  
    }   
	notOk <- which(!apply(A[ok,2:3,drop=FALSE],1,all))
	if(is.numeric(notOk) && length(notOk)>0L){
		if(length(notOk) < length(ok)) {
			ok <- ok[-notOk]				
		} else { message(.makeMessage("Conditions 'ftol_abs' and 'score_tol' cannot be met.")) }
	}
		
	nm <- apply(dm[ok,-ki,drop=FALSE],2,which.min)							# 'minor' and 'det' removed as criteria because Iobs not pos. def.
	id <- sapply(1:length(nm),function(x) sum(nm[1:x] == nm[x]))
	maxid <- cbind(nm, id)[,2]
	id <- as.numeric(nm[which(maxid == max(maxid),arr.ind=TRUE)])
	if(length(id) > 1L){
		for(i in id)		  	
		 dimnames(dm)[[1]][ok[i]] <- paste0(c(row.names(dm)[ok[i]],"*"),collapse=" ")	   
		
	 	mi <- which.min(dm[ok[id],"value"])									# index of id best roots
		dimnames(dm)[[1]][ok[id[mi]]] <- paste0(c(row.names(dm)[ok[id[mi]]],"*"),collapse="")
	} else {
		mi <- 1L
		dimnames(dm)[[1]][ok[id]] <- paste0(c(row.names(dm)[ok[id]],"*"),collapse=" ")
	}	
	return( structure(row.names(dm), "id"=ok[id[mi]], "valid"=all((A[mi,])))) 
}

.evalRoots <- function(QD, par = NULL, opts = NULL)
{	   
	if(.isError(QD))	  	
	 return(.qleError(message="Evaluation of roots failed.",call=sys.call(),error=QD) )
	
	options <-list("ftol_abs"=1e-6, "score_tol"=1e-4)
	opts <-
	  if(is.null(opts))
	     options
	  else {	  
		  id <- which(!is.na(pmatch(names(opts),names(options))))
		  if (length(id)>0L)
			  options[names(opts[id])] <- opts[id]	
		  options
	}    
	if(is.null(par)){
	   par <- try(do.call(rbind,lapply(QD,"[[","par")),silent=TRUE)
	   stopifnot(is.matrix(par))
    }
		
	# some criteria for selection of the 'best' root
	# select the root for which most of the criteria below are minimized
	X <- lapply(QD,
			function(qd) {
			   if(.isError(qd)) return (qd)			  								
			   try(
			    c("minor"=.isPosDef(qd$Iobs) || as.integer(!all(sign(qd$Iobs) == sign(qd$I))),
				  "det"= abs(1-det(qd$I)/det(qd$Iobs)), 
				  "value"=-log(det(qd$I))+qd$value,			  
				  "|score_max|"=max(abs(qd$score),na.rm=TRUE),				
				  "lamI_min"=-min(eigen(qd$I)$values),
				  "lamIm_max"=abs(max(geneigen(qd$varS,qd$I,only.values=TRUE),na.rm=TRUE)-1)
  				  ), silent=TRUE)
			})
	
	ok <- which(sapply(X, function(x) !.isError(x) ))	
	if(length(ok) == 0L){
		msg <- .makeMessage("Cannot find or compare plausible solutions.")
		message(msg)
		return(.qleError(message=msg,call=sys.call(),error=X))
  	}		
	id <- 1L
	dm <- NULL
	M <- as.data.frame(do.call(rbind,X[ok]))
	row.names(M) <- if(!is.null(par)) row.names(par)	      		   
	
	# return data frame of X[ok]
	dm <- cbind(par[ok,,drop=FALSE],M)
	isErr <- which(!(1:length(X) %in% ok))
	if(length(isErr)>0L) {
		attr(dm,"X") <- X[isErr]
		attr(dm,"isErr") <- isErr
	}
	
	nms <- try(.evalBestRoot(M,opts),silent=TRUE)
	if(!.isError(nms)){
		row.names(dm) <- nms
		id <- attr(nms,"id")
		attr(dm,"id") <- id 
		attr(dm,"valid") <- attr(nms,"valid")		
	} else {
		attr(dm,"par") <- par[1L,]
		msg <- .makeMessage("Failed to select best solution.")
		message(msg)
		attr(dm,"error") <- structure(list(message = msg, call = match.call()), id=id, nms=nms)
		return (dm)
	}   	
	
	if(anyNA(id) || length(id) != 1L){	 
	 msg <- .makeMessage("Cannot select any parameter as a root of quasi-score.")
	 message(msg)
	 attr(dm,"par") <- par[1L,]
	 attr(dm,"error") <- structure(list(message = msg, call = match.call()),id=id)
 	} else	{
		attr(dm,"par") <- par[id,]
	}	
	return (dm)	  
}


#' @name		checkMultRoot
#' 	
#' @title       Assess plausibility of parameter estimates
#' 
#' @description		Check out and compare estimated roots of the quasi-score vector
#' 
#' @param est 		object of class \code{qle}, the estimation results from function \code{\link{qle}}
#' @param par	    list or matrix of estimated parameters as roots of the quasi-score vector
#' @param opts		list of upper bounds for a root of the quasi-score vector, see details
#' @param verbose   logical, \code{TRUE} for intermediate output
#' 
#' @return A data frame with columns named corresponding to each component of the investigated parameter,
#'  `\code{minor}`, `\code{det}`, `\code{value}`, `\code{score_max}`, `\code{logdetI}`, `\code{lamI_min}`,
#'  `\code{lamIm_max}` and `\code{varS_trace}` (see vignette for details). The first column shows the leading minor of
#'   the observed QI matrix which is not positive definite. Estimated model parameters for which the observed QI matrix
#'   is not positive definite are excluded from further root selection. 
#' 
#' @details Only for the quasi-likelihood approach the function inspects the (numerical) consistency of the found
#'  parameters in `\code{par}` by comparing each observed quasi-information matrix with the expected one.
#'  The degree of dissimilarity of both matrices is measured by certain scalar equivalent criteria (see vignette)
#'  and the parameter for which these are smallest is chosen. 
#'  
#' @examples 
#'  data(qleresult)
#' 
#'  # and just a single root 
#'  checkMultRoot(OPT,verbose=TRUE)
#' 
#' @author M. Baaske
#' @rdname checkMultRoot
#' @export
checkMultRoot <- function(est, par = NULL, opts = NULL, verbose = FALSE)
{			
   if(est$qsd$criterion != "qle")
	  stop("Consistency check of multiple roots only for criterion `qle`.")
   if(.isError(est))
	  stop("The estimation result from function `qle` has errors. Please check the argument `est`.")
   options <- list("ftol_abs"=1e-6, "score_tol"=1e-4)
   
   opts <-
    if(is.null(opts))
	 options	
    else {		 
	  # check defaults	  
	  id <- which(!is.na(pmatch(names(opts),names(options))))
	  if (length(id)>0L)
		options[names(opts[id])] <- opts[id]	
	  options
   }  
   # always use estimate from est first
   if(!is.null(par)){
	   par <- .LIST2ROW(par)
	   stopifnot(NCOL(par)==length(est$par))	   
   }
   par <- rbind("par"=est$par,par) 
   info <- attr(est,"optInfo")  			   
   QD <- try(quasiDeviance(par,est$qsd,W=info$W,
			  theta=info$theta,cvm=est$cvm,verbose=verbose),
		  silent=TRUE)
   
   if(.isError(QD)) {
	   msg <- paste0("Failed to get quasi-deviance: ","\n")
	   message(msg)
	   return(.qleError(message=msg,call=match.call(),error=QD))
   }
   dm <- try(.evalRoots(QD,par,opts),silent=TRUE)	
   if(.isError(dm)) {
	   msg <- .makeMessage("Could not check consistency of roots.")
	   message(msg)
	   return(.qleError(message=msg,call=match.call(),error=dm))	   
   }
   return (dm)  
}

# intern use only:'value' is modified QD/Mahalanobis distance as test statistic
.rootTest <- function(par, value, I, obs, alpha, test, ...,
		       			multi.start = 0L, Npoints = 10, cl = NULL,
			   			  na.rm = TRUE, pl = 0L, verbose = FALSE){	
	aiqm <- NULL
	mScore <- NULL	
	xdim <- length(par)
	opt.args <- list(...)
	hasError <- integer(0)	
	stopifnot(is.numeric(value))
	
	RES <- 
	 if(multi.start > 0L){
		if(verbose)
		   cat("Re-estimate parameters (possibly use multi-start approach):","\n")
		# multi start root finding, no nested parallel execution
		# including restart if more than one method is given 
		if(is.null(opt.args$nstart))
			opt.args$nstart <- (xdim+1L)*Npoints
		do.call(doInParallel,
				c(list(X=obs,
					FUN=function(obs,...) {
						# not in parallel!
						multiSearch(x0=par,...,obs=obs,inverted=TRUE,check=FALSE,
							multi.start=(multi.start > 1L),cl=NULL,verbose=FALSE,cores=1L)
					},
					cl=cl), opt.args))
	 } else {
		#including restart if more than one method is given
		if(verbose)
		  cat("Re-estimate parameters:","\n")
		do.call(doInParallel,
				c(list(X=obs,
						FUN=function(obs,...) {
							searchMinimizer(x0=par,...,obs=obs,
								inverted=TRUE,check=FALSE,verbose=verbose)
						},
					cl=cl), opt.args))   
	 }

	if(.isError(RES))
	  return(RES)
	
    # check results again
	ok <- which(sapply(RES,function(x) !.isError(x) && x$convergence >= 0L))
	if(length(ok) == 0L){
		stop(.makeMessage("All re-estimations failed or did not converge.","\n"))						
	} else if(length(ok) < length(RES)){
		message(paste0("A total of ",length(RES)-length(ok)," re-estimations failed or did not converge."))							
	}	
	# average inverse QI
	invI <- 
		lapply(RES[ok],
			function(x) {						
				try(gsiSolve(x$I),silent=TRUE) 
			})
	badInv <- sapply(invI,function(x) inherits(x,"try-error") || anyNA(x))
	if(any(badInv))
	 message(paste0("A total of ",sum(badInv)," inversions of quasi-information matrices failed. Check attribute `info`."))    
	else if(all(badInv)){
	 msg <- paste0("All iversions of quasi-information matrices failed.")
	 message(msg)
	 return(.qleError(message=msg,call=match.call(),error=badInv))	 
	}
	# average matrix of inverse qi matrices
	aiqm <- matrix(
			colMeans(
				do.call(rbind,
					lapply(invI[!badInv],as.numeric)
				)
			),ncol=xdim)				
	
	# estimates	
	mpars <- do.call(rbind,lapply(RES[ok],function(x) x$par))
	mScore <- do.call(rbind,lapply(RES[ok],function(x) x$score))
	has.na <- (rowSums(is.na(cbind(mScore,mpars))) > 0L)	
	
	if(na.rm && any(has.na)) {
		ok <- ok[-which(has.na)]
		mpars <- mpars[ok,,drop=FALSE]
		mScore <- mScore[ok,,drop=FALSE]
		warning("Removed `NA` values from quasi-scores.")
	}
	# average quasi-score
	mScore <- try(colMeans(mScore),silent=TRUE)
	# some (empirical) measures	
	msem <- .MSE(mpars,par)	
	# value of test statistic at re-estimated parameters			
	tvals <- sapply(RES[ok],"[[","value") 
	stopifnot(is.numeric(tvals))
	
	# invert QI for predicted std. error (asymptotic) at estimated theta 
	qi <- try(gsiInv(I),silent=TRUE)
	if(inherits(qi,"try-error") || anyNA(qi))
		message("Inversion of quasi-information matrix failed")
	
	# get efficient score test (with MC parameters)
	B <- structure(
			data.frame(
					cbind("par"=par,
						  "se"=apply(mpars,2,sd),					
						  "rmse"=sqrt(diag(msem)),
						  "bias"=colMeans(t(t(mpars)-par)),
						  "mean"=colMeans(mpars))),
			"sb"=value, "Sb"=tvals,
			"test"=test)
	
	relEF <-
	 if(!anyNA(c(msem,qi)) && is.matrix(qi) && is.matrix(msem)) {
		 #try(abs(1-sqrt(diag(msem))/sqrt(diag(qi))),silent=TRUE)
		  try(abs(1 - sqrt(diag(qi))/sqrt(diag(msem))),silent=TRUE)
	} else {
		message("Failed to compute relative difference of empirical and predicted error.")
		NULL
	}
	
	# had errors
	hasError <- which(!(1:length(RES) %in% ok))
	if(length(hasError) > 0L)	
		message(paste0("A total of ",length(hasError)," re-estimations failed."))
		
	res <- .qleTest(B,alpha)					# test results
	res$par <- par
	
	# results
	structure(res,				
			msem=msem,							# mean square error matrix
			aiqm=aiqm,							# average inverse QI (re-estimated parameters)
			qi=qi,								# inverse QI at estimated theta
			relEF=relEF,
			obs=NULL,							# (MC) observations
			mean.score=mScore,					# average score/gradient
			criterion=test,
		info=list(badInv=which(badInv), 			# inversion errors
					hasNa=which(has.na), 			# indices of NA parameters 
					hasError=hasError,
					iseed=NULL),																
	class=c("qleTest"))	
}

#' @name 	qleTest
#'
#' @title	Monte Carlo testing
#'
#' @description Monte Carlo hypothesis testing 
#'
#' @param est			object of class \code{qle} or \code{mahal}, estimation results after calling function \code{\link{qle}}
#' @param par0			optional, vector of parameter for the null hypothesis 
#' @param obs0			optional, vector of observed statistics corresponding to `\code{par0}`
#' @param ...			arguments passed to the simulation function `\code{sim}`, \code{\link{searchMinimizer}} and \code{\link{multiSearch}}
#' @param sim			user supplied simulation function, see \code{\link{qle}}
#' @param criterion		optional, \code{NULL} (default), name of the test statistic, either "\code{qle}" or "\code{mahal}" which overwrites the function criterion used for estimation of the model parameter
#' @param nsim			numeric, number of (initial) simulation replications for each new sample point
#' @param fnsim 		optional, a call returning the number of simulation replications applied to a new
#' 						sample point with the current environment of calling function \code{qle},
#' 						default is the initial value `\code{qsd$nsim}`, respectively `\code{nsim}`
#' @param obs			optional, \code{NULL} (default), simulated statistics at the hypothesised parameter, if not given, these are generated at `\code{par0}` or at `\code{est$par}` 
#' @param alpha			significance level for testing the hypothesis
#' @param multi.start   integer, \code{=0,1,2}, level of multi start root finding (see details)
#' @param na.rm 		logical, \code{TRUE}  (default), whether to remove `NA` values from the matrix of re-estimated parameters
#' @param cores 		number of cores for multistart searches for each given/generated observation, only if \code{multi.start>0} enabled and ignored otherwise
#' @param cl			cluster object, \code{NULL} (default), of class \code{MPIcluster}, \code{SOCKcluster}, \code{cluster}
#' @param iseed			integer, the seed for initializing the cluster workers for parallel computations 
#' @param verbose   	logical, \code{TRUE} for intermediate output
#'
#' @details The function tests the null hypothesis \eqn{H_0:\,\hat{\theta}=\theta_0}, that is, whether the statistical
#'  model w.r.t. to the estimated parameter explains the observed statistics, against the alternative \eqn{H_1:\,\hat{\theta}\neq\theta_0} based
#'  on a Monte Carlo approach (see vignette). Due to the approximate nature of the assumed statistical model for the observed statistics the
#'  exact distribution of the test statistics, that is, the Mahalanobis distance or quasi-deviance, is generally unknown and therefore
#'  its asymptotic distribution might be an unrealistic assumption for the null hypothesis. For this reason, and in order to retrieve an empirical
#'  p-value for testing, we generate bootstrap observations and re-estimate the model parameter for each observation in the same way as done before
#'  when estimating the model parameter. This includes all possible types of variance approximations available (by kriging or average approximations)
#'  and types of prediction variances (by kriging or cross-validation).
#' 
#'  The function expects an estimation result `\code{est}` as returned from the main estimation function \code{\link{qle}}. If any simulated statistics
#'  are available at the final parameter estimate or at `\code{par0}`, then these can be passed by `\code{obs}` and used as bootstrapped observations of the
#'  summary statistics. Otherwise the function first generates those using `\code{nsim}` model replications. The criterion function approximations are used as
#'  specified in the object `\code{qsd}` and will not be further augmented by additional samples or simulations during the test procedure.
#'  The value of the test statistic is either chosen as the current criterion function value at the estimated parameter or it is re-computed at the
#'  given parameter `\code{par0}` using, if given, the `real` observed statistics `\code{obs0}`. The user can also select a different criterion function
#'  as a test statistic compared to the estimation before which can be set by `\code{criterion}`. Apart from the quasi-deviance as a test statistic, in
#'  principle, any supported type of a least squares criterion, more generally, the Mahalanobis distance, can be used which only depends on the prefered type
#'  of variance matrix approximation, see \code{\link{covarTx}}.
#' 
#'  In order to efficiently find the roots of the quasi-score vector we implement a multi start concept for minimizing the criterion function.
#'  Option `\code{multi.start=0}` starts a single root finding from the estimated parameter (as a starting point) for each newly generated observation.
#'  Using `\code{multi.start=1}` starts a multi start root finding only in case the local optimization gets stuck into a local minimum or does not
#'  converge and setting `\code{multi.start=2}` always triggers a multi start local search for each simulated observation. Practically, the re-estimations
#'  of the parameters might still fail to converge. However, the user can control the convergence conditions of the local solvers
#'  (including the quasi-scoring iteration) by the corresponding control parameters (see \code{\link{searchMinimizer}}). Any failed re-estimation is
#'  excluded from the test results and stored in the attribute `\code{info}`. In addition, as part of the returned data frame `\code{param}`
#'  the empirical standard error, predicted standard error (based on the average inverse quasi-information matrix), the root mean square error,
#'  the bias and sample mean value of the re-estimated parameters are also available. For a full example we refer the reader to the package vignette.
#' 
#' @return An object of class \code{qleTest} as a list of:
#'  \item{param}{ data frame of estimated parameters and error measures}
#' 	\item{test}{ the test result}
#'  \item{Stest}{ name of the test} 
#' 
#' with attributes:
#' 
#' 	 \item{msem}{ mean square error matrix of re-estimated parameters}
#'   \item{aiqm}{ average inverse quasi-information matrix over all re-estimated parameters}
#' 	 \item{qi}{ inverse quasi-information matrix at the parameter to be tested `\code{est$par}`}
#'   \item{relEF}{ relative difference of the empirial and predicted standard error of the parameter to be tested} 
#'   \item{obs}{ list of simulated statistics either at the estimated parameter or at the optional parameter `\code{par0}`}
#'   \item{optRes}{ results from re-estimating the model parameters for each simulated observation `\code{obs}`}
#'	 \item{mean.score}{ average quasi-score, respectively, average gradient of the MD at the re-estimated parameters}
#' 	 \item{criterion}{ name of criterion function used as a test statistic: "\code{qle}" or "\code{mahal}"}  
#' 	 \item{info}{ list of the following elements: indices of re-estimation results where the inversion of the quasi-information matrix failed,
#'       the re-estimated parameters have `NA`s, criterion function minimizations failed or did not converge numerically,
#'       the integer seed value `\code{iseed}`}
#' 
#' @author M. Baaske
#' @rdname qleTest
#' @export
qleTest <- function(est, par0 = NULL, obs0=NULL, ..., sim, criterion = "qle",
		             nsim = 100, fnsim = NULL, obs = NULL, alpha = 0.05, multi.start = 0L,
					  na.rm = TRUE, cores = 1L, cl = NULL, iseed = NULL, verbose = FALSE)
{				  
	if(.isError(est))
	  stop("Estimation result has errors. Please see attribute `error`.")    
    # last evaluation of criterion function  	
	if(.isError(est$final))
	  stop("Final criterion function evaluation failed. Please check attribute `error`.")
	# simulation increase function `nsim` 
	if(missing(nsim))
	  nsim <- attr(est$qsd$qldata,"nsim")  	  
    Fnsim <-
	  if(is.null(fnsim)){
		  as.call(list(function(n) n, quote(nsim)))
	  } else if(is.call(fnsim)) {		 
		  fnsim[[1]] <- match.fun(fnsim[[1]])		 
		  # make current environment available in function call 
		  environment(fnsim[[1]]) <- environment()		 
		  fnsim
	  }
	  else {
		  stop("Expected numeric value or an object of class `call` in argument `nsim`.")
	 }
    args <- list(...)	
	# basic checks
	stopifnot(class(est) == "qle")
  	stopifnot(class(est$qsd)=="QLmodel")
	Npoints <- nrow(est$qsd$qldata)							# number of multi-start points
	xdim <- attr(est$qsd$qldata,"xdim")						# dimension of the parameter
	criterion <- match.arg(criterion,c("qle","mahal"))  	# test statistic 
			
	# check arguments for local searches
	id <- pmatch(names(args),names(formals(searchMinimizer)))
	# check `searchMinimizer' args 
	opt.args <- args[which(!is.na(id))]	
	
	# check with qsd hidden	
	.checkfun(searchMinimizer,opt.args,
			hide.args=c("x0","qsd"),check.default=FALSE)
	
	# use last Sigma (unless kriging variance matrix)  
	info <- attr(est,"optInfo")
	if(est$qsd$var.type != "kriging")		
	  opt.args <- c(opt.args,list(W=info$W,theta=info$theta))
     
	# test at estimated parameter: simply overwrite `est$final`
	# with something of class QSResult for testing another parameter
	if(is.null(par0)){		
		local <- est$final
		if(.isError(local) || !attr(est,"optInfo")$minimized)
		 stop("Final optimization failed. Please check attribute `final` and `optInfo`.")
	 	else if(local$convergence < 0)
		  warning(paste0("Last local search did not converge by method: ",local$method)) 
	  	# use estimated parameter as default:
		# which might render the test meaningless unless
		# a different test statistic, i.e. "mahal" or "qle" (set by `criterion`),
		# is used as before for estimation!
	  	par0 <- local$par
    }
		
	if(est$qsd$criterion != criterion) {	
		est$qsd$criterion <- criterion
		# check input observed statistics
		# default: use original data (statistics)
	    if(!is.null(obs0)){
		  obs0 <- unlist(obs0)
		  if(anyNA(obs0) | any(!is.finite(obs0)))
			  warning("`NA` or `Inf` values detected in `obs0`.")
		  if(!is.numeric(obs0) || length(obs0) != length(est$qsd$covT))
			  stop("`obs0` must be a (named) `numeric` vector or list of length equal to the number of caoariance models `qsd`.")
		  # overwrite observed statistics
		  est$qsd$obs <- obs0
	    }		
		
        local <-
		 tryCatch({				
			 if(est$qsd$criterion == "qle"){
				 quasiDeviance(par0,est$qsd,Sigma=attr(est$final,"Sigma"),
						 W=info$W,theta=info$theta,cvm=est$cvm,verbose=verbose)				 
			 } else {
				 mahalDist(par0,est$qsd,Sigma=attr(est$final,"Sigma"),
						 W=info$W,theta=info$theta,inverted=TRUE,cvm=est$cvm,
						  verbose=verbose)
			 }
		 }, error = function(e) {
			  msg <- .makeMessage("Error in criterion function evaluation due to ",conditionMessage(e))				 
			 .qleError(message=msg,call=sys.call(),error=e)		
		 })
		 if(!.isError(local)){			
			 local <- local[[1]]
		 } else {
			 message(paste0("Cannot continue testing due to ",local$message))
			 return(local)
		 }			
	}
		
    # MC simulation of observed statistics
	# if no observations supplied	
	if(is.null(obs)){		
		nid <- which(!is.na(id))
		if(length(nid) > 0L)
		 args <- args[-nid]		
		sim <- match.fun(sim)
		# check `sim` input values
		.checkfun(sim,args,remove=TRUE)
		nsim <- try(eval(Fnsim,envir=environment()),silent=TRUE)
		stopifnot(is.numeric(nsim) || nsim > 0)
				
		simFun <- function(x) try(do.call(sim,c(list(x),args)))		
		sim.args <-
			list(sim=simFun,X=par0,nsim=nsim, mode="list",
					cl=cl,iseed=iseed,verbose=verbose)
		if(verbose)
		 cat("Simulate observed statistics...","\n")
 		
		obs <- tryCatch(
				do.call(simQLdata,sim.args),
				error = function(e) {
					msg <- paste0("Simulating observed statistics failed: ",conditionMessage(e))
					message(msg) 
					.qleError(message=msg,call=match.call(),error=e)		   
				}
			)
		if(.isError(obs))			
		 return(obs)
	 	
	} else {
		if(.isError(obs))
		  stop("Argument `obs` has errors.")
		else if(class(obs) != "simQL" || !is.list(obs))
		  stop("Argument `obs` must be of class `simQL` and `list` type.")	   
	}
		
 	RES <- 
 	 if(multi.start > 0L){																		# multi-start only in case of non-convergence  
		if(verbose)
		 cat("Re-estimate parameters (possibly use multi-start approach):","\n")
		# multi start root finding, no nested parallel execution
		# including restart if more than one method is given 
		if(is.null(opt.args$nstart))
		  opt.args$nstart <- (xdim+1L)*Npoints
		do.call(doInParallel,
			c(list(X=obs[[1]],
					FUN=function(obs,...) {
						# no cluster but local parallel execution
						multiSearch(x0=par0,qsd=est$qsd,...,   		
						 cvm=est$cvm,obs=obs,inverted=TRUE,check=FALSE,
						   multi.start=(multi.start>1L),cl=NULL,						# multi.start = 2L, then always multi-start
						   cores=cores,verbose=FALSE) 
					},
					cl=cl), opt.args))
	} else {																			# never multi-start but restart if provided another routine
		if(verbose)
		 cat("Re-estimate parameters:","\n")
		#including restart if more than one method is given
		do.call(doInParallel,
			  c(list(X=obs[[1]],
					FUN=function(obs,...) {
						searchMinimizer(x0=par0,qsd=est$qsd,...,    
							cvm=est$cvm,obs=obs,inverted=TRUE,
							  check=FALSE,verbose=verbose)
				},
				cl=cl), opt.args))   
	}
    if(.isError(RES)) {
	  msg <- paste0("Could not find MC replicated parameters.")
	  message(msg)
	  return(.qleError(message=msg,call=match.call(),error=RES))	
	}
	# check results again
	ok <- which(sapply(RES,function(x) !.isError(x) && x$convergence >= 0L ))
	if(length(ok) == 0L){
		msg <- paste("All re-estimations failed or did not converge: ")
		warning(msg)
		return(.qleError(message=msg,call=match.call(),error=RES))	   
	} else if(length(ok) < length(RES))
		message("A total of ",length(ok), " failures in re-estimating the parameters. Check attribute `optRes` and `info`.")
  
	# also check H_0?
	# because sampling MC theta
	# must be done under H0	
	aiqm <- NULL
	mScore <- NULL	    	  
	# average quasi-information matrix
	invI <- 
		lapply(RES[ok],
			function(x) {						
				try(gsiInv(x$I),silent=TRUE) 
			})
	badInv <- sapply(invI,function(x) inherits(x,"try-error") || anyNA(x))
	if(any(badInv))
		message(paste0("A total of ",sum(badInv)," inversions of quasi-information matrices failed. Check attribute `info`."))    
	else if(all(badInv)){
		msg <- paste0("All iversions of quasi-information matrices failed.")
		message(msg)
		return(.qleError(message=msg,call=match.call(),error=badInv))	 
	}
	# average matrix of inverse qi matrices
	aiqm <- matrix(
			 colMeans(
				do.call(rbind,
					lapply(invI[!badInv],as.numeric)
				)
			),ncol=xdim)				
 	
	# estimates	
	mpars <- do.call(rbind,lapply(RES[ok],function(x) x$par))
	mScore <- do.call(rbind,lapply(RES[ok],function(x) x$score))	
	has.na <- (rowSums(is.na(cbind(mScore,mpars))) > 0L)	
	if(na.rm && any(has.na)) {
		ok <- ok[-which(has.na)]
		mpars <- mpars[ok,,drop=FALSE]
		mScore <- mScore[ok,,drop=FALSE]
		warning("Removed `NA` values from quasi-scores.")
  	}
	mScore <- try(colMeans(mScore),silent=TRUE)		
	if(nrow(mpars) <= 10L)
	 warning("Only a small number of 10 or less parameters could be re-estimated.")	
    
 	# some (empirical) measures	
	msem <- .MSE(mpars,local$par)
		
	# value of test statistic at re-estimated parameters			
	tvals <- sapply(RES[ok],"[[","value")	
	stopifnot(is.numeric(c(local$value,tvals)))
	
	# invert QI for predicted std. error (asymptotic) at estimated theta 
	qi <- try(gsiInv(local$I),silent=TRUE)
	if(inherits(qi,"try-error") || anyNA(qi))
	  message("Inversion of quasi-information matrix failed")
	
	# get efficient score test (with MC parameters)
	B <- structure(
			 data.frame(
			  cbind("par"=local$par,
				    "se"=apply(mpars,2,sd),					
					"rmse"=sqrt(diag(msem)),
					"bias"=colMeans(t(t(mpars)-local$par)),
					"mean"=colMeans(mpars))),
			"sb"=local$value, "Sb"=tvals,
			"test"=est$qsd$criterion)
	
	# had errors
	hasError <- which(!(1:length(RES) %in% ok))
	if(length(hasError) > 0L)	
	  message(paste0("A total of ",length(hasError)," re-estimations failed."))
  
	relEF <-
	 if(!anyNA(c(msem,qi)) && is.matrix(qi) && is.matrix(msem)) {
		 #try(abs(1-sqrt(diag(msem))/sqrt(diag(qi))),silent=TRUE)
		 try(abs(1 - sqrt(diag(qi))/sqrt(diag(msem))),silent=TRUE)
	 } else {
		message(.makeMessage("Detected `NAs` values (for relative differences) in quasi-information or MSE matrix while testing."))		
		NA
	 }
	res <- .qleTest(B,alpha)					# test results
	res$par <- est$par
	
	# results
	structure(res,				
		    msem=msem,							# mean square error matrix
	    	aiqm=aiqm,							# average inverse QI (re-estimated parameters)
			qi=qi,								# inverse QI at estimated theta
	    	relEF=relEF,
			obs=if(verbose) obs else NULL,		# (MC) observations
	    	optRes=if(verbose) RES else NULL,	# optimization results
			mean.score=mScore,					# average score/gradient
			criterion=est$qsd$criterion,			
	  info=list(badInv=which(badInv), 			# inversion errors
			  	hasNa=which(has.na), 			# indices of NA parameters 
		 		hasError=hasError,
				iseed=iseed),																
	 class=c("qleTest"), call=match.call())	
}

# printing function
#' @name print.qleTest
#' 
#' @title print \code{qleTest} results 
#' 
#' @description print the results from \code{\link{qleTest}}
#' 
#' @param x      object of class \code{qleTest} from \code{\link{qleTest}}
#' @param pl	 ignored
#' @param digits number of (significant) digits
#' @param format format character(s), see \code{\link{formatC}}
#' @param ... 	 ignored
#' 
#' @rdname print.qleTest
#' @method print qleTest
#' @export 
print.qleTest <- function(x, pl = 1, digits = 5, format="e", ...) {
	if(.isError(x)){
	 	print(x)
		invisible(return(NULL))
 	}
	
	if(!is.null(attr(x,"call"))){
	  cat("\nCall:\n\n")
	  cat(paste(deparse(attr(x,"call")), sep="\n", collapse = "\n"), "\n\n", sep="")
	}
	# consistency check of solution
	chk <- attr(x,"solInfo")
	if(!is.null(chk)){
		cat("Consistency check - the smaller the better: \n\n")
		print(format(signif(chk,digits=digits)),print.gap=2,right=TRUE,quote=FALSE)
		cat("\n\n")		
	}	
	cat("Coefficients:\n\n")	
	print(format(x$param, digits=digits), print.gap = 2, right=TRUE, quote = FALSE)	
	cat("\n\n")
	cat(x$Stest,"\n\n")
	vals <- c(format(x$test[1], digits=digits),
			  formatC(signif(x$test[2], digits=digits), digits=digits, format=format, flag="#"))
	names(vals) <- colnames(x$test)
	print(vals, print.gap = 2, right=TRUE, quote = FALSE)
	
	cat("\n\n")
	if(!is.null(attr(x,"mean.score"))) {
	  if(attr(x,"criterion") == "mahal")
		cat("Average gradient: \n\n")
	  else cat("Average quasi-score: \n\n")
	  print(format(attr(x,"mean.score"), digits=digits), print.gap = 2, right=TRUE, quote = FALSE)
	  cat("\n\n")
	  qi <- attr(x,"qi")
	  aiqm <- attr(x,"aiqm")
	  if(!is.null(aiqm) && !.isError(aiqm) &&
		 !is.null(qi) && !.isError(qi) &&
		 !is.null(attr(x,"relEF"))) {
	 
		  pse <- as.data.frame(cbind(
			      formatC(signif(sqrt(diag(aiqm)),digits=digits), digits=digits, format=format, flag="#"),
				  formatC(signif(sqrt(diag(qi)),  digits=digits), digits=digits, format=format, flag="#"),
				  formatC(signif(attr(x,"relEF"), digits=digits), digits=digits, format=format, flag="#")))
	  	  dimnames(pse) <- list(row.names(x$param),c("Average", "Estimate", "EF"))
		  cat("Predicted Std. Errors: \n\n")
		  print(format(pse, digits=digits), print.gap = 2, right=TRUE, quote = FALSE)
	  } else {
		  message("Cannot show some of the error matrices (see results).")
	  }
    }	
	invisible(x)	
}
