# Copyright (C) 2018 Markus Baaske. All Rights Reserved.
# This code is published under the GPL (>=3).
#
# File: 	krige.R
# Date:  	14/03/2018
# Author: 	Markus Baaske
#
# Implements kriging prediction, variance approximation methods,
# quasi-deviance and Mahalanobis distance function as a criterion
# for QL estimation

#' @name estim
#'
#' @title Kriging prediction and numerical approximation of derivatives
#'
#' @description 	
#'
#' @param models 	object of class \code{krige} either as a list of covariance models or
#' 	 				class \code{covModel} as a single covariance model, see \code{\link{setCovModel}}
#' @param points 	matrix or list of points to predict the sample means of statistics
#' @param Xs		matrix of sample points
#' @param data		data frame of sample means of statistics at sampled points
#' @param krig.type name of kriging type, either "\code{dual}" (default) or "\code{var}"
#'
#' @return
#'  \item{estim}{ list of predicted values of sample means of statistics (including prediction
#' 		 variances if `\code{krig.type}` equals to "\code{var}")}
#'  \item{jacobian}{ list of Jacobians at predicted values of sample means of statistics}  
#'
#' @details The function can be used to predict any values by kriging given a covariance model. Each covariance model is given as an
#'  element of the list `\code{models}` including its own trend model and covariance function name. There are two types of kriging predictors
#'  available. First, the \emph{dual kriging} predictor, set by `\code{krig.type}`="\code{dual}" or the one based on the calculation of prediction
#'  variances, setting `\code{krig.type}` to "\code{var}". Both types result in exactly the same predicted values and only differ by whether or not
#'  kriging variances should be computed. The data, e. g. sample means of each statistic, must be given as column vectors where each row corresponds
#'  to a sample point in the data frame given in the argument `\code{data}`. 
#'
#' @examples 
#' data(normal) 
#' 
#' X <- as.matrix(qsd$qldata[,1:2])
#' p <- c("mu"=2,"sd"=1)
#' 
#' # get simulated statistics at design X
#' Tstat <- qsd$qldata[grep("^mean[.]",names(qsd$qldata))]
#' 
#' # low level prediction, variances and weights
#' estim(qsd$covT,p,X,Tstat,krig.type="var")
#' 
#' # Jacobian 
#' jacobian(qsd$covT,p,X,Tstat)
#'  
#'    
#' @author M. Baaske
#' @rdname estim
#' @export
estim <- function(models, points, Xs, data,krig.type=c("dual","var")) {   
	UseMethod("estim",models)
}

#' @method estim krige 
#' @export
estim.krige <- function(models, points, Xs, data, krig.type = c("dual","var")) { 
	krig.type <- match.arg(krig.type)		
	if(!is.list(data) || length(data)!=length(models))
		stop("Expected 'data' as  a list of same length as 'models'.")	 		 	 
	
	if(!is.matrix(points))
		points <- .LIST2ROW(points)
	
	.Call(C_kriging,Xs,data,points,models,krig.type)
}

#' @method estim covModel
#' @export 
estim.covModel <- function(models, points, Xs, data,krig.type=c("dual","var")) {	
	krig.type <- match.arg(krig.type)		
	if(!is.matrix(points))
		points <- .LIST2ROW(points)
	
	.Call(C_kriging, Xs, as.data.frame(data), points, list(models), krig.type)
}

#' @title Jacobian
#'
#' @description Jacobian of mean values of statistics 
#' 
#' @inheritParams estim
#'
#' @details
#'   The function `\code{jacobian}` computes the partial derivatives of sample means of the statistics
#'   as columns and for each component of the parameter vector as rows by forward difference approximations.
#' 
#' @rdname estim
#' 
#' @export
jacobian <- function(models, points, Xs, data, krig.type=c("dual","var")) {
	UseMethod("jacobian",models)
}

#' @method jacobian krige 
#' @export
jacobian.krige <- function(models, points, Xs, data,krig.type=c("dual","var")) {
	krig.type <- match.arg(krig.type)	
	if(!is.list(data) || length(data)!=length(models))
		stop("Expected \'data\' as  a list of same length as 'models'.")	
	if(!is.list(points))
		points <- .ROW2LIST(points)	
	.Call(C_estimateJacobian,Xs,data,points,models,krig.type)
}

#' @method jacobian covModel 
#' @export
jacobian.covModel <- function(models, points, Xs, data,krig.type=c("dual","var")) {
	krig.type <- match.arg(krig.type)	
	if(!is.list(points))
	  points <- .ROW2LIST(points)
  
	.Call(C_estimateJacobian,Xs,as.data.frame(data),points,list(models),krig.type)
}

#' @title Kriging the sample means of statistics
#'
#' @description
#'  \describe{
#'    \item{\code{predictKM},}{wrapper for kriging the sample means of statistics} 	  
#'	  \item{\code{varKM},}{ calculate the kriging prediction variances}
#'	  \item{\code{extract},}{ extract the results of kriging} 
#'  } 
#' 
#' @details For a list of fitted covariance models the function \emph{predictKM} predicts the
#' 	 sample means of statistics at (unsampled) points, calculates the prediction variances, if applicable,
#'   at these points and extracts the results. Note that, since we aim on predicting the "error free" value of the sample means,
#'   we use the \emph{smoothing} kriging predictor as described in [2, Sec. 3.7.1].
#'
#' @param models   list of covariance models, see \code{\link{setCovModel}}
#' @param ... 	   further arguments passed to function \code{\link{estim}}
#' 
#' @return \item{predictKM}{ list of kriging predicted values}
#' 
#' @examples 
#' data(normal)
#' X <- as.matrix(qsd$qldata[,1:2])
#' p <- c("mu"=2,"sd"=1)
#' 
#' # get simulated statistics at design X
#' Tstat <- qsd$qldata[grep("^mean[.]",names(qsd$qldata))]
#' 
#' # predict and extract 
#' predictKM(qsd$covT,p,X,Tstat)
#' 
#' # prediction variances
#' varKM(qsd$covT,p,X,Tstat)
#' 
#' @rdname krige
#' @export
predictKM <- function(models,...) {
	km <- estim(models,...)
	if(.isError(km)) {
	  message("Kriging estimation (mean) failed.\n")
	  return(km)
  	}
	extract(km, type="mean")
}

#' @return \item{varKM}{ list of kriging prediction variances}
#' 
#' @rdname krige
#' @export
varKM  <- function(models, ...) {
	km <- estim(models, ..., krig.type="var")
	if(.isError(km)) {
	  message("Kriging prediction error estimation (var) failed.\n")
	  return(km)
	}
 	extract(km, type="sigma2")
}


#' @param X	 	kriging result
#' @param type  return type of results, see details below
#'
#' @details
#'  The function \emph{extract} either returns the predicted values, the
#'  prediction variances or the kriging weights for each point. 
#'
#' @return \item{\code{extract}}{ matrix of corresponding values (see details)}  
#'
#' @rdname krige
#' @export
extract <- function(X, type = c("mean","sigma2","weights")) UseMethod("extract",X)

#' @method extract krigResult 
#' @export
extract.krigResult <- function(X, type=c("mean","sigma2","weights")) {
	type <- match.arg(type)
	switch(type,
	   "weights" = { 
		   if(!is.null(attr(X,"weights")))		# only for dual kriging weights!
			 return( t(do.call(rbind,(attr(X,"weights")))) )
		   lapply(X,function(x) {
			 w <- matrix(unlist(x$weights),nrow=length(x$weights[[1]]))
			 colnames(w) <- names(x$weights)
			 w
		 	})	
	   },
	   do.call(rbind,lapply(X,"[[",type))
	)
}

# intern
#' @importFrom expm logm
varLOGdecomp <- function(L) {
	vmats <- try(lapply(1L:nrow(L), function(i) .chol2var(unlist(L[i,]))), silent=TRUE)
	if(inherits(vmats,"try-error")) {
	  return(simpleError(.makeMessage("Matrix logarithm of covariance matrices failed.")))
	}
	decomp <- lapply(vmats,
			   function(Xs) {
				   m <- try(expm::logm(Xs,method="Eigen"),silent=TRUE)
				   if(inherits(m,"try-error"))
					return (m)
				   as.vector(m[col(Xs)>=row(Xs)])
				 }
		       )
   err <- sapply(decomp, function(x) .isError(x))
   if(any(err)) {
	   msg <- paste0("Matrix logarithm decomposition failed: ")
	   message(msg)
	   return(.qleError(n=msg,call=match.call(), hasError=which(err), error=decomp))
   }
   as.data.frame(unlist(do.call(rbind,decomp)), ncol=length(decomp[[1L]] ) )
}

#' @title Variance matrix approximation
#'
#' @description Approximating the variance-covariance matrix of statistics
#'
#' @param qsd		object of class \code{\link{QLmodel}}
#' @param W			weight matrix for weighted average approximation
#' @param theta	    parameter vector for weighted average approximation 
#' @param cvm		list of fitted cross-validation models, see \code{\link{prefitCV}}
#' @param doInvert	if \code{TRUE}, return the inverse of the variance matrix approximation
#'
#' @return
#' 	List of variance matrices with the following structure:
#' 	\item{VTX}{ variance matrix approximation}
#'  \item{sig2}{kriging prediction variances of statistics at \code{theta}}
#'  \item{inv}{ if applicable, the inverse of either \code{VTX}}
#'
#' @details	The function estimates the variance matrix of statistics at some (unsampled) point \code{theta} by either
#'  averaging (the \emph{Cholesky} decomposed terms or matrix logarithms) over all simulated variance matrices
#'  of statistics at previously evaluated points of the parameter space or by a kriging approach which treats the Cholesky
#'  decomposed terms of each variance matrix as a data vector for kriging.
#' 
#'  In addition, a Nadaraya-Watson kernel-weighted average approximation can also be applied in order to bias the variance
#'  estimation towards a more locally weighted estimation, where smaller weights are assigned to points being more
#'  distant to an estimate of the model parameter \code{theta}. A reasonable symmetric weighting matrix 
#'  \code{W} of size equal to the problem dimension, say \code{q}, can be freely chosen by the user. In addition, the user can select
#'  different types of weighted variance averaging methods such as "\code{wcholMean}", "\code{wlogMean}",
#'  "\code{kriging}" (kriging the variance matrix of statitics) defined by \code{qsd$var.type}. If \code{theta} is given, then prediction variances
#'  of the statistics are either derived from the kriging results of the statistics at \code{theta} or based on a (possibly more robust)
#'  cross-validation (CV) approach, see the vignette for details.
#'    
#' @examples 
#'  data(normal)
#'  # average approximation of variance matrices
#'  covarTx(qsd,theta=c("mu"=2,"sd"=1))
#' 
#' @author M. Baaske
#' @rdname covarTx
#' @export
covarTx <- function(qsd, W = NULL, theta = NULL, cvm = NULL, doInvert = FALSE)
{		
	xdim <- attr(qsd$qldata,"xdim")
	Xs <- as.matrix(qsd$qldata[seq(xdim)])
	
	nstat <- length(qsd$covT)
   	var.type <- qsd$var.type
	krig.type <- qsd$krig.type
	Tnames <- names(qsd$obs)
	dataT <- qsd$qldata[(xdim+1L):(xdim+nstat)]
	dataL <- qsd$qldata[grep("^L[^b]",names(qsd$qldata))]					# Cholesky decomposed terms
	nc <- ncol(dataL)
		
	sig2 <-
	 if(!is.null(theta)) {
		 if(!is.null(cvm) && class(cvm) == "cv") {		  
		 	tryCatch({							
				Y <- estim(qsd$covT,theta,Xs,dataT,krig.type="var")
				# cross-validation variance/RMSE of statistics
				cverrorTx(theta,Xs,dataT,cvm,Y,"cve")		
			 }, error = function(e) { e })
		 
		 } else if((krig.type == "var")) {
		       try(varKM(qsd$covT,theta,Xs,dataT),silent=TRUE)
		 } 
	} else NULL 	
	 
	if(.isError(sig2)) {
		message(.makeMessage("Failed to get prediction variances. "),
		  if(inherits(sig2,"error")) conditionMessage(sig2))		
		sig2 <- NULL
	}	
	if(!is.null(W)) {
		if(!is.matrix(W))
		 stop("`W` has to be a matrix.")
		stopifnot(nrow(W)==ncol(W) && nrow(W)==xdim)
	}
	
	tryCatch({
		## get list of covariance matrices
		## one for each prediction point
		if(var.type %in% c("logMean","wlogMean")) {
			mlogV <- try(varLOGdecomp(dataL[1:length(qsd$covL)]),silent=TRUE)
			if(.isError(mlogV)) {
				msg <- paste0("Matrix logarithm failed.")
				message(msg)
				return(.qleError(message=msg,call=match.call(), error=mlogV ) )
			}			
			if(var.type=="logMean" || is.null(W) || is.null(theta))
			  return (varCHOLmerge(rbind(colMeans(mlogV)),sig2, var.type, doInvert))
		    		  	
			d <- try(exp(-.distX(Xs,rbind(theta),W)*0.5),silent=TRUE)
			if(inherits(d,"try-error") || !is.numeric(d) || any(is.na(d))) {
			   msg <- paste0("Weighted distances calculation error: NAs values possibly produced.")		
			   message(msg)
			   return(.qleError(message=msg,call=match.call(),error=d))
			}
			varCHOLmerge(rbind(colSums(mlogV*matrix(rep(d,nc),nrow(dataL),nc))/sum(d)),
					sig2,var.type,doInvert,Tnames)
			
		} else if(var.type %in% c("cholMean","wcholMean")) {
			if(var.type=="cholMean" || is.null(W) || is.null(theta))
			  return (varCHOLmerge(rbind(colMeans(dataL)),sig2,var.type,doInvert,Tnames))	
			
			# weighting matrix has to be inverted!
			# but inverse of QI as variance of theta inverted is QI, thus W = I.
			d <- try(exp(-.distX(Xs,rbind(theta),W)*0.5),silent=TRUE)
			if(inherits(d,"try-error") || !is.numeric(d) || any(is.na(d))) {
				msg <- paste0("Weighted distances calculation error: NAs values possibly produced.")
				message(msg)
				return(.qleError(message=msg,call=match.call(),error=d))
			} 				
			varCHOLmerge(rbind(colSums(dataL*matrix(rep(d,nc),nrow(dataL),nc))/sum(d)),
					sig2,var.type,doInvert,Tnames)
		} else if(var.type %in% c("kriging","logKrig")) {		
			## the following krigin of variance matrix is only used for
			## functions in R code, any other code calling C does its own
			## kriging more efficiently
			if(is.null(qsd$covL) || is.null(theta))
			  stop("Kriging the variance matrix requires arguments `qsd$covL` and `theta` to be given.")	 
			
			# nCovL: number of Cholesky decomposed terms (excluding bootstrap variances)		
			L <- estim(qsd$covL,theta,Xs,dataL[1:length(qsd$covL)],krig.type="var")			
			Lm <- do.call(rbind,sapply(L,"[","mean"))
			Lsig <- try(sqrt(do.call(rbind,sapply(L,"[","sigma2"))),silent=TRUE)
			if(inherits(Lsig, "try-error") || anyNA(Lsig))
			  stop("Could not extract Kriging variances of variance matrix interpolation models.")		 	
			
		    VTX <- varCHOLmerge(Lm,NULL,var.type,doInvert,Tnames)		
			lapply(seq_len(NROW(Lm)),
					function(i) {
						Sig2 <- if(var.type == "kriging")
						  try(.chol2var(as.numeric(Lsig[i,])),silent=TRUE)
					  	else try(expm::expm(.mergeMatrix(as.numeric(Lsig[i,]))),silent=TRUE)
						if(inherits(Sig2,"try-error")){
						    msg <- paste0("'chol2var' failed to compute the kriging variance of (kriging) variance matrix models.")
							message(msg)
							return(.qleError(message=msg,error=Sig2))
					  	}
					    V <- VTX[[i]]
						res <- 
						 structure(
						  list("VTX"=V$VTX, #V$VTX+diag(diag(Sig2),nstat), 	# with kriging variance of statistics mean 'sig2'				   				  
							   "sig2"=diag(Sig2))					# a matrix: kriging variances of variance matrix models
			   			  )					  
						  if(doInvert) {
							  res$inv <- try(do.call("gsiInv",list(res$VTX)),silent=TRUE)
							  if (inherits(res$inv,"try-error") || !is.numeric(res$inv) || any(is.na(res$inv)))
							   message("Variance matrix inversion failed for kriging the variance matrix.")			
						  }
						  res
				}
			)		  
		} else {			
			stop("Variance matrix of statistics cannot be computed.")
		}	
	}, error = function(e) {
		msg <- paste0("Covariance matrix estimation failed: ", conditionMessage(e),".\n") 
		message(msg)
	    return(.qleError(message=msg,call=match.call(),error=e))
	   }
	)
}

# helpers (intern)
.chol2var <- function(Xs) {
	n <- (-1 + sqrt(1 + 8*length(Xs)))/2;
	m <- matrix(0,n,n)
	m[col(m)>=row(m)] <- Xs
	return( crossprod(m) )
}

.chol2Upper <- function(Xs) {
	n <- (-1 + sqrt(1 + 8*length(Xs)))/2;
	m <- matrix(0,n,n)
	m[col(m)>=row(m)] <- Xs
	return( m )
}

.mergeMatrix <- function(Xs) {
	n <- (-1 + sqrt(1 + 8*length(Xs)))/2;
	m <- matrix(0,n,n)
	m[col(m)>=row(m)] <- Xs
	m[lower.tri(m)] = t(m)[lower.tri(m)]
	return( m )
}

# sig2 is matrix: rows are kriging variance => put as diagonal matrix
varCHOLmerge <- function(Xs, sig2=NULL,var.type="wcholMean",doInvert=FALSE,Tnames = NULL){
	UseMethod("varCHOLmerge",Xs)
} 

varCHOLmerge.matrix <- function(Xs,sig2=NULL,var.type="wcholMean", doInvert=FALSE,Tnames = NULL) {
	if(!is.null(sig2) && is.matrix(sig2) )
		structure(lapply(seq_len(NROW(sig2)),
			function(i) varCHOLmerge(Xs[1L,],sig2[i,],var.type,doInvert,Tnames) ),"var.type"=var.type)
	else structure(list(varCHOLmerge(Xs[1L,],NULL,var.type,doInvert,Tnames)),"var.type"=var.type)
}

# intern
#' @importFrom expm expm
varCHOLmerge.numeric <- function(Xs, sig2=NULL, var.type="wcholMean", doInvert=FALSE, Tnames = NULL) {
   err <- 
	  function(e) {
			message(paste0(.makeMessage("try to invert again...\n")))
			tmp <- try(do.call(gsiInv,list(varMat)),silent=TRUE)
			if (!is.numeric(tmp) || !is.matrix(tmp) || any(is.na(tmp))) {
				msg <- .makeMessage(paste0("Matrix inversion error: "),conditionMessage(e))
				message(msg)
			    return (.qleError(message=msg, call=sys.call(), error=e))
		 	}
			return (tmp)
	  }	
	VTX <- try({
	 if(var.type %in% c("cholMean","wcholMean","kriging")) 
		.chol2var(Xs)
	 else expm::expm(.mergeMatrix(Xs))		
	}, silent=TRUE)		
	if(inherits(VTX,"try-error"))
	  stop(paste0("Failed to merge covariance matrix by: ",var.type))
  	if(!is.null(Tnames))
     dimnames(VTX) <- list(Tnames,Tnames)
    res <- list("VTX"=VTX)

	if(!is.null(sig2)) {
		n <- length(sig2)
		stopifnot(nrow(VTX)==n)
		res$sig2 <- sig2
		if(doInvert) {
			res$inv <- try(do.call("gsiInv",list(res$VTX)),silent=TRUE)
			if (inherits(res$inv,"try-error") || !is.numeric(res$inv) || any(is.na(res$inv)))
			  return(.qleError(message="Variance matrix inversion failed: ",call=sys.call(),error=res))			
		}
	} else {
		if(doInvert) {
			if(var.type %in% c("logMean","wlogMean")) {
				minv <- try( expm::expm(-.mergeMatrix(Xs)),silent=TRUE)
				if(.isError(minv)){
					msg <- paste0("Merge matrix (logarithm) failed: ")
					message(msg)
					return(.qleError(message=msg,call=sys.call(),error=minv,Xs))
				}
				res$inv <- minv
		 	} else {
				minv <-
					tryCatch({
					  varMat <- .chol2Upper(Xs)
					  do.call("chol2inv",list(varMat))
					}, error = err)
				if(.isError(minv))
				  return(.qleError(message="Not a matrix: ",call=sys.call(),error=minv))
				res$inv <- minv
			}
		}
	}	
	return(res)
}

#' @name quasiDeviance
#'
#' @title Quasi-deviance computation
#'
#' @description
#'  The function computes the quasi-deviance (QD) at parameters (called points) of the parameter
#'  search space including the computation of the quasi-score vector and optionally its variance matrix.   
#'
#' @param points		list or matrix of points where to compute the QD (a numeric vector is considered to be a (multidimensional) point)
#' @param qsd		    object of class \code{\link{QLmodel}} 
#' @param Sigma		    variance matrix estimate of statistics (see details)
#' @param ...		    further arguments passed to \code{\link{covarTx}}
#' @param cvm			optional, ether list of covariance models of the statistics for cross-validation based estimation of prediction variances of the statistics or
#' 					    of class \code{cv} or list of cross-validation models of class \code{cvfull} of the QL model \code{qsd} for computation of the quasi-deviance
#' 						and error estimation of the quasi-score approximation w.r.t the kriging prediction models (see \code{\link{prefitCV}})
#' @param obs	 	    numeric vector of observed statistics, this overwrites `\code{qsd$obs}`, if supplied
#' @param w			    numeric value, \code{=0.5} (default) as scalar weight, \code{0<=w<=1}, for evaluation of candidate points
#' @param criterion		name of criterion function, either "\code{qle}" (default) or "\code{mahal}" which overwrites stored type name in `\code{qsd}`
#' @param check			logical, \code{TRUE} (default), whether to check input arguments
#' @param value.only  	if \code{TRUE} only the values of the QD are returned
#' @param na.rm 		logical, if \code{TRUE} (default) remove `Na`s from the result
#' @param cl			cluster object, \code{NULL} (default), of class \code{MPIcluster}, \code{SOCKcluster}, \code{cluster}
#' @param verbose   	logical, \code{TRUE} for intermediate output
#'
#' @return Numeric vector of QD values, if values only, or a list with elements:
#' \item{value}{ either quasi-deviance value or Mahalanobis distance value }
#' \item{par}{ parameter estimate}
#' \item{I}{ quasi-information matrix}
#' \item{score}{ quasi-score vector}
#' \item{jac}{ Jacobian of sample average statistics}
#' \item{varS}{ estimated variance of quasi-score, if applicable}
#' \item{Iobs}{ observed quasi-information}
#' 
#'  The matix `\code{Iobs}` is called the \eqn{\emph{observed quasi-information}} (see [2, Sec. 4.3]),
#'  which, in our setting, can be calculated at least numerically as the Jacobian of the quasi-score vector.
#'  Further, `\code{varS}` denotes the approximate variance-covariance matrix of the quasi-score given the observed
#'  statistics and serves as a measure of estimation precision (see [1] and the vignette, Sec. 3.2).
#'   
#' @details The function calculates the QD (see [1]). It is the primary function criterion to be minimized
#'   for estimating the unknown model parameter by \code{\link{qle}} and involves the computation of the quasi-score
#'   and quasi-information matrix at a particular parameter. From a statistical point of view, the QD can be seen as
#'   a generalization to the \emph{efficient score statistic} (see [3] and the vignette) and is used as a decision
#'   rule in the estimation function \code{\link{qle}} in order to hypothesise about the true model parameter. A modified value of
#'   the QD, using the inverse of the variance of the quasi-score w.r.t. the kriging approximation models of the statistics 
#'   is stored in the result `\code{qval}`.
#'    
#'   Quasi-deviance values which are relatively small (compared to the empirical quantiles of its approximate chi-squared
#'   distribution) suggest a solution to the quasi-score equation and hence could identify the unknown model parameter
#'   in some probabilistic sense. Estimated parameters including different observed statistics can be investigated by 
#'   a MC goodness-of-fit test, see \code{\link{qleTest}}.
#' 
#'   Further, if we use a weighted variance average approximation of statistics (see \code{\link{covarTx}}),
#'   then the QD value is calculated rather locally w.r.t. to an estimated parameter `\code{theta}`. A constant variance
#'   matrix is also applicable to the computation of the QD. However, if supplied, `\code{Sigma}` is used
#'   as is with kriging variances added at the `\code{points}` as diagonal terms.    
#' 
#' 	 \subsection{Use of prediction variances}{ 
#' 	 In order to not only account for the simulation variance but additionally for the approximation error of the
#'   quasi-score vector we include the prediction variances of the involved statistics either based on
#'   a cross-validation or kriging approach. If `\code{cvm}` is not given, then the prediction variances are obtained based
#'   on the kriging procedure applied to the statistics. The variance matrix `\code{varS}` of the
#'   quasi-score vector is part of the return list. Besides the quasi-information matrix the observed quasi-information matrix
#'  (as a numerically derived Jacobian, given by `\code{Iobs}`, of the quasi-score vector) is also returned. A good match between
#'   those two matrices suggests an approximate root if the corresponding QD value is relatively small. This can be further investigated
#'   using the function \code{\link{checkMultRoot}}.
#' 
#'   Alternatively, also CV-based prediction variances (which involve additional covariance models given by `\code{cvm}` for each left out
#'   sample point) for each single statistic can be used to produce relatively robust parameter estimation results but for the price of
#'   much higher computational costs. In practice this might overcome the general tendency inherent to the plug-in kriging predictor underestimating
#'   the prediction variances of the sample means of the statistics. In particular, the CV approach is recommended in case one favours kriging type for
#'   the approximation of the variance matrix. 
#'   }
#' 
#'  \subsection{Mahalanobis Distance}{ 
#' 	 The function computes the Mahalanobis distance of the given statistics \eqn{T(X)\in R^p} with different options
#'  for the approximation type of the variance matrix. The Mahalanobis distance can be used as an alternative criterion function
#'  for estimating the unknown model parameter during the main estimation function \code{\link{qle}}.
#'  
#'  There are several options how to estimate or choose the variance matrix of the statistics \eqn{\Sigma}.
#'  First, in case of a given constant variance matrix estimate `\code{Sigma}`, the Mahalanobis distance reads
#' 	\deqn{ (T(x)-E_{\theta}[T(X)])^t\Sigma^{-1}(T(x)-E_{\theta}[T(X)]) }
#'  and `\code{Sigma}` is used as given.
#' 
#'  As a second option, the variance matrix \eqn{\Sigma} can be estimated by the average approximation 
#'  \deqn{\bar{V}=\frac{1}{n}\sum_{i=1}^n V_i  }
#'  based on the simulated variance matrices \eqn{V_i=V(\theta_i)} of statistics over all sample points
#'  \eqn{\theta_1,...,\theta_n} (see vignette).
#'  Unless `\code{qsd$var.type}` equals "\code{const}" additional prediction variances are added as diagonal terms in order
#'  to account for the kriging approximation error of the statistics using kriging. A weighted version of these average approximation
#'  types is also available (see \code{\link{covarTx}}).
#'  
#'  As a continuous variance approximation we use a kriging approach (see [1]). Then \deqn{\Sigma(\theta) = Var_{\theta}(T(X))}
#'  denotes the variance matrix which depends on the parameter \eqn{\theta\in R^q} and corresponds to the
#'  formal function argument `\code{points}`. Each time a value of the criterion function is calculated at any parameter
#'  `\code{points}` the variance matrix is estimated by the correpsonding approach with added prediction variances as explained above.
#'  Note that in this case the argument `\code{Sigma}` is ignored.
#'  }
#' 
#' @examples
#' data(normal)
#' quasiDeviance(c(2,1), qsd)
#' # (weighted) least squares
#' quasiDeviance(c(2,1), qsd, Sigma=diag(2), criterion = "mahal")
#' 
#' @author M. Baaske
#' @rdname quasiDeviance
#' @export
quasiDeviance <- function(points, qsd, Sigma = NULL, ..., cvm = NULL, obs = NULL, 
					check = TRUE, w = 0.5, criterion = NULL, value.only = FALSE,
						na.rm = TRUE,  cl = NULL, verbose = FALSE)
{		
	if(check)
	 .checkArguments(qsd,Sigma=Sigma) 
 	stopifnot(is.numeric(w)) 		
	if(!is.list(points))
	 points <- .ROW2LIST(points) 	
	xdim <- attr(qsd$qldata,"xdim")	
	
	# Unless Sigma is given always continuously update variance matrix.
	# If using W, theta or Sigma for average approximation, then
	# at least update added kriging prediction variances at each point.
	tryCatch({ 
	  Sigma <- varMatrix(qsd,Sigma,...,cvm=cvm)
	  if(.isError(Sigma)){
		  msg <- .makeMessage("Could not compute variance matrix.")
		  message(msg)
		  return(.qleError(message=msg,error=Sigma))
	  }  			
	  qlopts <- list("varType"=qsd$var.type, "useCV"=!is.null(cvm) && class(cvm)=="cv")
	  # overwrite (observed) statistics	
	  if(!is.null(obs)) {
		  obs <- unlist(obs)
		  if(anyNA(obs) | any(!is.finite(obs)))
		   warning("`NA` or `Inf` values detected in `obs`.")
		  if(!is.numeric(obs) || length(obs) != length(qlm$covT))
		   stop("`obs` must be a (named) `numeric` vector or list of length equal to the given statistics in `qsd`.")
		  qsd$obs <- obs
	  }
	  # overwrite criterion if given
	  if(!is.null(criterion)) {
		  criterion <- match.arg(criterion,c("qle","mahal"))
		  qsd$criterion <- criterion
	  }	  
	  # full QL model
	  args.qsd <- list("qlm"=qsd,"VTX"=Sigma,"X"=as.matrix(qsd$qldata[seq(xdim)]))
	  
	  # LOO CV models patch
      if(!is.null(cvm) && class(cvm) == "cvfull") {
		 cvm <- doInParallel(cvm,
				  function(qlm, criterion, obs,...) {				
					  Sigma <- varMatrix(qlm,...)
					  if(.isError(Sigma))						  
					    return(.qleError(message=.makeMessage("Could not compute variance matrix for CV models."),error=Sigma))					  
					  qlm$obs <- obs
					  qlm$criterion <- criterion				
					  list("qlm"=qlm,"VTX"=Sigma,"X"=as.matrix(qlm$qldata[seq(attr(qlm$qldata,"xdim"))]))				
				  }, criterion=criterion,obs=obs,Sigma=Sigma,...) 
		 hasErr <- which(sapply(cvm,function(x) .isError(x)))
		 if(length(hasErr)>0) {
			 msg <- .makeMessage("Failed to compute variance matrix for cross-validation models.")
			 message(msg)
			 return(.qleError(message=msg,error=cvm))
		 }		  
	  }	  
	 
	  ret <-
		  if(length(points) >= 100 &&
			(length(cl) > 1L || getOption("mc.cores",1L) > 1L))
	      {
				m <- if(!is.null(cl)) length(cl) else getOption("mc.cores",1L)		
				M <- .splitList(points, m)
				names(M) <- NULL
			    fargs <- 
 				 if(qsd$criterion == "qle") {
					c(list(X=M,
						FUN=function(points, args.qsd, qlopts, cvm, value.only, w) {
							.Call(C_quasiDeviance,points,args.qsd,qlopts,cvm,value.only,w)	 
						}, cl=cl),
					list(args.qsd, qlopts, cvm, value.only, w))
				 } else {
					 c(list(X=M,
						 FUN=function(points, args.qsd, qlopts, cvm, value.only, w) {
							 .Call(C_mahalanobis,points,args.qsd,qlopts,cvm,value.only,w)	 
						 }, cl=cl),
					 list(args.qsd, qlopts, cvm, value.only, w))
				}				
				unlist( do.call(doInParallel,fargs), recursive = FALSE)
 							
		} else {
			switch(qsd$criterion,
			 "qle" = .Call(C_quasiDeviance,points,args.qsd,qlopts,cvm,value.only,w),
			 "mahal" = .Call(C_mahalanobis,points,args.qsd,qlopts,cvm,value.only,w) )			
		}		
		# check for NAs
		if(na.rm){
			has.na <- as.numeric(which(is.na(sapply(ret,"[",1))))	
			if(length(has.na) == length(ret)){
				stop("All quasi-deviance computations have `NA` values.")
			}
			if(length(has.na > 0L)){		
				message("Removing `NA` values from results of quasi-deviance calculation.")				
				return( structure(ret[-has.na], "hasNa"=has.na))				
			} else if(length(has.na) == length(points)) {				
				stop(.makeMessage("All return values are `NA` while computing quasi-deviance."))
			}
		}
		return(ret)
	}, error = function(e) {
		message(.makeMessage("Calculation of quasi-deviance failed: ",conditionMessage(e)))
	 	stop(e)  # re-throw error
	})
}

#' @name multiDimLHS
#'
#' @title Multidimensional Latin Hypercube Sampling (LHS) generation
#'
#' @description The function generates or augments a multidimensional design in a hyperbox.
#'
#' @param N		    number of points to randomly select or augment an existing sample set
#' @param lb	    lower bounds defining the (hyper)box of the parameter search space
#' @param ub 		upper bounds defining the (hyper)box of the parameter search space
#' @param method    type of sampling, `\code{randomLHS}`, `\code{maximinLHS}` or `\code{augmentLHS}` 				    
#' @param X 		optional, matrix of existing sample points, \code{NULL} (default), for augmentation only
#' @param type 		either "\code{list}" or "\code{matrix}" as return type
#'
#' @return Either return a list or matrix of sampled vectors or newly generated points if an existing sample set
#'  was to be augmented.
#'
#' @rdname multiDimLHS
#' @importFrom lhs randomLHS maximinLHS augmentLHS
#'
#' @examples
#' data(normal)
#' # generate a design
#' X <- multiDimLHS(N=5,qsd$lower,qsd$upper,type="matrix")
#' 
#' # augment design X 
#' rbind(X,multiDimLHS(N=1,qsd$lower,qsd$upper,X=X,
#' 				method="augmentLHS",type="matrix"))
#' 
#' 
#' @author M. Baaske
#' @importFrom lhs randomLHS maximinLHS augmentLHS
#' @export
#' @seealso \code{\link[lhs]{randomLHS}}, \code{\link[lhs]{maximinLHS}}, \code{\link[lhs]{augmentLHS}}
multiDimLHS <- function(N, lb, ub, method = c("randomLHS","maximinLHS","augmentLHS"),
						X = NULL, type = c("list","matrix"))
{	
	if(!is.null(X) && !is.matrix(X))
	  stop("`X` must be a matrix.")
	dimX <- length(ub)
	if(dimX != length(lb))
	  stop("`lb` and `ub` bounds vector do not match.")
	
    method <- get(match.arg(method))	
	# back to [0,1]
	lhs.grid <- 
		if(!is.null(X)) {				  	
			for(i in 1L:dimX ) {
			 stopifnot(ub[i]!=lb[i])
			 X[,i] <- (X[,i]-lb[i])/(ub[i]-lb[i])		 
		 	}
			# augment X, only return newly generated points
			if( N < 1L)
			  stop("Number of points to augment must be positive (N > 0).")
			rbind(do.call(method,list("lhs"=X,"m"=N))[(nrow(X)+1L):(nrow(X)+N),])		
		} else do.call(method,list("n"=N,"k"=dimX))
	for(i in 1L:dimX )
	 lhs.grid[,i] <- lhs.grid[,i]*(ub[i]-lb[i])+lb[i]	
 	type <- match.arg(type)
 	switch(type,
		"list" = {
			colnames(lhs.grid) <- names(ub)
			return (lapply(seq_len(nrow(lhs.grid)), function(i) as.list(lhs.grid[i,])))
		},
		"matrix" = {
			colnames(lhs.grid) <- names(ub)
			lhs.grid
		}
     )
}

#' @name Subset of statistics
#' 
#' @title Optimal subset selection of statistics
#' 
#' @description The function finds a subset of at most \eqn{kmax <= p} statistics, where \code{p} is the number of available statistics
#' in the list `\code{qsd$covT}` (and at least of size equal to the length \code{q} of the parameter `\code{theta}`) and thus minimizes the expected
#' estimation error of the parameter when this subset is used for estimation. Based on the eigenvalue decomposition of the
#' variance-covariance matrix of the statistics this subset is chosen among all subsets of size at most equal to `\code{kmax}` or for
#' which all proportional contributions to each parameter component are greater than or equal to `\code{cumprop}` whatever happens first.
#' 
#' Since both matrices depend on `\code{theta}` so does the chosen subset of statistics. However, using a list of parameters as `\code{theta}`
#' returns a list of corresponding subsets. One can then easily choose the most frequent subset among all computed ones given either
#' a sample of parameters distributed over the whole parameter space or an appropriate smaller region, where, e.g., the
#' starting point is chosen from or the true model parameter is expected to lie in. 
#' 
#' @param theta 	list or matrix of points where to compute the criterion function
#' 				 	and to choose `\code{kmax}` statistics given the QL model `\code{qsd}`
#' @param qsd		object of class \code{\link{QLmodel}} 
#' @param kmax   	number of statistics to be selected with q <= \code{kmax} <= p
#' @param cumprop	numeric vector either of length one (then replicated) or equal to the length of `\code{theta}` which sets the
#' 				    proportions (0 < \code{cumprop} <= 1) of minimum overall contributions to each parameter component given the statistics 					
#' @param ...		further arguments passed to \code{\link{quasiDeviance}}
#' @param cl		cluster object, \code{NULL} (default), of class \code{MPIcluster}, \code{SOCKcluster}, \code{cluster}
#' @param verbose  	logical, \code{FALSE} (default) for intermediate output
#' 
#' @return A list which consists of 
#' 	\item{id}{ indices of corresponding statistics}
#' 	\item{Tmax}{ named vector of \code{id} corresponding maximum values of statistics}
#'  \item{cumprop}{ cumulated proportions of contributions of selected statistics to each of the parameter components} 
#'  \item{sorted}{ list of statistics (for each parameter) sorted in decreasing order of proportional contributions to the quasi-information}  
#' 
#' @rdname optStat
#' 
#' @examples
#'  data(normal)
#'  # must select all statistics and thus using the
#'  # full information since we only have to statistics available 
#'  optStat(c("mu"=2,"sigma"=1),qsd,kmax=2)[[1]]
#' 
#' @author M. Baaske
#' @export 
optStat <- function(theta, qsd, kmax = p, cumprop = 1, ..., cl = NULL, verbose=FALSE) 
{	
	p <- length(qsd$covT)
	q <- attr(qsd$qldata,"xdim")
	stopifnot(length(cumprop)>0L)
	if(kmax > p || q > kmax)	
	 stop("`kmax` must be at most equal to the number of available statistics and at least equal to the number of model parameter.")	
 	if(length(cumprop) > 1L){
		if(length(cumprop) != q)
		  stop("`cumprop` must be of length equal to number of parmater components or a scalar value.")
	 	else { stopifnot(all(cumprop<=1) && all(cumprop>0)) }
	} else cumprop <- rep(cumprop,q)
	
	# quasi-deviance
	QD <- quasiDeviance(theta,qsd,...,value.only=FALSE,cl=cl,verbose=verbose) 
	
	# evaluate statistics at theta
	ret <- doInParallel(QD,
			FUN=function(x,q,p,kmax,prop) {				
				### TODO: kriging variances for average approximation?	Unterscheidung kriging wcholMean?			
				V <- attr(x,"Sigma")	  
				###
				nms <- colnames(V)
				if(is.null(nms)){
					nms <- paste0("V",1:p)
					attr(V,"dimnames") <- list(nms,nms)
				}
				S <- try(eigen(V),silent=TRUE)
				if(inherits(S,"try-error"))
				  return(S) 
				L <- (x$jac %*% S$vectors %*% diag(1/sqrt(S$values)))^2
				# relative to total contribution
				L <- L/rowSums(L)
				colnames(L) <- nms
				# sort each row as contribution of each statistic to each parameter
				M <- lapply(1:nrow(L),function(i) sort(L[i,], decreasing=TRUE))
				# find either kmax best statistics or until cumulative proportions
				# for each parameter component are greater than prop				
			    T <- sapply(M,"[",1)
				dup <- duplicated(names(T))
				# use only non duplicated names
				if(any(dup)){
				 T[names(T[dup])[1]] <- max(T[dup]) 
				 T <- T[!dup] 
				}
				if( kmax > q && min(T) < min(prop) ) {
					stp <- FALSE
					for(i in 2:p){
						B <- sapply(M,"[",i)
						ix <- order(B,decreasing=TRUE)					
						for(k in ix){
						  if(names(B[k]) %in% names(T))
							next
						  else {
							 T <- c(T,B[k])
							 if(length(T) == kmax || 
								all(sapply(M, function(x) sum(x[names(T)])) >= prop) ){
							   stp <- TRUE
							   break;					     
						   	 }
						  }
					    }
						if(stp) break
					}		
				}
			    names(M) <- names(theta)
				id <- na.omit(pmatch(names(T),colnames(V)))
				# rank matrix
				rankM <- matrix(0,nrow=q,ncol=p)
				dimnames(rankM) <- list(names(theta),nms)
				for(i in 1:nrow(rankM))
				 rankM[i,] <- pmatch(colnames(rankM),names(M[[i]]))
				
				structure(
				  list("id"=id,
					"Tmax"=apply(do.call(rbind,M)[,names(T)],2,max),
					"cumprop"=apply(L, 1, function(x) sum(x[id])),
					"sorted"=M, "rankMat"=rankM))		
			}, q=q, p=p, kmax=kmax, prop=cumprop,
			cl=cl)	
			
	return( ret )
}
