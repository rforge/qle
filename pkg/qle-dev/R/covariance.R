# Copyright (C) 2017 Markus Baaske. All Rights Reserved.
# This code is published under the GPL (>=3).
#
# File: 	covariance.R
# Date:  	12.04.2017
# Author: 	Markus Baaske
#
# Define covariance structures: SIRF-k and Matern

#' @name
#'  setCovModel
#' 
#' @title 	
#' 	Set a covariance model
#' 
#' @description 
#' 	Set a covariance model for kriging the statistics or for the variance matrix of statistics.
#' 
#' @param model			name of covariance model
#' @param param			numeric vector, \code{NULL} (default), starting values of covariance parameters 
#' @param npoints		number of sample points already evaluated for covariance parameter estimation
#' @param var.sim		numeric vector, \code{NULL} (default), local simulation variances used as local nugget variances
#' @param as.nugget	    logical, \code{FALSE} (default), whether to treat `\code{var.sim}` as
#' 						a fixed local nugget variance
#' @param nugget		starting value for (global) nugget value estimation
#' @param trend			trend order number: linear=1 or quadratic=2 for polynomial trend terms 					    
#' @param fixed.param	vector of names, corresponding to `\code{param}` of covariance parameters, hold fixed
#' @param lower			lower bounds for REML estimation 
#' @param upper			upper bounds for REML estimation
#' @param ...			additional arguments which can be stored
#' 
#' @return Object of class \code{covModel}
#' 
#' @details The function sets the a covariance model for kriging the sample mean values of statistic. The covariance model
#'  (including a polynomial trend function) defines the spatial dependence between different locations (points) of the
#'  parameter space. Currently, the function provides the generalized covariance models
#'  (`\code{sirfk}`, see \code{\link{fitSIRFk}}) of order \eqn{k=1,2} 
#'  and the Mat\eqn{\textrm{\'{e}}}rn covariance model with scale (or sill) parameter `\code{scale}`, smoothness parameter
#'  `\code{alpha}`, respectively, `\code{nu}`, and the range parameter `\code{rho}` defined only for the latter. Also, 
#'  a power exponential covariance model, "\code{powexp}", is supported. 
#'  
#'  \subsection{Use of simulation variance}{
#'  If a vector of simulation variances is set by `\code{var.sim}` for each sampling location these are used as local
#'  nugget variance approximations of the sampling variability due to the repeated simulations. The length should match
#'  the number of points `\code{npoints}` otherwise the given components are replicated. Setting `\code{as.nugget}` treats
#'  these values as fixed for all subsequent covariance model estimations. In addition, a global scalar valued nugget, which
#'  captures the variance of the underlying random function, can be set by `\code{nugget}` as a starting value for the REML
#'  estimation. Clearly, both types of nugget variances directly influence the REML estimated covariance parameters.
#'  } 
#' 
#'  \subsection{Default parameters}{ 
#'  The default starting parameters are set to \code{("scale"=0.001,"alpha"=1.5)} for the `\code{sirfk}` model. The
#'  Mat\eqn{\textrm{\'{e}}}rn model uses the following parameters \code{("scale"=1.0,"nu"=2.5,"rho"=3.0)}. The default
#'  parameters for the power exponential covariance model are \code{"scale"=1.0}, (isotropic) range parameter
#'  \code{"phi"=1.0} and power \code{"kappa"=1.5} with \eqn{0<\kappa\leq 2}. 
#'  The corresponding lower and upper bounds are chosen such that the underlying (intrinsic) random function remains
#'  twice continuously differentiable. Further, setting the names of the covariance parameters in `\code{fixed.param}`,
#'  excludes these parameters from subsequent REML estimations such that they remain unchanged. The above settings could
#'  be applicable for a wide range of statistics but, however, generally depend on the kind of statistics to be interpolated
#'  and thus have to be chosen carefully. 
#'  Note that a valid (generalized) covariance model for kriging requires at least \eqn{q+2} design points
#'  for trend order \eqn{k=1} and \eqn{1+(q+1)(q+2)/2} for \eqn{k=2} where \eqn{q} is the problem dimension.
#'  }   
#'  
#' @examples 
#'  # set the standards sirf-2 covariance model
#'  setCovModel("sirfk",npoints=12)
#' 
#' @author M. Baaske
#' @rdname setCovModel
#' @export
setCovModel <- function(model = "sirfk", param = NULL, npoints = 0, var.sim = NULL,
						 as.nugget = FALSE, nugget = 1e-4, trend = 2, fixed.param = NULL, 
						   lower = NULL, upper = NULL,...)
{	
	stopifnot(npoints>0)
	trend.nr <- pmatch(trend,c(1,2))
	if(anyNA(trend.nr)) {
	  msg <- paste0("Invalid drift number. Use either linear=1 or quadratic=2 drift.")
	  message(msg)
	  stop(msg)
  	}
 	
    cov.list <- c("sirfk","matern","powexp")
 	cov.nr <- pmatch(model, cov.list)
	if (is.na(cov.nr)) {
		msg <- paste("Unknown covariance model!",
				"Possible values are", paste(cov.list,collapse = ","))
		message(msg)
		stop(msg)
	}		
	
	fix.nugget <- NULL	
	if(!is.null(var.sim)) {		
	  if(!is.numeric(var.sim) || is.matrix(var.sim))
		 stop("`var.sim` has to be a numeric vector.")
	  fix.nugget <-
		# treated as local/global nugget
		if(length(var.sim) != npoints)
			as.numeric(unlist(lapply(list(var.sim), rep, length.out = npoints)))			
		else var.sim	 				
	}
	
	switch(model,
		"sirfk" = {			
			i <- min(which(trend < c(1,2,3), arr.ind=TRUE))
			param <-
			 if(is.null(param))
			  c("scale"=0.001,"alpha"=1.5,"nugget"=nugget)
			 else {			   
			   if(anyNA(match(names(param),c("scale","alpha"))))
				 stop("Invalid parameter vector for model `sirfk`.")
			   p <- c("scale"=0.001,"alpha"=1.5)  
			   p[names(param)] <- param 
			   c(p,"nugget"=nugget)
			 }
		 	
	        if(is.null(lower)) {
			  lower <- rep(-Inf,length(param))	
			}
			lower <- if(anyNA(lower) | any(!is.finite(lower))) c(1e-6,1.0,0) else lower
			upper <- if(is.null(upper)) {			  
				c(5,i-0.001,3)						
			} else c(upper[1],min(upper[2],i-1e-6),upper[3])			 				
		},		
		"matern" = {
			param <-
			 if(is.null(param))
				c("scale"=1.0,"nu"=2.5,"rho"=3.0,"nugget"=nugget)
			 else {
				 if(anyNA(match(names(param),c("scale","nu","rho"))))
					 stop("Invalid parameter vector for model `matern`.")
				 p <- c("scale"=1.0,"nu"=2.5,"rho"=3.0)  
				 p[names(param)] <- param
			     c(p,"nugget"=nugget)
			 }
			if(is.null(lower) || is.null(upper)) {
				lower <- c(1e-6,2+1e-10,0.1,0)
				upper <- c(10,3,10,10)
			}			
		},
		"powexp" = {
			param <-
			 if(is.null(param))
				param <- c("scale"=1.0,"phi"=1.0,"kappa"=1.5,"nugget"=nugget)
			 else {
				 if(anyNA(match(names(param),c("scale","phi","kappa"))))
					 stop("Invalid parameter vector for model `powexp`.")
				 p <- c("scale"=1.0,"phi"=1.0,"kappa"=1.5)  
				 p[names(param)] <- param				 
				 c(p,"nugget"=nugget)
			 }
	 
			if(is.null(lower) || is.null(upper)) {
				lower <- c(1e-6,1e-4,1e-4,0)
				upper <- c(3,3,1.9999,10)
			}			
		}
	)	
	start <- param
	if(any(start<lower) || any(start>upper)) {
		msg <- paste0("At least one parameter does not match constraints. Check `lower` and `upper`.")
		message(msg)
		stop(msg)
	}
	 
		
	fixed <-  
	 if(is.null(fixed.param)){
		 seq(length(start))
	 } else {
		id <- which(is.na(match(names(start),fixed.param)))
		if(length(id) == 0) {
		  cat("All parameters are excluded estimation.")
		} else {
		  start <- start[id]	    # new start
		  lower <- lower[id]		# new bounds
		  upper <- upper[id]		  
		} 
		id
	 }  	
		
	structure(list("model"= cov.nr,
				   "param"= param,
				   "start"= start,				  				   				   
				   "trend"= trend,
				   "fix.nugget"= fix.nugget,	
				   "as.nugget"= as.nugget,
				   "fixed"= fixed,
				   "lower"= lower,
				   "upper"= upper,...),
			class="covModel"
	)
	
}

#' @name 		reml
#' 
#' @title 		Restricted maximum likelihood (REML)
#' 
#' @description Calculate the REML function value for given covariance parameters 
#' 
#' @param models   	 object of class \code{krige} (list of covariance models) or class
#'					 	 \code{covModel} (a single covariance model), see \code{\link{setCovModel}}
#' @param pars 	     covariance parameter vector (including global scalar nugget value)
#' @param data	  	 data frame of simulated statistics, each column corresponds to a 
#' 					 single covariance model in `\code{models}`
#' @param Xs	  	 matrix of sample points
#' @param verbose	 logical, if \code{TRUE}, print intermediate output
#' 
#' @return List of REML function values
#' 
#' @details Given a list of covariance models the function calculates the REML function values at
#'   parameters `\code{pars}`.
#' 
#' @author M. Baaske
#' @rdname reml
#' @export
reml <- function(models, pars, data, Xs, verbose = FALSE) {	
	stopifnot(is.data.frame(data))
	stopifnot(length(models) == ncol(data))
	
	lapply(1:length(models),
		function(i) {
			tryCatch({		
				 F <- .Call(C_Fmat,Xs,models[[i]]$trend)	
				 P <- .Call(C_Pmat,F)	
				 y <- crossprod(P,data[[i]]) # t(P)%*%z
				 fnREML(pars,y,Xs,P,models[[i]],verbose=verbose)
				}
			  ,error = function(e) {
				  msg <- .makeMessage("REML function evaluation failed: ",conditionMessage(e))
			  	  message(msg) 
				  structure( list( message=msg,  call = sys.call() ), class=c("error","condition"), error=e )
			 	}	  
			)			
		})
}

# intern, covModel has to be initialized! 
fnREML <- function(p, y, Xs, P, model, fixed = seq(length(p)), verbose = FALSE)
{	
	 tryCatch({			
		model$param[fixed] <- p		
		# covariance matrix
		Cmat <- .Call(C_covMatrix,Xs,model)		
		if (!is.matrix(Cmat) || anyNA(Cmat) || !is.numeric(Cmat))
		 stop("Covariance matrix has `NA`s.")								
			
	 	W <- crossprod(P, Cmat %*% P)
		rc <- rcond(Cmat)
		msg <- NULL	
		if(rcond(W) < 1e-10 || rc < 1e-10) {		  
		  msg <- paste(c("Covariance matrix reciprocal condition number: ", rc,"\nnear zero at parameter: \n\t ",
						 format(p, digits=6, justify="right"),"\n"),
				  collapse = " ")
		  warning(msg)
	    }
	  	
		Wc <- chol(W)
		w  <- backsolve(Wc,y,transpose=TRUE)	
	
	    structure(
		   as.vector(.5*(t(y)%*%w) + sum(log(diag(Wc)))),
			  "info"=list("rcond"=rc,"msg"=msg, "p"=p) )
	  
	} ,error = function(e) {
			msg <- .makeMessage("Error in function 'fnREML': ",
					conditionMessage(e))
		    message(msg)						
			stop(msg)									
		}
	)		
	
}  

## TODO: Numerical derivatives of log likelihood
##  -> use gradient of loglik
## REML ->alpha
## alpha <- 1.5
## phi <- 1.0
## D <- as.matrix(dist(X))
## diag(D)<- (fix.nugget+nugget)
## H <- log(phi*D)
## exp(2*alpha*H)*2*H
#' @importFrom nloptr nl.grad
fnGradREML <- function(p, y, Xs, P, model, fixed = NULL, verbose = FALSE) {
	list("objective"=fnREML(p,y,Xs,P,model,fixed,verbose),
		 "gradient"=nloptr::nl.grad(p, fnREML, heps = .Machine$double.eps^(1/3),y,Xs,P,model,fixed)) 
}  

## TODO add data as parameter
# arguments '...' manipulate global options for nloptr 
#' @importFrom nloptr nl.opts
doREMLfit <- function(model, Xs, opts, verbose = FALSE )
{
	# return if all parameters are fixed
	if(!is.null(model$fixed) && length(model$fixed)==0L) {
	  return(
		structure(
		  list(model = model,convergence = 1L),
		optres=NULL, class = "reml") )
	}	
	err <- NULL	
	fn <- fnREML
	if(!is.null(opts$algorithm) &&
	   opts$algorithm == "NLOPT_LD_LBFGS") {
	   message("Caution: using `LBFGS` in REML function might not be stable.")
	   fn <- fnGradREML
 	}
	 
 	tryCatch({		
        nms <- names(model$param)				
		Fmat <- .Call(C_Fmat,Xs,model$trend)		
		P <- .Call(C_Pmat,Fmat)	
		
		# transform to mean zero data				
		y <- crossprod(P,as.numeric(model$dataT))
		
		model$dataT <- NULL		
		p0 <- .PROJMED(model$start,model$lower,model$upper)
		
		res <- nloptr::nloptr(p0, fn, lb = model$lower, ub = model$upper, opts = opts,
						y = y, Xs = Xs, P = P, model = model, fixed = model$fixed,
						verbose = verbose)
		msg <- "Normal convergence."
		if(inherits(res,"error") || is.null(res) || anyNA(res$solution)){
			msg <- .makeMessage("Function call to 'nloptr' failed.")				
			message(msg)
			return(.qleError(message=msg,call=match.call(),error=res))
		}
		
		converged <- FALSE
		sol <- res$solution			
		if(res$status >= 0L) {
			converged <- TRUE
			model$param[model$fixed] <- sol						   			   
	    } else {
		   verbose <- TRUE
		   msg <- .makeMessage("Estimation of covariance parameters did not converge.")		   
		   message(msg)
	    }		
	    structure(
		    list(model=model,
				 convergence=converged,
				 message=msg),
		  optres = if(verbose) res else NULL, class = "reml") 
				 
	 }, error = function(e){
			 msg <- .makeMessage("Nloptr error fitting covariance parameters: ",
					  conditionMessage(e))
			 message(msg)
			 structure(
				list(model=model,
					 convergence=FALSE,
					 message=msg,
					 call=sys.call()),
			  error=e)			  		
		}
	) 	 
}


#' @name fitCov
#' 
#' @title Fitting covariance models by REML
#' 
#' @description The function estimates the (hyper)parameters of the covariance models by
#' 	  the \emph{Restricted Maximum Likelihood} (REML) method .
#' 
#' @param models  	 object of class \code{krige} (list of covariance models) or
#' 	 				 class \code{covModel} (a single covariance model), see \code{\link{setCovModel}}
#' @param Xs	 	 matrix of sample points
#' @param data		 data frame of simulated sample means of statistics
#' 					 first column corrspond to first model in list `\code{models}` and so forth
#' @param controls	 list of control arguments, see \code{\link[nloptr]{nloptr}}
#' @param cl		 cluster object, \code{NULL} (default), see \code{\link[parallel]{makeCluster}}
#' @param verbose 	 logical, \code{TRUE} for intermediate output
#' 
#' @return An object of class \code{reml} which consists of a list of named lists
#'  (`\code{model}`, `\code{convergence}`) each storing a fitted covariance model itself
#'  together with the optimization results from \code{\link[nloptr]{nloptr}} as an attribute
#'  named `\code{optres}`. The default method for estimating the covariance parameters is \code{\link[nloptr]{mlsl}}.  
#' 
#' @details The function fits a list of covariance models using the REML method. In order to avoid singularities
#'  of the so-called trend matrices make sure to use at least the minimum required number of sample points in
#'  `\code{Xs}`, which depends on the defined trend order (see \code{\link{setCovModel}} for details).
#' 
#' @examples 
#'   
#' data(normal)
#' 
#' # fit 1st statistic and get REML results
#' fitCov(qsd$covT[1],
#'        Xs=as.matrix(qsd$qldata[1:2]),
#'        data=qsd$qldata["mean.T1"],verbose=TRUE)
#'   
#' @author M. Baaske
#' @rdname fitCov 
#' @export 
fitCov <- function(models, Xs, data, controls = list(),
			     	  cl = NULL, verbose = FALSE) {
		
	if(!is.data.frame(data))
		stop("Expected argument `data` of class `data.frame`.")	
	if(!is.matrix(Xs))
		stop("Expected argument `Xs` to be  a matrix of sample locations.")
	
	opts <- nloptr::nl.opts()	
	if(length(controls)>0L) {		
		opts[names(controls)] <- controls
	} else {
		opts <- list("algorithm" = "NLOPT_GN_MLSL",
				"local_opts" = list("algorithm" = "NLOPT_LN_COBYLA","ftol_rel" = 1.0e-6,
						"xtol_rel" = 1.0e-6,"maxeval" = 1000),
				"maxeval" = 200, "xtol_rel" = 1.0e-6, "ftol_rel" = 1.0e-6, "population"=0)	
	}	
	for(i in 1:length(models))
	 models[[i]]$dataT <- as.numeric(data[[i]])
		
 	mods <- doInParallel(models, doREMLfit, Xs=Xs, opts = opts,
			  cl=cl, verbose=verbose)
	  
	if(inherits(mods,"error")) {
		msg <- paste0("REML estimation failed: ",conditionMessage(mods),"\n")
		message(msg)
		return(.qleError(message=msg,
				call=match.call(),error=mods))
	}	
	errId <- which(sapply(mods,function(x) .isError(x)))
	if(verbose) {	  
	  if(any(errId))
		message(paste(c("Failed fitting covariance models with index: ",as.character(errId)), collapse=" ")) 
	  else {
		id <- which(sapply(mods,function(x) x$convergence))
		if(!all(id)) {
		 message(paste(c("REML failed to converge: ",as.character(id)), collapse=" ")) 
		} else
		 message("Successfully fitted covariance parameters.","\n")		
	  }
	}	
	structure(mods,
		opts = opts,
		error = if(length(errId) > 0L) errId else NULL)
}

#' @name QLmodel
#' 
#' @title Construct quasi-likelihood approximation 
#' 
#' @description Aggregate and construct the data for quasi-likelihood estimation
#' 
#' @param qldata		simulation data (see \code{\link{setQLdata}})
#' @param lb		    numeric vector of lower bounds defining the (hyper)box
#' @param ub 			numeric vector of upper bounds defining the (hyper)box
#' @param obs	    	numeric vector of (real data) statistics
#' @param mods			list of (fitted) covariance models (see \code{\link{fitSIRFk}}) 
#' @param nfit			number of cycles, \code{nfit=1} (default), after which covariance
#' 						parameters are re-estimated otherwise re-used 
#' @param cv.fit 		logical, \code{TRUE} (default), whether to re-fit CV models (re-estimate covariance parameters)	
#' @param var.type  	name of the variance approximation method (see \code{\link{covarTx}})
#' @param useVar    	logical, \code{TRUE} (default), whether to use prediction variances (see details)
#' @param criterion 	global criterion function for sampling and minimization, either "\code{qle}" or "\code{mahal}"				    	
#' @param verbose       logical, \code{FALSE} (default), whether to give further output 
#' 
#' @return Object of class \code{\link{QLmodel}} which stores the data frame of simulation results, bounds on
#'  the parameter space, covariance models for kriging, observations vector as well as options for kriging and fitting. 
#' 
#' @details The function aggregates all information for estimation storing the fitted
#'   covariance models of statistics and the type of variance matrix approximation. For advanced users
#'   this function explicitly offers the data structure to construct individual covariance models for
#'   each statistic as defined by \code{\link{setCovModel}}. The user has the choice whether or not to
#'   make use of of kriging prediction variances by `\code{useVar}` later for approximation
#'   of the variance matrix of statstics. If this is \code{TRUE}, then a kriging procedure calculating prediction 
#'   variances is used. Otherwise we use the so-called \emph{dual} approach which has some computational advantage when
#'   no prediction variances are needed.
#' 
#' 	 See \code{\link{getQLmodel}} for an example. 
#' 
#' @author M. Baaske
#' @rdname QLmodel
#' @export 
QLmodel <- function(qldata, lb, ub, obs, mods, nfit = 1, cv.fit = TRUE,
		    var.type = c("wcholMean","cholMean","wlogMean","logMean","kriging","const"),
				useVar = TRUE, criterion = c("qle","mahal"), verbose = FALSE)
{	
	if(!inherits(qldata,"QLdata"))
	 stop("expected argument `qldata` of class `QLdata`.")
	if(missing(lb) || missing(ub))
	 stop("Arguments `lb` and `ub` are missing.")
 	
 	dx <- attr(qldata,"xdim")
    if(!is.numeric(lb) || !is.numeric(ub) ||
	   length(lb)!=length(ub) || dx != length(ub))
  	  stop("Dimensions of `lb` or `ub` do not match.")
 	
 	obs <- unlist(obs)
	if(anyNA(obs) | any(!is.finite(obs)))
	 stop("`NA`,`NaN` or `Inf`values detected in argument `obs.")
	if(!is.numeric(obs))
	  stop("Argument `obs` must be a (named) numeric vector or list.")
  	stopifnot(!is.null(mods$covT))
	
	if(length(mods$covT) != length(obs))
	  stop("Number of covariance models `covT` and length of observations vector `obs` must equal.")
  	var.type <- match.arg(var.type)
	criterion <- match.arg(criterion)
		
	if(is.null(mods$covL) && var.type == "kriging")
	  stop("Covariance models for variance matrix interpolation must be set for argument \'var.type\'.")
	if(!is.numeric(nfit) || length(nfit)>1L )
	 stop("Argument 'nfit must be numeric of length one.")	
	
 	covT <- .extractCovModels(mods$covT,verbose)
	stopifnot(class(covT)=="krige")
	
	covL <- NULL
	if(!is.null(mods$covL)){		
		covL <- .extractCovModels(mods$covL,verbose)
		stopifnot(class(covL)=="krige")
	}
	# reml optimization options 
	opts <- attr(mods,"opts")
	if(is.null(opts) || length(opts) == 0L){
		opts <- list("algorithm" = "NLOPT_GN_MLSL",
				"local_opts" = list("algorithm" = "NLOPT_LN_COBYLA","ftol_rel" = 1.0e-6,
						"xtol_rel" = 1.0e-6,"maxeval" = 1000),
				"maxeval" = 200, "xtol_rel" = 1.0e-6, "ftol_rel" = 1.0e-6, "population"=0)			  
	}
	# minimum required sample size
	minN <- ifelse(min(sapply(covT,	function(x) x$trend)) < 2, dx+2, (dx+1)*(dx+2)/2+1)
	if(nrow(qldata)<minN) {
	 stop(paste0("Choose the size of the initial sample for parameter dimension ",dx,
			" at least of size: ",minN))
 	}	
	
	structure(
	    list("qldata" = qldata,
			 "lower" = lb,
			 "upper" = ub,
			 "covT" = covT,
			 "covL" = covL,			
			 "obs" = obs,		 
			 "var.type" = var.type,			 
			 "krig.type" = if(useVar) "var" else "dual",
			 "criterion" = criterion,			 
			 "minN" = minN,
			 "nfit" = nfit,
			 "cv.fit"=cv.fit),
	  opts = opts,
	  class="QLmodel"
	)
}

# intern
.extractCovModels <- function(covs, verbose = FALSE) {
	if(is.null(covs))
	  return (NULL)
    # which one procudces errors in fitting?
	errId <- which(sapply(covs,function(x) .isError(x)))
	
	if(length(errId)>0L)
	  message(.makeMessage("A total of ",length(errId)," errors detected in fitted covariance models.")) 
	
    structure(
		lapply(covs,
		 function(x) {
			if(verbose)
			 structure(x$model,"optres"=attr(x,"optres"))
			else x$model
		 }
	    ), error = if(length(errId)>0L) errId else NULL,
	  class = "krige"
 	)	
}

#' @name fitSIRFk
#' 
#' @title Covariance model fitting
#' 
#' @description Fit a generalized covariance model to simulation data
#' 
#' @param qldata		object of class \code{QLdata}, a data frame from \code{\link{setQLdata}}
#' @param set.var 		if \code{TRUE} (default), set simulation variances as local nugget variances 
#' @param var.type      name of variance matrix approximation type (see \code{\link{covarTx}})  
#' @param var.opts	    list of arguments passed to \code{\link{setCovModel}}
#' 						(only for `\code{var.type}` equal to "\code{kriging}")
#' @param intrinsic 	logical, if \code{TRUE}, use a nugget variance estimate for kriging
#' 						approximations of the variance matrix
#' @param ...			arguments passed to \code{\link{setCovModel}}
#' @param cl			cluster object, \code{NULL} (default), see \code{\link[parallel]{makeCluster}}
#' @param controls		list of control parameters passed to \code{\link[nloptr]{nloptr}} minimization function
#' @param verbose		if \code{TRUE}, print intermediate information
#' 
#' @return List of fitted covariance models for kriging the 
#'  sample means of statistics named `\code{covT}` and optionally
#'  the variance matrix of statistics, `\code{covL}`. The object also stores
#'  the reml optimization parameters `\code{controls}`. 
#' 
#' @details The function estimates the parameters of the covariance model `\code{sirfk}` by default using REML for kriging
#'   the statistics and kriging the variance matrix of statistics only if `\code{var.type}` does not equal "\code{const}".
#'   We use a (self-similar) intrinsic random function of order \eqn{k} (see, e.g. [1]) with \eqn{k=1,2} for all statistics
#'   (including a default quadratic drift term, \eqn{k=2}). The user can also define different covariance models for each statistic
#'   separately (see \code{\link{setCovModel}}). If the user prefers to fit all statistics by other covariance models, then these
#'   can be specified by their names in `\code{model}`. Further arguments are passed to \code{\link{setCovModel}}
#'   by `\code{...}`.
#'   		
#'   The argument `\code{var.opts}` only sets the options for the covariance models for kriging the variance matrix.
#'   The optional arguments `\code{var.sim}` and `\code{var.opts$var.sim}` set the local or global
#'   \dfn{nugget} values for each sample point depending on `\code{set.var}`. Both arguments (passed to \code{\link{setCovModel}})
#'   must be data frames of lengths corresponding to the number of covariance models of statistics and, respectively, to the number of
#'   \emph{Cholesky} decompositions in case of kriging the variance matrix. If `\code{set.var}` is \code{TRUE}, then the values in
#'  `\code{var.sim}` are used as fixed `nugget` values and replicated to match the number of sample points if required. Otherwise
#'   they are considered as simulation variances and hence scaled by 1/\code{nsim}, which is meaningful only for kriging the statistics.
#'   The same principle applies in case of kriging the variance matrix. Here the values in `\code{var.opts$var.sim}` (of length one or
#'   equal to the number of corresponding sample points) are used as scale factors for all Cholesky decomposed terms
#'   if `\code{intrinsic}` is \code{TRUE} and otherwise considered as local nugget variances. A global nugget value is estimated during
#'   the REML covariance parameter estimation for both cases.
#' 	 Note that the returned object can also be constructed manually and passed as an input argument to
#'   \code{\link{QLmodel}} in case the user prefers to set up each covariance model separately. Also see function
#'   \code{\link{getQLmodel}} for an example. The default method for estimating the covariance parameters is \code{\link[nloptr]{mlsl}}.
#'  
#' @author M. Baaske
#' @rdname fitSIRFk
#' @export
fitSIRFk <- function(qldata, set.var = TRUE, var.type = "wcholMean",
						var.opts = list("var.sim"=1e-6), intrinsic = FALSE, ...,
						 controls = list(), cl = NULL, verbose = FALSE)
{	
	args <- list(...)
	stopifnot(is.data.frame(qldata))
	stopifnot(is.logical(set.var) && is.logical(intrinsic))
	
	if(length(args)>0L) {		 
		nms <- formals(setCovModel)
		.checkfun(setCovModel, c(nms,args))			
 	}
	
	xdim <- attr(qldata,"xdim") 									# dimension of parameter to estimate!
	nsim <- attr(qldata,"nsim")										# number of simulation replications
	Xs <- data.matrix(qldata[seq(xdim)])	
	dataT <- qldata[grep("^mean.",names(qldata))]					# simulated statistic data
		
	np <- nrow(Xs)	
	nstat <- ncol(dataT)	
	
	# all cov models are equal	
	useVarSim <- !is.null(args$var.sim)	
	# if not used as a fixed nugget: set.var == FALSE	
	set.var <- rep(set.var,length.out=ncol(dataT))
	dfvar <-
	 if(useVarSim) {		
	  rep(as.data.frame(args$var.sim),length.out=ncol(dataT))		
	} else NULL
	
	covT <- 
		lapply(1:ncol(dataT),
		     function(i){
			   args$var.sim <-
				 if(set.var[i]) {					  
					  args$as.nugget <- FALSE
					  qldata[[xdim+nstat+i]]/nsim
				 } else if(useVarSim) {					  
					  args$as.nugget <- TRUE
					  dfvar[[i]]
				  } else NULL
				 fnargs <- c(list("dataT"=dataT[[i]],	      		    # temporarly add the data				  				  			  
								  "npoints"=np,
								  "type"="statcov"), args)
				 do.call(setCovModel,fnargs)					   
			}
		)

	 covL <- NULL	 
	 if(var.type == "kriging"){
		 # check input
		 args <- var.opts
		 if(length(args)>0L) {
			 nms <- formals(setCovModel)			
			.checkfun(setCovModel, c(nms,args))				 
	 	 }		 	   		   		 
  		 useVarSim <- !is.null(args$var.sim)
		 Lvec <- qldata[grep("^L+",names(qldata))]		 
		 # individually set intrinsic noise terms as nugget		
		 intrinsic <- rep(intrinsic,length.out = ncol(Lvec))		 
		 dfvar <- 
		   if(useVarSim) {
			   rep(as.data.frame(args$var.sim),length.out=ncol(Lvec))			   
		   } else NULL		   
		 
		 covL <- lapply(1:ncol(Lvec),
				  function(i)  {
				   args$var.sim <-
					 if(intrinsic[i]) {
						 args$as.nugget <- FALSE
						 ## TODO: improve! estimate nugget variance
						 args$ptol <- if(useVarSim) as.numeric(dfvar[[i]]) else 1
						 args$ptol*Lvec[[i]]						 
					  } else if(useVarSim) {						 
						 args$as.nugget <- TRUE
						 dfvar[[i]]
					 } else NULL
					 	 	    
					 fnargs <- c(list("dataT"=Lvec[[i]],									  		  
									  "npoints"=np,
									  "type"="kriging"), args)							   	  
					 do.call(setCovModel,fnargs)					 
				  }
		 )	 
	 }		 	 
	 # (default) reml optimization options 
	 if(length(controls)>0L) {		
		 opts <- nloptr::nl.opts()
		 opts[names(controls)] <- controls
	 } else {
		 opts <- list("algorithm" = "NLOPT_GN_MLSL",
				  "local_opts" = list("algorithm" = "NLOPT_LN_COBYLA","ftol_rel" = 1.0e-6,
						 "xtol_rel" = 1.0e-6,"maxeval" = 1000),
				 "maxeval" = 200, "xtol_rel" = 1.0e-6, "ftol_rel" = 1.0e-6, "population"=0)		  
	 }
	 	 
	 # REML fit covariance models (statistics and variance matrices)
	 mods <- doInParallel(c(covT,covL), doREMLfit, Xs=Xs, opts = opts,
			 	cl=cl, verbose=verbose)
		
	 if(inherits(mods,"error")) {
		msg <- paste0("REML estimation failed: ",conditionMessage(mods),"\n")		
		message(msg)
		return(.qleError(message=msg,
				 call=match.call(),error=mods))
	 }
	 errId <- which(sapply(mods,function(x) .isError(x)))
	 if(any(errId)) {
		 msg <- paste(c("Failed fitting covariance parameters: ",
				   as.character(errId)), collapse=" ")
   		 message(msg)
		 return(.qleError(message=msg,error=mods[errId]))
	 } else {
		 if(verbose)
		   cat("Successfully fitted covariance parameters.\n")		 
	 }
	ret <- structure(
			  list("covT" = mods[1:nstat],
			  	   "var.type" = var.type),
		     opts = opts,
			 error = if(length(errId)>0L) errId else NULL)	
	 
	if(!is.null(covL))
	  ret$covL <- mods[(nstat+1):length(mods)]
	
	return ( ret )	
}

#' @name getQLmodel
#' 
#' @title Setup the quasi-likelihood estimation model  
#' 
#' @description  Setup the quasi-likelihood model data
#'
#' @param runs   	object of class \code{simQL}, the simulation results (see \code{\link{simQLdata}})
#' @param lb		lower bounds defining the (hyper)box
#' @param ub 		upper bounds defining the (hyper)box
#' @param obs		numeric vector of observed statistics
#' @param X   		matrix of sample locations (model parameters)
#' @param useVar   	logical, \code{TRUE} (default), whether to use prediction variances
#' @param criterion the criterion function to be minimized for parameter estimation (see \code{\link{qle}})
#' @param ...		arguments passed to \code{\link{fitSIRFk}} for fitting covariance models
#'  
#' @return Object of class \code{\link{QLmodel}}
#' 
#' @details The function is a wrapper to \code{\link{setQLdata}} and \code{\link{fitSIRFk}}
#'  in order to setup the data required for estimating the model parameters.
#' 
#' @example /inst/examples/normal.R
#' 
#' @author M. Baaske
#' @rdname getQLmodel
#' @export
getQLmodel <- function(runs, lb, ub, obs, X = NULL, useVar = TRUE, criterion = "qle", ...)
{			 
	# default function:
	# using variances for REML and prediction
	# for all covariance models
	args <- list(...)
	verbose <- isTRUE(args$verbose)
	useChol <-
		if(!is.null(args$var.type) && args$var.type == "const") {
		  FALSE
	    }  else TRUE		
			
	tryCatch({
        if(.isError(runs))
		  stop("Simulations have errors. Please check the input argument `runs`.")
		if(verbose)
		  cat("Collect data for fitting covariance models of statistics.\n")
		qldata <- setQLdata(runs,X,chol=useChol,na.rm=TRUE,verbose=verbose)
		if(.isError(qldata)) {
			return(qldata)
		}
		# fitting statistics
		if(verbose)
		 cat("Fitting covariance models...\n")	 	
	    id <- which(is.na(pmatch(names(args),names(c(formals(fitSIRFk),formals(setCovModel))))))
	 	args.tmp <- if(length(id)>0L) {
			args[-id]
		} else args	    
		
		# fitting			 	
	    mods <- do.call(fitSIRFk,c(list(qldata),args.tmp))	
		if(.isError(mods)) {
			message(.makeMessage("Failed fitting covariance models: ","\n"))			
			return(mods)
		}		  
		if(verbose)
		  cat("Setup QL approximation model...\n")
	  	id <- which(is.na(pmatch(names(args),names(formals(QLmodel)))))
	  	if(length(id) > 0L) {
		  args <- args[-id]		  
	 	}
	    do.call(QLmodel,
			 c(list(qldata,
					lb,ub,
			  	    obs,
				    mods,				    
					useVar=useVar,
				    criterion=criterion),
		     args))
		
	}, error = function(e) {
		msg <- .makeMessage("Failed to setup QL model: ",conditionMessage(e))
		message(msg)
		return(.qleError(message=msg,
				 call=match.call(),error=e))	
	   }
	)
}

#' @name updateCovModels
#' 
#' @title Update covariance models
#' 
#' @description The function updates the current covariance models
#'  stored in `\code{qsd}`.
#' 
#' @param qsd			object of class \code{\link{QLmodel}} (to be updated)
#' @param nextData		object of class \code{QLdata} with new simulation results
#' @param fit 			logical, if \code{TRUE} (default), re-estimate covariance parameters
#' @param cl			currently ignored
#' @param controls	    list of control parameters passed to \code{\link[nloptr]{nloptr}}
#' @param verbose 		logical, whether to print intermediate information
#' 
#' @return Object of class \code{\link{QLmodel}} as a list of updated covariance models
#' 
#' @details The function updates both, the covariance models for kriging the statistics, and, if applicable,
#'  the ones for kriging the variance matrix of statistics. In practice, the user hardly needs to call this
#'  function except for empirical studies of how additional sample points might influence the overall predictive
#'  quality of the quasi-score and/or criterion function approximations.
#' 
#'  If `\code{fit}` is \code{TRUE}, then the function re-estimates the covariance parameters for each statistic separately
#'  each time a total of `\code{qsd$nfit}` new sample points have been added. Thus, we can choose whether to fit the updated
#'  covariance models each time, `\code{qsd$nfit}`=1, or after each 2nd, 3rd, and so on newly added point in order to save some
#'  computational time.
#' 
#' @examples
#' 
#' data(normal)
#' # old design
#' X <- as.matrix(qsd$qldata[c(1,2)])
#' 
#' # augment design
#' Xnew <- multiDimLHS(N=10,qsd$lower,qsd$upper,X=X,
#'            method="augmentLHS",type="matrix")
#' 
#' # new simulations
#' Xsim <- simQLdata(sim=qsd$sim,nsim=100,X=Xnew)
#' 
#' # prepare data
#' Xdata <- setQLdata(Xsim,Xnew,chol=TRUE)
#' 
#' # re-estimate covariance parameters
#' qsd2 <- updateCovModels(qsd,Xdata,fit=TRUE) 
#'  
#' 
#' @author M. Baaske
#' @rdname updateCovModels
#' @export
updateCovModels <- function(qsd, nextData, fit = TRUE,
						cl = NULL, controls = list(), verbose=FALSE)
{		
	stopifnot(class(qsd) == "QLmodel")
	stopifnot(inherits(nextData,"QLdata"))	
	
	nnew <- NROW(nextData)	
	xdim <- attr(qsd$qldata,"xdim")
	nsim <- attr(qsd$qldata,"nsim")
	nsim.new <- attr(nextData,"nsim")
	
	nstat <- length(qsd$covT)	
	stid <- (xdim+1):(xdim+nstat)
	vid <- c((xdim+nstat+1):(xdim+2*nstat))
	
	vars <- qsd$qldata[vid]
	vars.new <- nextData[vid]
	
	# combine old and new data and
	# check columns also because L+ might not be given
	if(ncol(nextData) != ncol(qsd$qldata))
	  stop("The number of columns of argument `nextData` does not match with `qldata`.")
	
	# merge to new data, one (sample point) added
	qsd$qldata <- rbind(qsd$qldata,nextData)	
	Xs <- as.matrix(qsd$qldata[seq(xdim)])
	np <- nrow(Xs)
		
	if(length(controls)>0L) {		
		opts <- nloptr::nl.opts()
		opts[names(controls)] <- controls
	} else {
		opts <- attr(qsd,"opts")	
	}		 
	# update function 
	fitit <- (fit && !(nrow(Xs) %% qsd$nfit))
	.update <- function(covT, data, vars.new=NULL){
		mod <- lapply(1:length(covT),		
				function(i) {		   
					xm <- covT[[i]]
					# set starting point
					xm$start <- xm$param[xm$fixed]									   
					
					if(!is.null(xm$fix.nugget)) {					  
					  xm$fix.nugget <-
						if(xm$as.nugget) {
						  c(xm$fix.nugget,rep(xm$fix.nugget[1],nnew))						
					  	} else if(!is.null(vars.new)) {							
						  c(xm$fix.nugget,vars.new[[i]]/nsim.new)
						} else if(!is.null(xm$ptol)) { # type == "kriging"						   						 
						  c(xm$fix.nugget, xm$ptol*data[[i]][-(1:(np-nnew))] )							
						} else NULL
		   			} # else (not using simulation variance for REML)
					if(fitit) {
						# temporarly store statistics (removed in 'doREMLfit')
						xm$dataT <- data[[i]]										      			
					}
					xm
				})	
		if(fitit) {  				
			res <- doInParallel(mod, doREMLfit, Xs=Xs, opts=opts,
					 cl=cl, verbose=verbose)		
			if(!inherits(res,"error")) {
				.extractCovModels(res,verbose)
			} else {
				msg <- paste0("Error fitting covariance parameters.")
				message(msg)
				structure("message"=msg,"error"=res)				
			}	
		} else structure(mod, class = "krige")
	}
	
	tryCatch({
	  qsd$covT <- .update(qsd$covT,qsd$qldata[stid],nextData[vid])	
	  # update kriging VARIANCE models
 	  # Cholesky terms are the data
	  if(qsd$var.type == "kriging"){
		if(is.null(qsd$covL))
		  stop("A covariance model for kriging the variance matrix must be defined but is `Null`.")
		qsd$covL <- .update(qsd$covL,qsd$qldata[grep("^L+",names(qsd$qldata))])
  	  }
 	}, error = function(e) {
	     msg <- .makeMessage("Failed to update covariance models: ",
				 conditionMessage(e))		
		return(.qleError(message=msg,call=match.call(),error=e))
	})

	return( qsd )
}
