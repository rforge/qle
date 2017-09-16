# Estimation of MC based error measures
# 	MC hypothesis testing with QD as test statistic
# 	Chi-squared goodness-of-fit test

## TODO:
## Restricted (conditioned) MC testing
## (linearly constrained)

# collect test results
.qleTest <- function(obj,...) {	
	degf <- nrow(obj)
	sb <- attr(obj,"sb")
	Sb <- attr(obj,"Sb")	
	pb <- (1+sum(Sb>=sb))/(length(Sb)+1)  # similar to : 1-ecdf(Sb)(sb)
	
	ans <- list()
	ans$param <- cbind(obj[c("par","se","rmse","bias","mean")])
	dimnames(ans$param) <- list(row.names(obj),
			c("Estimate", "Std. Error", "RMSE", "Bias", "Mean"))
	
	tname <- attr(obj,"test")
	nm <- ifelse(tname=="mahal", "LS criterion", "Score-test")
	Stest <- noquote(paste("Monte Carlo Test (",nm,"):\n\n  Df: ",degf,sep=""))	
	
	sb <- noquote(cbind(sb, pb))# , ifelse(degf>0,pchisq(sb, degf, lower.tail = FALSE),"****")))	
	dimnames(sb) <- list("Test E(Q)=0:  ", c("S > 0", "P-value")) #, "Pr(>Chisq)"))
	ans$Stest <- Stest	
	ans$test <- sb
	ans$tname <- tname
	attr(ans$test,"Sb") <- Sb
	class(ans) <- c("qleTest")
	return(ans)
}

#' @name		checkMultRoot
#' 	
#' @title       Inspect estimated parameters
#' 
#' @description		Check out and compare estimated roots of the quasi-score vector
#' 
#' @param est 		object of class "\code{qle}", the estimation results from function \code{\link{qle}}
#' @param par	    list or matrix of estimated parameters as roots of the quasi-score vector
#' @param verbose   logical, \code{TRUE} for intermediate output
#' 
#' @return A data frame with columns `\code{quasi-deviance}`, `\code{Minor}`, `\code{det}`, `\code{max}`,
#'  `\code{trace}` and `\code{score_max}` (see vignette). The second column shows the leading minor of
#'   the observed QI matrix which is not positive definite. If so, then the corresponding parameter estimate
#'   cannot be consistent at least in theory. The rows show the corresponding values
#'   for each parameter `\code{par}`. Further, the results of function \code{\link{quasiDeviance}} for each
#'   parameter in `\code{x}` are returned as `attr(,"qdev")`. In case of any errors the corresponding (failed)
#'   results are returned as a named attribute `\code{X}`. 
#' 
#' @details Ony in case of the quasi-likelihood approach using the criterion function "\code{qle}" we can compare the
#'  observed quasi-information matrix with the expected one by some (scalar equivalent) criteria, see the vignette.
#'  We measure the dissimilarity of both matrices evaluated at approximate roots given in `\code{par}` for which these
#'  criteria should be neglectably small if they correspond to a plausible root. By the smallest value of any of these
#'  criteria we extract the best root.
#'  
#' @author M. Baaske
#' @rdname checkMultRoot
#' @export
checkMultRoot <- function(est, par = NULL, verbose = FALSE){			
   if(est$qsd$criterion != "qle")
	  stop("Consistency check of multiple roots only for criterion `qle`.")
   if(.isError(est))
	  stop("The estimation result from function `qle` has errors. Please check the argument `est`.")
  
   if(is.null(par))
	 par <- est$par   
   info <- attr(est,"optInfo")  			   
   QD <- try(quasiDeviance(par,est$qsd,
				Sigma=NULL,
				W=info$W,theta=info$theta,
				cvm=est$cvm,
				cl=cl,
				verbose=verbose),
		   silent=TRUE)
   
   if(.isError(QD)) {
	   msg <- paste0("Failed to get quasi-deviance: ","\n")
	   message(msg)
	   return(.qleError(message=msg,call=match.call(),error=QD))
   }
   xdim <- attr(est$qsd$qldata,"xdim")
   
   X <- lapply(QD,
		  function(qd) {
			  qIinv <- 
				  if(!.isError(qd)) {
					try(solve(qd$I),silent=TRUE)
				  } else NA	  	
			  
			  if(.isError(qIinv) || !is.numeric(qIinv) || anyNA(qIinv)) {		
				  warning("Failed to invert quasi-information.")
				  return (qIinv)
			  }
			  M <- qIinv%*%qd$Iobs
			  
			  c("quasi-deviance"=qd$val,
				"Minor"=.isPosDef(qd$Iobs),
		        "|det|"=abs(1-det(M)),
			    "|max|"=max(abs(M-diag(xdim))),
			    "|trace|"=abs(1-sum(diag(M))/xdim), 
			    "|score_max|"=max(abs(qd$score)))                						   
		  })
  
	isErr <- sapply(X,function(x) .isError(x) || !is.numeric(x) || anyNA(x))	
	M <- do.call(rbind,X[which(isErr == FALSE)])
	# check names
	row.names(M) <-
	 if(is.matrix(par) && !is.null(row.names(par)))
	  row.names(par)
     else if(is.list(par) && !is.null(names(par))) 		
	  names(par) 
     else paste0("par",1:nrow(M)) 
	
	structure(
	  as.data.frame(M),
	  "qdev" = QD,
	  "X" = if(any(isErr)) X[isErr] else NULL,
	  "error" = if(any(isErr)) isErr else NULL)		   	 
}


#' @name 	qleTest
#'
#' @title	Monte Carlo testing
#'
#' @description Monte Carlo hypothesis testing 
#'
#' @param est			object of class \code{qle}, the estimation results from function \code{\link{qle}}
#' @param local	    	optional, object of class \code{QSResult}, \code{NULL} (default), local estimation results 
#' @param sim			user supplied simulation function (see \code{\link{qle}})
#' @param ...			arguments passed to the simulation function `\code{sim}`
#' @param nsim			number of model replications to generate the simulated statistics
#' @param obs			optional, \code{NULL} (default), simulated statistics at the hypothesized parameter
#' @param check.root    logical, \code{FALSE} (default), whether to check consistency of estimated parameter (see \code{\link{checkMultRoot}})  						
#' @param na.rm 		logical, \code{TRUE}  (default), whether to remove `NA` values from the matrix of
#' 						re-estimated parameters
#' @param cl			cluster object, \code{NULL} (default), see \code{\link[parallel]{makeCluster}}
#' @param iseed			integer, the seed, \code{NULL} (default) for no seeding of the RNG stream for each worker
#' @param verbose   	logical, \code{TRUE} for intermediate output
#'
#' @return An object of class "\code{qleTest}" as a list of
#'  (predicted) standard errors as follows:
#' 
#' 	\item{param}{ data frame of estimated parameters and error measures}
#' 	\item{test}{ the test result}
#'  \item{Stest}{ name of the test}
#' 
#'   with attributes:
#' 	 \item{msem}{ mean square error matrix of re-estimated parameters}
#'   \item{aiqm}{ average inverse quasi-information matrix over all re-estimated parameters}
#'   \item{obs}{ list of simulated observed statistics}
#'   \item{optRes}{ results from re-estimating the model parameters for each simulated observation from `\code{obs}`}
#'	 \item{mean.score}{ average quasi-score or average gradient of MD at the re-estimated parameters}
#'   \item{info}{ list of indices of re-estimation results where the inversion of the quasi-information matrix failed,
#'       the re-estimated parameters have NA values, and criterion function minimizations have errors or did not converge,
#'       the integer seed `\code{iseed}`}
#' 
#'  @details 
#'  
#' 	The function expects an object of class \code{\link{qle}}. Simulated statistics which are already available at the estimated
#'  parameter can be passed by the argument `\code{obs}`. Otherwise the function first generates those realizations using `\code{nsim}`
#'  model replications. The criterion functions are used as specified in the object `\code{qsd}` and will not be further improved by
#'  additional samples during the test since this would result in a full estimation procedure again. The test statistic is either chosen
#'  as the current criterion function in `\code{OPT}` (see  \code{\link{getQLmodel}}, argument `\code{criterion}`) or taken from the
#'  optional argument `\code{local}` if supplied. Given optimization results by argument `\code{local}` of class \code{QSResult}, the
#'  user can select a different criterion function as a test statistic than used before for estimating the parameter. Apart from
#'  the QD as a test statistic, in principle, any supported type of least squares criterion (or more general Mahalanobis distance) can
#'  be used depending on the choice of variance matrix approximation (see \code{\link{covarTx}}). In some situations, the re-estimations
#'  might fail to converge due. However, the user can control the convergence of local solvers (including quasi-scoring) by the
#'  corresponding control parameters passed to the solvers, see \code{\link{searchMinimizer}}. Failed re-estimation results are extracted
#'  and stored in the attribute `\code{info}`. In addition, we return the standard error, predicted standard error (based on the average
#'  inverse quasi-information matrix), root mean square error, bias and sample mean value of the re-estimated parameters in order to assess
#'  the accuracy of the estimated parameter.    
#'
#' 	\subsection{Background}{
#' 
#'  The function tests the null hypothesis \eqn{H_0:\,\hat{\theta}=\theta_0}, that is, whether the statistical
#'  model is true, against the alternative \eqn{H_1:\,\hat{\theta}\neq\theta_0} by a Monte Carlo approach (see vignette).
#'  Due to the approximate nature of the assumed statistical model for the observed data the exact distribution of the test
#'  statistics (Mahalanobis distance or quasi-deviance) is generally unknown and therefore its asymptotic distribution might
#'  be an unrealistic assumption for the null hypothesis. For this reason, and in order to retrieve an empirical P-value for
#'  testing, we (artifically) generate new observations from the outcome of the model replications and re-estimate the model
#'  parameter for each realization in the same way as done before when estimating the model parameter. This includes all versions
#'  of variance approximations (by kriging or average approximations) and types of prediction variances (by kriging or a CV-based
#'  approach) used to construct the statistics for estimating the parameter. 
#' 
#' }
#' 
#' For an in depth example see the package vignette.
#' 
#' @author M. Baaske
#' @rdname qleTest
#' @export
qleTest <- function(est, local = NULL, sim, ...,
			 		 nsim = 100, obs = NULL, check.root=FALSE,
					  na.rm = TRUE, cl = NULL, iseed = NULL, verbose = FALSE)
{				  
	if(.isError(est))
	  stop("Estimation result has errors. Please see attribute `error`.")    
    # last evaluation of criterion function  	
	if(.isError(attr(est,"final")))
	  stop("Final criterion function evaluation has errors. Please check attribute `error`.")
		
    args <- list(...)
	# basic checks
	stopifnot(class(est) == "qle")
  	stopifnot(class(est$qsd)=="QLmodel")
	xdim <- attr(est$qsd$qldata,"xdim")
	
	# estimated parameter	 
	if(is.null(local)){		
		local <- attr(est,"final")
		if(.isError(local) || !attr(est,"optInfo")$minimized)
		 stop("Final optimization result has errors. Please check attribute `final` and `optInfo`.")
	 	else if( local$convergence < 0)
		  warning(paste0("Local search did not converge by method: ",local$method))
    } else {
       stopifnot(class(local)=="QSResult")	 
	   # set current test statistic
	   est$qsd$criterion <- local$criterion		
	}	
	# argument ids
	id <- pmatch(names(args),names(formals(searchMinimizer)))
	# check `searchMinimizer' args 
	opt.args <- args[which(!is.na(id))]	
	
	# check with qsd hidden	
	.checkfun(searchMinimizer,opt.args,
		hide.args=c("x0","qsd"),check.default=FALSE)
		
	# use last Sigma (unless kriging variance matrix)  
	if(est$qsd$var.type != "kriging"){
		info <- attr(est,"optInfo")
		opt.args <- c(opt.args,list(W=info$W,theta=info$theta)) 
  	}
	# use final (local) method by default
	# if not given as an input argument
	if(!("method" %in% names(opt.args)))
	 opt.args$method <- local$method
    # invert QI for predicted std. error (asymptotic) at estimated theta 
	qi <- try(solve(local$I),silent=TRUE)
	if(inherits(qi,"try-error") || anyNA(qi) )
	  qi <- NULL
	
	# MC simulation of observed `data`
	# if no observations supplied	
	if(is.null(obs)){		
		nid <- which(!is.na(id))
		if(length(nid) > 0L)
		 args <- args[-nid]		
		sim <- match.fun(sim)
		# check `sim` input values
		.checkfun(sim,args,remove=TRUE)
		stopifnot(is.numeric(nsim) || nsim > 0)
				
		simFun <- function(x) try(do.call(sim,c(list(x),args)))		
		sim.args <-
			list(sim=simFun,X=local$par,nsim=nsim,
				 mode="list",cl=cl,iseed=iseed,verbose=verbose)
		if(verbose)
		 cat("Simulating observed statistics...","\n")
 		
		obs <-
			tryCatch(
				do.call(simQLdata,sim.args),
				error = function(e) {
					msg <- paste0("Monte Carlo simulating observed statistics failed: ",
							 conditionMessage(e))
					message(msg) 
					.qleError(message=msg,call=match.call(),error=e)		   
				}
			)
		if(.isError(obs)) {			
			message("Simulating the observed statistics failed.")
			return(obs)
	 	}
		
	} else {
		if(.isError(obs))
		  stop("Argument `obs` has errors.")
		else if(class(obs) != "simQL" || !is.list(obs))
		  stop("Argument `obs` must be of class `simQL` and `list` type.")	   
	}
		
	# search minimizers	
	hasError <-  badInv <- 0
	if(verbose)
	 cat("Estimating parameters...","\n")
	
 	RES <- do.call(doInParallel,
			  c(list(X=obs[[1]],
					FUN=function(x,...)	{
						searchMinimizer(x0=local$par,qsd=est$qsd,...,    
							cvm=est$cvm,obs=x,info=TRUE,inverted=TRUE,
								check=FALSE,verbose=verbose)
					},
				cl=cl), opt.args))    		
	
	if(inherits(RES,"error")) {
	  msg <- paste0("Could not find MC replicated parameters: ",
			  conditionMessage(RES))
	  message(msg)
	  return(.qleError(message=msg,call=match.call(),error=RES))	
	}
	# check results again
	ok <- which(sapply(RES,function(x) !.isError(x) && x$convergence>0 ))
	if(length(ok) == 0L){
		msg <- paste("All re-estimations have errors or did not converge: ")
		message(msg)
		return(.qleError(message=msg,call=match.call(),error=RES))	   
	} else if(length(ok) < length(RES))
		warning("Errors in re-estimating the parameters. Check attribute `optInfo` and `info`.")
  
	# also check H_0?
	# because sampling MC theta
	# must be done under H0	
	aiqm <- NULL
	mScore <- NULL
	badInv <- integer(0)
    	  
	# average quasi-information matrix
	invI <- 
		lapply(RES[ok],
			function(x) {						
				try(solve(x$I),silent=TRUE) 
			})
	badInv <- ok[sapply(invI,function(x) inherits(x,"try-error") || anyNA(x))]
	if(length(badInv) == length(ok)) {
	  warning("Inversion of quasi-information matrices failed. Check attribute `info`.")
    } else {
	  if(length(badInv) > 0L)
		message(paste0("A total of ",length(badInv),
		 " inversions of quasi-information matrices failed. Check attribute `info`.")) 
	   
		# average matrix of inverse qi matrices
		aiqm <- matrix(
		  	      colMeans(
					do.call(rbind,
				 	 lapply(invI[ifelse(length(badInv) > 0L,-badInv,TRUE)],as.numeric)
					)
			 	  ),ncol=xdim)				
 	}
	# estimates	
	mpars <- t(sapply(RES[ok],"[[","par"))
	# mean zero distribution of score?	
	mScore <- t(sapply(RES[ok],"[[","score"))	
	has.na <- rowSums(is.na(cbind(mScore,mpars)))>0	
	if(na.rm && any(has.na)) {
		ok <- ok[-which(has.na)]
		mpars <- mpars[ok,,drop=FALSE]
		mScore <- mScore[ok,,drop=FALSE]		
  	}
	mScore <- try(colMeans(mScore),silent=TRUE)		
	if(nrow(mpars) <= 10L)
	 warning("Only a number of 10 or less parameters could be re-estimated.")	
    
 	# some (empirical) measures	
	Mse <- function(x) {(t(x)%*%x)/nrow(x)}
	msem <- Mse(mpars-do.call(rbind,rep(list(local$par),nrow(mpars))))
		
	# value of test statistic at re-estimated parameters			
	tvals <- sapply(RES[ok],"[[","value") 
	stopifnot(is.numeric(tvals))
	
	# get efficient score test (with MC parameters)
	B <- structure(
			 data.frame(
			  cbind("par"=local$par,
				    "se"=apply(mpars,2,sd),					
					"rmse"=diag(msem)^0.5,
					"bias"=colMeans(t(t(mpars)-local$par)),
					"mean"=colMeans(mpars))),
			"sb"=local$val, "Sb"=tvals,
			"test"=est$qsd$criterion)
	
	# had errors
	hasError <- which(!(1:length(RES) %in% ok))
	if(length(hasError) > 0L)	
	  message(paste0("A total of ",length(hasError)," re-estimations failed."))
  
    chk <- NULL
    if(check.root && est$qsd$criterion=="qle") {
		chk <- checkMultRoot(est,local$par,verbose=verbose)
		if(.isError(chk))
		 message(.makeMessage("Consistency check for the estimated model parameter failed."))		
	}
  
	# results
	structure(.qleTest(B),					# test results
		    msem=msem,						# mean square error matrix
	    	aiqm=aiqm,						# average inverse QI (re-estimated parameters)
			qi=qi,							# inverse QI at estimated theta
	    	obs=obs,						# (MC) observations
	    	optRes=RES,					    # optimization results
			mean.score=mScore,				# average score/gradient
			criterion=est$qsd$criterion,
			solInfo=chk,					# values of consistency criteria
	  info=list(badInv=badInv, 									# inversion errors
			  	hasNa=which(has.na), 							# indices of NA parameters 
		 		hasError=hasError,
				iseed=iseed),																
	 class=c("qleTest"), call=sys.call())	
}

# printing function
#' @name print.qleTest
#' 
#' @title print \code{qleTest} results 
#' 
#' @description print the results of call to \code{\link{qleTest}}
#' 
#' @param x      object of class `\code{qleTest}` from call to \code{\link{qleTest}}
#' @param pl	 not used yet
#' @param digits number of (significant) digits to display
#' @param ... 	 ignored, additional arguments
#' 
#' @rdname print.qleTest
#' @method print qleTest
#' @export 
print.qleTest <- function(x, pl = 1, digits = 5,...) {
	cat("\nCall:\n\n")
	cat(paste(deparse(attr(x,"call")), sep="\n", collapse = "\n"), "\n\n", sep="")
	# consistency check of solution
	chk <- attr(x,"solInfo")
	if(!is.null(chk)){
		cat("Consistency check - the smaller the better: \n\n")
		print(chk,print.gap=2,quote=FALSE)
		cat("\n\n")
		cat("----\n\n")
	}	
	cat("Coefficients:\n")	
	print(format(x$param, digits=digits),
			print.gap = 2, quote = FALSE)	
	cat("\n\n")
	cat(x$Stest,"\n\n")
	vals <- c(format(x$test[1], digits=digits),
			   formatC(signif(x$test[2], digits=digits), digits=digits,format="fg", flag="#"))
	names(vals) <- colnames(x$test)
	print(vals,print.gap = 2, quote = FALSE)
	
	cat("\n----\n\n")
	if(!is.null(attr(x,"mean.score"))) {
	  if(attr(x,"criterion") == "mahal")
		cat("Average gradient: \n\n")
	  else cat("Average quasi-score: \n\n")
	  print(format(attr(x,"mean.score"), digits=digits),
			  print.gap = 2, quote = FALSE)
	  cat("\n\n")
	  if(!is.null(attr(x,"aiqm")))
		  pse <- as.data.frame( cbind(diag(attr(x,"aiqm"))^0.5,diag(attr(x,"qi"))^0.5) )
	  	  dimnames(pse) <- list(row.names(x$param),c("Average","Estimate"))
		  cat("Predicted std. errors (asymptotic): \n\n")
		  print(format(pse, digits=digits),
				  print.gap = 2, quote = FALSE) 
    }	
	invisible(x)	
}