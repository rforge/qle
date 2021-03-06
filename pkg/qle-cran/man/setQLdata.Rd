% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simQLData.R
\name{setQLdata}
\alias{setQLdata}
\title{Setup of quasi-likelihood data for estimation}
\usage{
setQLdata(runs, X = NULL, var.type = "wcholMean", Nb = 0,
  na.rm = TRUE, verbose = FALSE)
}
\arguments{
\item{runs}{list or matrix of simulation results obtained from \code{\link{simQLdata}}}

\item{X}{list or matrix of (design) points, i.e. model parameters}

\item{var.type}{name of variance matrix approximation type, "\code{wcholMean}" (default)}

\item{Nb}{number of bootstrap samples, \code{=0} (default, no bootstrap used), to be generated for kriging the variance matrix,
only if `\code{var.type}`=`\code{kriging}`}

\item{na.rm}{if \code{TRUE} (default), remove `NA` values from simulation results}

\item{verbose}{if \code{TRUE}, print intermediate output}
}
\value{
An object of class \code{QLdata} as a data frame with columns:
  \item{X}{ Model parameters (\code{n=1},...,\code{q}) }
  \item{mean}{ Results of simulation runs (\code{m=1},...,\code{p}) }
  \item{var}{ Simulation variances of statistics (\code{m=1},...,\code{p}) }
  \item{L}{ if applicable, Cholesky decomposed terms of variance matrices of statistics (k=1,...,(m*(m+1)/2))}
  \item{Lb}{ if applicable, bootstrap variances of covariances}
	where `\code{p}` denotes the number of user defined statistics and `\code{q}` the problem dimension, that is,
 the number of statistical model parameters to be estimated.	

	The following items are stored as attributes:

  \item{type}{ see above}
	 \item{nsim}{ number of simulations spent at each (design) point}
	 \item{xdim}{ length of model parameter}
	 \item{nWarnings}{ Number of warnings during simulation runs}
 \item{nErrors}{ Number of errors during simulation runs}
	 \item{nIgnored}{ List of parameters ignored due to errors}
}
\description{
Aggregate the data for quasi-likelihood estimation
}
\details{
The function aggregates all neccessary data for quasi-likelihood estimation storing the
	sample points and the corresponding simulation results of the statistics. If `\code{X}` equals \code{NULL},
 then the sample points are taken from the object `\code{runs}`.

 The most critical part is the decomposition of variance matrices for each sample point unless `\code{var.type}`
 equals "\code{const}" in which case a constant variance matrix approximation is expected to be given by the user in function \code{\link{qle}}.
 The Cholesky decompositions are used for an average approximation of the variance matrices of the statistics when calculating the
 quasi-score vector or any type of function criterion. If these fail for any reason we try to ignore, if possible, the corresponding sample
 points and exclude these from all subsequent computations. Unless a constant variance matrix estimate is used, the default is to approximate the
 variance matrix at any model parameter by either a kriging approximation of the \emph{Cholesky} terms (kriging the variance matrix) or as an average
 over all sampled variance matrices (variance matrix average approximation) also based on the decomposed Cholesky terms (see vignette).
}
\examples{
# simulate model statistics at LHS design
sim <- simQLdata(sim =
         function(x,cond) {
           X <- rlnorm(cond$n,x[1],x[2])
           c("MED"=median(X),"MAD"=mad(X))
         },
         cond=list("n"=10),
         nsim=10, N=10, method="maximinLHS",
         lb=c("mu"=-1.5,"sd"=0), ub=c("mu"=2,"sd"=1))

# setup the QL data model using defaults
qldata <- setQLdata(sim,verbose=TRUE) 
  
}
\author{
M. Baaske
}
