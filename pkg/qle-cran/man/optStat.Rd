% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/krige.R
\name{Subset of statistics}
\alias{Subset of statistics}
\alias{optStat}
\title{Optimal subset selection of statistics}
\usage{
optStat(theta, qsd, kmax = p, cumprop = 1, ..., cl = NULL,
  verbose = FALSE)
}
\arguments{
\item{theta}{list or matrix of points where to compute the criterion function
and to choose `\code{kmax}` statistics given the QL model `\code{qsd}`}

\item{qsd}{object of class \code{\link{QLmodel}}}

\item{kmax}{number of statistics to be selected with q <= \code{kmax} <= p}

\item{cumprop}{numeric vector either of length one (then replicated) or equal to the length of `\code{theta}` which sets the
proportions (0 < \code{cumprop} <= 1) of minimum overall contributions to each parameter component given the statistics}

\item{...}{further arguments passed to \code{\link{quasiDeviance}}}

\item{cl}{cluster object, \code{NULL} (default), of class \code{MPIcluster}, \code{SOCKcluster}, \code{cluster}}

\item{verbose}{logical, \code{FALSE} (default) for intermediate output}
}
\value{
A list which consists of 
	\item{id}{ indices of corresponding statistics}
	\item{Tmax}{ named vector of \code{id} corresponding maximum values of statistics}
 \item{cumprop}{ cumulated proportions of contributions of selected statistics to each of the parameter components} 
 \item{sorted}{ list of statistics (for each parameter) sorted in decreasing order of proportional contributions to the quasi-information}
}
\description{
The function finds a subset of at most \eqn{kmax <= p} statistics, where \code{p} is the number of available statistics
in the list `\code{qsd$covT}` (and at least of size equal to the length \code{q} of the parameter `\code{theta}`) and thus minimizes the expected
estimation error of the parameter when this subset is used for estimation. Based on the eigenvalue decomposition of the
variance-covariance matrix of the statistics this subset is chosen among all subsets of size at most equal to `\code{kmax}` or for
which all proportional contributions to each parameter component are greater than or equal to `\code{cumprop}` whatever happens first.

Since both matrices depend on `\code{theta}` so does the chosen subset of statistics. However, using a list of parameters as `\code{theta}`
returns a list of corresponding subsets. One can then easily choose the most frequent subset among all computed ones given either
a sample of parameters distributed over the whole parameter space or an appropriate smaller region, where, e.g., the
starting point is chosen from or the true model parameter is expected to lie in.
}
\examples{
 data(normal)
 # must select all statistics and thus using the
 # full information since we only have to statistics available 
 optStat(c("mu"=2,"sigma"=1),qsd,kmax=2)[[1]]

}
\author{
M. Baaske
}
