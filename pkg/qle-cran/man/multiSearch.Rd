% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qsOpt.R
\name{multiSearch}
\alias{multiSearch}
\title{A multistart version of local searches for parameter estimation}
\usage{
multiSearch(x0 = NULL, qsd, ..., nstart = 10, xstart = NULL,
  optInfo = FALSE, multi.start = FALSE, cores = getOption("mc.cores",
  0L), cl = NULL, pl = 0L, verbose = FALSE)
}
\arguments{
\item{x0}{numeric, \code{NULL} (default), list,  vector or matrix of starting parameters}

\item{qsd}{object of class \code{\link{QLmodel}}}

\item{...}{arguments passed to \code{\link{searchMinimizer}}}

\item{nstart}{number of random samples from which to start local searches (if `\code{x0}`=\code{NULL}, then ignored)}

\item{xstart}{list of starting points for local minimization or root finding by quasi-scoring}

\item{optInfo}{logical, \code{FALSE} (default), whether to store original local search results}

\item{multi.start}{logical, \code{FALSE} (default), whether to perform a multistart local search always otherwise only if first local search did not converge}

\item{cores}{integer, number of local CPU cores used, default is \code{options(mc.cores,0L)} and if \code{=0} then take all available cores}

\item{cl}{cluster object, \code{NULL} (default), of class \code{MPIcluster}, \code{SOCKcluster}, \code{cluster}}

\item{pl}{print level, use \code{pl}>0 to print intermediate results}

\item{verbose}{if \code{TRUE} (default), print intermediate output}
}
\value{
Object of class \code{QSResult} and attribute `\code{roots}`, i.e. the matrix of estimated parameters for which any of
 the available minimization methods has been successfully applied. If `code{optInfo}` is \code{TRUE}, then the originally estimtation reuslts
 are also returned. The best solution is stored as an attribute named `\code{par}` if found any.
}
\description{
The function is a multistart version of \code{\link{searchMinimizer}} which selects the best
	root of the quasi-score (if there is any) or a local minimum from all found minima according to the criteria described
 in the vignette.
}
\details{
The function performs a number of local searches depending which local method `\code{method}` was passed to
 \code{\link{searchMinimizer}}. Either the starting points are given by `\code{x0}` or are generated as an augmented 
 design based on the sample set stored in `\code{qsd}`. The function evaluates all found solutions and selects the one which 
 is best according to the criteria defined in the vignette. If none of the criteria match, then the parameter for which the lowest value
 of the criterion function was found is returned. Multistart searches can be done using a cluster object. Then for each generated/given obervation
 a number of \code{cores>1} multistart searches is performed in parallel if \code{fun="mclapply"} using the local cores of each cluster node.
}
\examples{
 data(normal)
 x0 <- c("mu"=3.5,"sigma"=1.5)
 S0 <- multiSearch(x0=x0,qsd,method=c("qscoring","bobyqa"),
           opts=list("ftol_stop"=1e-9,"score_tol"=1e-3),nstart=4,
            optInfo=TRUE,verbose=TRUE)

 roots <- attr(S0,"roots")
 id <- attr(roots,"id")
 stopifnot(!is.na(id)) 
 id  # index of best root found in matrix roots
 attr(roots,"par")  # the final parameter estimate w.r.t. id
 
}
\seealso{
\code{\link{checkMultRoot}}
}
\author{
M. Baaske
}
