\dontrun{
# Quasi-Scoring iteration with and
# without using Kriging variances
library(qle)
data(normal)

# starting point
x0 <- c("mu"=3.5,"sigma"=0.5)
opts <- list("pl"=10, "ftol_stop"=1e-9, "score_tol"=1e-6)
	 
## Scoring with added Kriging variances
QSE <- qscoring(qsd, x0, opts=opts, verbose=TRUE)
print(QSE)

## Scoring with avergae variance matrix approximation
## and without using prediction variances (dual kriging)
# qsd$krig.type <- "dual"
# qsd$var.type <- "cholMean"
# QSE0 <- qscoring(qsd, x0, opts=opts, verbose=TRUE)

# check results
stopifnot(QSE$status>0 && QSE$convergence)

## Check with quasi-deviance
QDtheta <- quasiDeviance(QSE$par, qsd)[[1]]
stopifnot(identical(QDtheta$score,QSE$score))
}
