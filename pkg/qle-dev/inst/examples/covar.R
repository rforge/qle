# Compare prediction variances to
# the ones obtained from cross-validation
library(qle)
data(normal)

## fit first covariance model
X <- as.matrix(qsd$qldata[,1:2])
T <- qsd$qldata["mean.T1"]
V <- qsd$qldata[["var.T1"]]
# get sample means of simulated statistics
Tstat <- qsd$qldata[grep("mean.",names(qsd$qldata))]

# sirf-2 covariance model (default)
# 'alpha' parameter is fixed and not reml estimated below,
# however, any parameter of the covariance model can be estimated
# by the reml method or excluded from covariance fitting	
cvm <- setCovModel(model="sirfk",  param = c("scale"=0.001,"alpha"=2.0),
		fixed.param = c("alpha"), nugget=1e-4, npoints = nrow(X), trend = 2,
		  var.sim = V)
  
# fitting by reml	
fit <- fitCov(list(cvm),X,T,verbose=TRUE)[[1]]

# original from the (already reml fitted) data
c(qsd$covT[[1]]$param,qsd$covT[[1]]$nugget)
# fitted again
(pnew <- c(fit$model$param,fit$model$nugget))
## test reml evaluation
val <- reml(list(cvm),pnew,T[1],X)[[1]]
stopifnot(attr(fit,"optres")$objective == val)

# set matern covariance for second model and fit  
T2 <- qsd$qldata["mean.T2"]
V2 <- qsd$qldata[["var.T2"]]
# covariance model for 2nd statistic
cm2 <- setCovModel(model="matern",
	    param = c("scale"=0.001,"nu"=2.5,"rho"=0.5),
		nugget = 1e-4, npoints = nrow(X), dim = ncol(X), var.sim = V2, trend = 2)
  
# 2nd statistic with Matern covariance
fit2 <- fitCov(list(cm2),X,T2,verbose=TRUE)[[1]]
# compare reml fits
c(fit$model$param,fit$model$nugget)
c(fit2$model$param,fit2$model$nugget)

# 2nd statistic with Matern covariance
qsd_new <- qsd
# replace 2nd covariance model (sirfk) by Matern 
qsd_new$covT <- structure(list(fit$model,fit2$model),
				  class="krige")

## Grid for MSE estimation
x <- seq(qsd$lower[1],qsd$upper[1],by=0.05)
y <- seq(qsd$lower[2],qsd$upper[2],by=0.05)
p <- as.matrix(expand.grid(x,y))
 
## Kriging MSE
## old fit kriging variances
kvar <- varKM(qsd$covT,p,X,Tstat)
## new fit with Matern, kriging variances
kvar_matern <- varKM(qsd_new$covT,p,X,Tstat)

## Empirical integrated MSE (by estimated kriging variances):
## Matern covariance suggests a better fit
colMeans(kvar)
colMeans(kvar_matern)

## show prediction variances 
## SIRF-k covariance
dev.new()
z1 <- matrix(kvar[,2],ncol=length(y))
plot(x = 0, y = 0, type = "n", xlim=,range(x), ylim=range(y),xlab = "", ylab = "")
contour(x, y, z1, col = "black", lty = "solid",
		nlevels = 50, add = TRUE,vfont = c("sans serif", "plain"))
try(points(X,pch=23,cex=0.8,bg="black"),silent=TRUE)

# Matern covariance
dev.new()
z2 <- matrix(kvar_matern[,2],ncol=length(y))
plot(x = 0, y = 0, type = "n", xlim=,range(x), ylim=range(y),xlab = "", ylab = "")
contour(x, y, z2, col = "black", lty = "solid",
		nlevels = 50, add = TRUE,vfont = c("sans serif", "plain"))
try(points(X,pch=23,cex=0.8,bg="black"),silent=TRUE)

## Cross-validation 
cvm <- prefitCV(qsd)
cv <- crossValTx(qsd, cvm, p, type = "cve")
## or with rmsd
# cv <- crossValTx(qsd, cvm, p, type = "rmsd")
colMeans(cv)
z3 <- matrix(cv[,2],ncol=length(y))

# show cross-validation based prediction errors
dev.new()
plot(x = 0, y = 0, type = "n", xlim=,range(x), ylim=range(y),xlab = "", ylab = "")
contour(x, y, z3, col = "black", lty = "solid",
		nlevels = 50, add = TRUE,vfont = c("sans serif", "plain"))
try(points(X,pch=23,cex=0.8,bg="black"),silent=TRUE)
