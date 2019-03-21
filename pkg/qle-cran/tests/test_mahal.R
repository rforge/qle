# Copyright (C) 2017 Markus Baaske. All Rights Reserved.
# This code is published under the L-GPL.
#
# File: 	test_mahal.R
# Date:  	12/04/2017
# Author: 	Markus Baaske
# 
# 	Test Mahalanobis distance with weighting matrix  
#	and average approximation of variance matrix

library(qle)
data(normal)

# options set when during initilization of `qsd`
qsd$var.type <- "wcholMean"
qsd$criterion <- "mahal"

# some parameters of the statistical model for evaluation
theta <- c("mu"=2,"sigma"=0.95) 

# compare criterion functions
MD <- mahalDist(theta,qsd)
QD <- quasiDeviance(theta,qsd)

# must be equal due to the same
# number of parameters (q) and statistics (p)
all.equal(MD[[1]]$value,QD[[1]]$value)

# least-squares (constant variance)
qsd$var.type <- "const"
qsd$criterion <- "mahal"
LQ <- mahalDist(theta,qsd,Sigma=diag(2))[[1]]
(S <- attr(LQ,"Sigma")) # actually already inverted
Xs <- as.matrix(qsd$qldata[c(1,2)])
Tstat <- qsd$qldata[c(3,4)]
pred <- estim(qsd$covT,theta,Xs,Tstat)[[1]]

# criterion value
t(qsd$obs-pred$mean)%*%S%*%(qsd$obs-pred$mean)
LQ$value

# same as above with constant Sigma (p==q)
LQ2 <- quasiDeviance(theta,qsd,Sigma=diag(2))[[1]]
LQ2$value
attr(LQ2,"Sigma")

# return the value of the sampling criterion 'logdet'. 
# A combination of minimizing estimation error (first term)
# pf model parameter 'theta' and prediction error of the
# quasi-score vector (second term)
crit <- function(qd,w=0.5) {													
	B <- solve(attr(qd,"Sigma"))%*%t(qd$jac)												
	varS <- t(B)%*%(attr(qd,"Sigma")+diag(qd$sig2))%*%B
	w*log(det(varS))+(1-w)*qd$value
}

# weighted average approximation
crit(MD[[1]],w=.5)
qsd$var.type <- "wcholMean"
mahalDist(theta,qsd,w=0.5,verbose=TRUE,value.only=2L)

# using a constant Sigma
qsd$var.type <- "const"
qd <- mahalDist(theta,qsd,Sigma=diag(2),w=0.5,verbose=TRUE,value.only=0L)[[1]]
qd$sig2 <- rep(0,2)
crit(qd,w=.5)
mahalDist(theta,qsd,Sigma=diag(2),w=0.5,verbose=TRUE,value.only=2L)
#quasiDeviance(theta,qsd,Sigma=diag(2),w=0.5,verbose=TRUE,value.only=2L)

# first term
crit(MD[[1]],w=1.0)
log(det(MD[[1]]$varS))
# second term
MD[[1]]$value
crit(MD[[1]],w=0.0)
