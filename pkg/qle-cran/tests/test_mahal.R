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
qsd$var.type <- "cholMean"
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
LQ <- mahalDist(theta,qsd,Sigma=diag(2))
(S <- attr(LQ,"Sigma")) # actually already inverted
Xs <- as.matrix(qsd$qldata[c(1,2)])
Tstat <- qsd$qldata[c(3,4)]
pred <- estim(qsd$covT,theta,Xs,Tstat)[[1]]
# criterion value
t(qsd$obs-pred$mean)%*%S%*%(qsd$obs-pred$mean)
LQ[[1]]$value

# same but now use Sigma with prediction variances
LQ2 <- quasiDeviance(theta,qsd,Sigma=diag(2))
LQ2[[1]]$value
attr(LQ2[[1]],"Sigma")