# Copyright (C) 2017 Markus Baaske. All Rights Reserved.
# This code is published under the L-GPL.
#
# File: 	test_scoring.R
# Date:  	12/04/2017
# Author: 	Markus Baaske
# 
# Testing quasi Fisher-scoring iteration

library(qle)
data(normal)

qsd$var.type <- "cholMean"

# Scoring with average variance approximation
(S0 <- qscoring(qsd,
	x0=c("mu"=3.5,"sigma"=1.5), 
	opts=list("pl"=10,
			  "ftol_stop"=1e-17,					# stopping value
			  "ftol_abs"=1e-3,
			  "score_tol"=1e-6,
			  "grad_tol"=1e-3,
			  "slope_tol"=1e-7),					# test for local minimum by < ftol_abs
	 verbose=TRUE))

# BUG! finit precision issues with `qscoring` (see NEWS file)
x0 <- S0$par
# check with quasi-deviance
(D <- quasiDeviance(x0,qsd,verbose=TRUE)[[1]])
(S <- searchMinimizer(c("mu"=3.5,"sigma"=1.5),qsd,
		method=c("bobyqa"),verbose=TRUE))

xt <- S$par
qle::jacobian(qsd$covT,xt,X,Tstat)[[1]]$mean
S$jac
(Dt <- quasiDeviance(xt,qsd,verbose=TRUE)[[1]])
rbind(S$score,S0$score)

## check derivatives
#library(numDeriv)
#X <- as.matrix(qsd$qldata[,1:2])
#Tstat <- qsd$qldata[grep("mean.",names(qsd$qldata))]     
#p <- predictKM(qsd$covT,x0,X,Tstat)
#fn <-function(x,qsd, X, T) as.numeric(predictKM(qsd$covT,x,X,T))
#
## compare derivatives
#(J0 <- jacobian(fn, x0, method="simple",qsd=qsd, X=X, T=Tstat))
#(J2 <- t(qle::jacobian(qsd$covT,x0,X,Tstat)[[1]]$mean))
