# Copyright (C) 2017 Markus Baaske. All Rights Reserved.
# This code is published under the L-GPL.
#
# File: 	test_quasiDeviance.R
# Date:  	27/07/2017
# Author: 	Markus Baaske
# 
# Testing quasi-deviance return values
# Shows how quasi-score, quasi-deviance, variances are computed

library(qle)
data(normal)

# type of variance matrix interpolation
qsd$var.type

# design, statistics, predictions
x0 <- c("mu"=2,"sigma"=1)
Xs <- as.matrix(qsd$qldata[c(1,2)])
Tstat <- qsd$qldata[c(3,4)]
pred <- estim(qsd$covT,x0,Xs,Tstat,krig.type="var")[[1]]

## dual kriging
# pred2 <- estim(qsd$covT,x0,Xs,Tstat,krig.type="dual")[[1]]

# compute quasi-deviance with use of kriging variances
D <- quasiDeviance(x0,qsd,verbose=TRUE)[[1]]

# remove kriging variances from variance matrix approximation
S <- attr(D,"Sigma") - diag(pred$sigma2)
invS <- solve(S)
B <- invS%*%t(D$jac)

# kriging prediction variances of sample means of statistics
stopifnot(D$sig2==pred$sigma2)

# quasi-score vector
D$score
(qs <- (qsd$obs-pred$mean)%*%B)

# quasi-information
D$I
t(B)%*%(S+diag(D$sig2))%*%B

# variance of quasi-score vector (only within the kriging model)
D$varS
(C <- t(B)%*%diag(D$sig2)%*%B)

# variance matrix of statistics Var(T(X))
print(S)
covarTx(qsd)			# no added kriging variances
covarTx(qsd,theta=x0)   # adding kriging variances at theta

# modified quasi-deviance value based on modified quasi-information matrix
D$value
qs%*%solve(D$I)%*%t(qs)

# modified quasi-deviance value based on kriging errot of quasi-score
D$qval
qs%*%solve(C)%*%t(qs)
