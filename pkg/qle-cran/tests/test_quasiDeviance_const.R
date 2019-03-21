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
qsd$var.type <- "const"
Sigma <- diag(1,2) 

# design, statistics, predictions
x0 <- c("mu"=2,"sigma"=1)
Xs <- as.matrix(qsd$qldata[c(1,2)])
Tstat <- qsd$qldata[c(3,4)]
pred <- estim(qsd$covT,x0,Xs,Tstat,krig.type="var")[[1]]

# compute quasi-deviance with use of kriging variances
X <- rbind(x0,x0+c(0.1,0.2),x0+c(1,1))
quasiDeviance(X,qsd,Sigma,value.only=2,verbose=TRUE)
quasiDeviance(X,qsd,Sigma,value.only=3,verbose=TRUE)

D <- quasiDeviance(X,qsd,Sigma=Sigma,verbose=TRUE)
# both are the same as above
sapply(D,function(x) sum(diag(x$I)))
sapply(D,function(x) sum(diag(x$varS)))

# use only first point next
D <- D[[1]]
S <- attr(D,"Sigma") #- diag(pred$sigma2)
invS <- solve(S)
B <- invS%*%t(D$jac)

# kriging prediction variances of sample means of statistics
stopifnot(D$sig2==pred$sigma2)

# quasi-score vector
D$score
(qs <- (qsd$obs-pred$mean)%*%B)

# (original) quasi-information equals modified quasi-information
D$I
(D$jac)%*%invS%*%t(D$jac)
(C <- t(B)%*%(S)%*%B)
D$varS

# variance matrix of statistics cannot be computed here
print(S)
# covarTx(qsd)	

# modified quasi-deviance value based on modified quasi-information matrix
D$value
qs%*%solve(D$varS)%*%t(qs)

# return the value of the modified quasi-deviance only 
quasiDeviance(x0,qsd,Sigma=Sigma,verbose=TRUE,value.only=TRUE)
