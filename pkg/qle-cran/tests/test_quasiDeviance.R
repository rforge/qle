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

# design, statistics, predictions
x0 <- c("mu"=2,"sigma"=1)
Xs <- as.matrix(qsd$qldata[c(1,2)])
Tstat <- qsd$qldata[c(3,4)]
pred <- estim(qsd$covT,x0,Xs,Tstat,krig.type="var")[[1]]
pred2 <- estim(qsd$covT,x0,Xs,Tstat,krig.type="dual")[[1]]

# compute quasi-deviance with use of kriging variances
D <- quasiDeviance(x0,qsd,verbose=TRUE)[[1]]

# quasi-score
S <- attr(D,"Sigma")
invS <- solve(S)
B <- invS%*%t(D$jac)

# prediction variances
# of sample mean of statistic Z=E[T(X)]
stopifnot(D$sig2==pred$sigma2)
D$score
(qs <- (qsd$obs-pred$mean)%*%B)

# quasi-information
D$I
t(B)%*%(S+diag(D$sig2))%*%B

# variance quasi-score vector
D$varS
(C <- t(B)%*%diag(D$sig2)%*%B)

# variance matrix of statistics Var(T(X))
print(S)
covarTx(qsd,theta=x0)

# value quasi-deviance
D$value
qs%*%solve(D$I)%*%t(qs)

# modified quasi-information: Mahalanobis distance of quasi-score
D$qval
qs%*%solve(C)%*%t(qs)
