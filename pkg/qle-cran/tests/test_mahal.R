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
mahalDist(theta,qsd,Sigma=diag(2))
