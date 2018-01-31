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

qscoring(qsd,
  x0=c("mu"=3.5,"sigma"=1.5), 
  opts=list("pl"=10, "slope_tol"=1e-7,"score_tol"=1e-4),
 verbose=TRUE)
