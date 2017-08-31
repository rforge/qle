# Copyright (C) 2017 Markus Baaske. All Rights Reserved.
# This code is published under the L-GPL.
#
# File: 	test_scoring.R
# Date:  	12.04.2017
# Author: 	Markus Baaske
# 
# Testing quasi Fisher-scoring iteration

library(qle)
data(normal)

# starting point
x0 <- c("mu"=2.5,"sigma"=1.5)

opts <- list("pl"=10,"ftol_stop"=1e-9,"score_tol"=1e-6)

# Scoring with average variance approximation
qscoring(qsd, x0, opts=opts, verbose=TRUE)
