# Copyright (C) 2017 Markus Baaske. All Rights Reserved.
# This code is published under the L-GPL.
#
# File: 	test_minimizer.R
# Date:  	12.04.2017
# Author: 	Markus Baaske
# 
# Test the criterion function minimization

library(qle)
data(normal)

searchMinimizer(c(2.5,1.5), qsd,verbose=TRUE) 
