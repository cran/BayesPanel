##########################################################################
## start-up and clean-up functions
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2010-present Siddhartha Chib and Chunhua Wu 
## 
##########################################################################

.onAttach <- function(...) {
 
   # echo output to screen
   cat("#\n# Bayesian Methods for Panel Data Modeling and Inference (BayesPanel)\n")
   cat("# Copyright (C) 2010-  Chunhua Wu and Siddhartha Chib \n")
   cat("#\n")
   require(Formula, quietly = TRUE)
   require(coda, quietly = TRUE)
   require(MASS, quietly = TRUE)
}

.onUnload <- function(libpath) {
    library.dynam.unload("BayesPanel", libpath)
}

