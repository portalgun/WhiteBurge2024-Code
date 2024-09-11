
library("gee")
setwd("/home/dambam/Documents/MATLAB/MVE/bin/DSP2/"


data <- read.csv("/Volumes/Data/.daveDB/Ext.txt")

formula <- R ~ D
# https://rdrr.io/rforge/alr/man/alr.html
# Gee = glm w/ correlation structure
#
# geepack
# geeglm
# gee4
#
# QIF
# GLS = special case of gee
# https://statmodeling.stat.columbia.edu/2006/12/27/generalized_est/
# repeated measures
mdl <-gee(
          formula,
          id=subj,
          data=data,
          family = binomial,
          # subset,     # what data to ommit
          # na.action,  # how to filter missing data
          # b = NULL,   # initial estimate
          tol = 0.001,
          maxiter = 25,
          silent = TRUE,
          corstr = "exchangeable",
          # Mv = 1,     # when corst=stat_M_dep non+stat_M_dep or AR-M
          # R = NULL,   # user specified correlation when corstr="fixed"
          contrasts = NULL,
          #
          scale.fix = FALSE,
          scale.value = 1,
          v4.4compat = FALSE
)
