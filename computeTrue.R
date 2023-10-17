library(lavaan)
library(cubature)
library(mvtnorm)
library(dplyr)

# compute true CATE and ATE
trueCATE = function(tag, refZ){
  load(paste("param_", tag, ".RData", sep = ""))
  
  psi_inv = solve(psi)
  H_covar_inv = solve(H_covar)
  lambda_psi = t(lambda) %*% psi_inv
  
  M = lambda_psi %*% lambda + H_covar_inv
  d = lambda_psi %*% t(refZ-Z_intercept)
  M_inv = solve(M)
  
  CATE.true = C1 %*% M_inv %*% d + C0
  
  return(CATE.true)
}

trueATE = function(tag){
  load(paste("param_", tag, ".RData", sep = ""))
  
  ATE.true = C0
  
  return(ATE.true)
}