library(lavaan)
library(cubature)
library(mvtnorm)
library(dplyr)

tag = "2Da"
load(paste("param_", tag, ".RData", sep = ""))

probAis1function = function(u){
  return(pnorm(b%*%u+A_intercept)*dmvnorm(u, sigma = H_covar))
}
probAis1 = cubintegrate(probAis1function, rep(-Inf,k), rep(Inf,k), fDim=1)$integral

expectUa1function = function(u){
  return(u%*%(dnorm(b%*%u+A_intercept)*dmvnorm(u, sigma = H_covar)))
}
m1 = cubintegrate(expectUa1function, rep(-Inf,k), rep(Inf,k), fDim=k)$integral/probAis1

IPWest.goal = C0 + (C1+C2) %*% m1 + (probAis1/(1-probAis1))*C2%*%m1

source("comparisonMethods.R")
source("datagen.R")
dataImport = function(tag){
  A = read.csv(paste("A_", tag, ".csv", sep = ""), header = FALSE)
  Z = read.csv(paste("Z_", tag, ".csv", sep = ""), header = FALSE)
  Y = read.csv(paste("Y_", tag, ".csv", sep = ""), header = FALSE)
  
  df = Z
  df$A = A$V1
  df$Y = Y$V1
  return(df)
}

ipw.ests = c()
for (run in 1:100){
  print(run)
  azGen(tag)
  df = dataImport(tag)
  ipw.ests[run] = IPWest(df)
}
hist(ipw.ests)
IPWest.goal
