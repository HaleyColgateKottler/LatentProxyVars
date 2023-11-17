library(lavaan)
library(cubature)
library(mvtnorm)
library(dplyr)

source("datagen.R")
source("latentMethod.R")

dataImport = function(tag){
  A = read.csv(paste("A_", tag, ".csv", sep = ""), header = FALSE)
  Z = read.csv(paste("Z_", tag, ".csv", sep = ""), header = FALSE)
  Y = read.csv(paste("Y_", tag, ".csv", sep = ""), header = FALSE)
  
  df = Z
  df$A = A$V1
  df$Y = Y$V1
  return(df)
}

CATE.errors = function(df, params, betas, tag){
  errors = c()
  
  Z = subset(df, select = -c(A, Y))
  for(i in 1:nrow(df)){
    z = Z[i,]
    cate = linearCATE(t(z), betas, params)
    cate.true = trueCATE(tag, z)
    errors[i] = cate-cate.true
  }
  return(errors)
}

tag = "2D"
k=2
p=6
sampSize = 10000
dataGen(2,6,2,.5,sampSize,.3,c(.5,.25, .3, -.1, -.2, .4),2, matrix(c(1,.5,.5,1), nrow = 2, byrow = TRUE), c(.7,.3),2,c(4,1), c(1,2),"2D")

d.ests = matrix(0,nrow = sampSize, ncol = p) 
d.trues = matrix(0,nrow = sampSize, ncol = p)
Md.ests = matrix(0,nrow = sampSize, ncol = k)
Md.trues = matrix(0,nrow = sampSize, ncol = k) 
bMdc.ests = matrix(0,nrow = sampSize, ncol = 1)
bMdc.trues = matrix(0,nrow = sampSize, ncol = 1)

azGen(tag)
rawData = dataImport(tag)
load(paste("param_", tag, ".RData", sep = ""))
  
# fit first model
model = '
  efa("efa1")*h1 + efa("efa1")*h2 =~ V1+V2+V3+V4+V5+V6
  A ~ h1 + h2
  A | 0*t1
  A ~ 1
  
  h1 ~ 0*1
  h2 ~ 0*1
  '
params = fitUZA(model, rawData, k, p)

# compute M true and M est
M.true = t(lambda)%*%solve(psi)%*%lambda + solve(H_covar)
Minv = solve(M.true)

M.est = t(params$lambda.est) %*% solve(params$psi.est) %*% params$lambda.est + solve(params$sigma.est)
M.inv.est = solve(M.est)

Z = subset(rawData, select = -c(Y,A))
for (i in 1:length(rawData$A)){
  z = Z[i,]
  
  # compute true denominator
  d.true = t(lambda) %*% solve(psi) %*% (t(z-Z_intercept))

  # compute est denominator
  d.est = t(params$lambda.est) %*% solve(params$psi.est) %*% (t(z - params$nu.est))
  
  
  d.ests[i,] = d.est 
  d.trues[i,] = d.true
  Md.ests[i,] = M.inv.est%*%d.est
  Md.trues[i,] = Minv%*%d.true 
  bMdc.ests[i] = t(params$b.est)%*%M.inv.est%*%d.est + params$c.est
  bMdc.trues[i] = t(B)%*%Minv%*%d.true + A_intercept
}

bMb1error = abs(t(B)%*%Minv%*%B+1 - (t(params$b.est)%*%M.inv.est%*%params$b.est+1))
Merror = norm(Minv - M.inv.est)
mseD = mean((d.trues - d.ests)^2)
mseMD = mean((Md.trues - Md.ests)^2)
msebMdc = mean((bMdc.trues - bMdc.ests)^2)

plot(d.trues[,2], abs(d.ests[,2]-d.trues[,2]))
plot(d.trues[,2], d.ests[,2])
