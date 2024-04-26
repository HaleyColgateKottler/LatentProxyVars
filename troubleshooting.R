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

tag = "2D1_29"
k=2
p=6
sampSize = 10000
dataGen(2,6,2,.5,sampSize,.3,c(.5,.25, .3, -.1, -.2, .4),2, matrix(c(1,0,0,1), nrow = 2, byrow = TRUE), c(.7,.3),2,c(4,1), c(1,2),tag)

num.trials = 500
c.ests = matrix(0,nrow = num.trials, ncol = 1) # should be .3
nu.ests = matrix(0,nrow=num.trials, ncol = p) # should be .5 .25 .3 -.1 -.2 .4
m.ests = matrix(0,nrow=num.trials, ncol = k*k)
sigma.ests = matrix(0,nrow=num.trials,ncol=k*k)
psi.ests = matrix(0,nrow=num.trials,ncol=p*p)
lambda.ests = matrix(0,nrow=num.trials,ncol=p*k)
b.ests = matrix(0,nrow=num.trials, ncol = k)

for (i in 1:num.trials){
  azGen(tag, sampSize)
  rawData = dataImport(tag)

  # fit first model
  model = '
    efa("efa1")*h1 + efa("efa1")*h2 =~ V1+V2+V3+V4+V5+V6
    A ~ h1 + h2
    A | 0*t1
    A ~ 1
    
    h1 ~ 0*1
    h2 ~ 0*1
    
    h1 ~~ 0*h2
    '

  params = fitUZA(model, rawData, k, p)
  c.ests[i] = params$c.est
  nu.ests[i,] = params$nu.est
  psi.ests[i,] = params$psi.est
  lambda.ests[i,] = params$lambda.est
  b.ests[i,] = params$b.est
  
  sigma.ests[i,] = params$sigma.est
  m.ests[i,] = t(params$lambda.est) %*% solve(params$psi.est) %*% params$lambda.est + solve(params$sigma.est)
}

load(paste("param_", tag, ".RData", sep = ""))

hist(c.ests)
mean(c.ests)
nu.errs = apply(nu.ests - Z_intercept, 1, mean)
hist(nu.errs)
mean(nu.errs) 

errs = t(apply(sigma.ests, 1, function(row) row - as.vector(H_covar)))
colMeans(errs)
hist(sigma.ests[,1], breaks=c(0,.5,1.1))
hist(sigma.ests[,2])
hist(sigma.ests[,4], breaks=c(0,.5,1.1))

# compute M true and M est
M.true = t(lambda)%*%solve(psi)%*%lambda + solve(H_covar)
m.errs = t(apply(m.ests, 1, function(row) row - M.true))
colMeans(m.errs)
colMeans(m.ests)

b.errs = t(apply(b.ests, 1, function(row) row - as.vector(B)))
hist(b.errs[,1])
hist(b.errs[,2])

psi.errs = t(apply(psi.ests, 1, function(row) row - as.vector(psi)))
colMeans(psi.errs)

lambda.errs = t(apply(lambda.ests, 1, function(row) row - as.vector(lambda)))
colMeans(lambda.errs)
