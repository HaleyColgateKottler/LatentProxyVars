# proximal causal inference
library(lavaan)
library(cubature)
library(mvtnorm)
library(dplyr)
library(pracma)
source("latentMethod.R")
source("datagen.R")

dataImport = function(tag){
  A = read.csv(paste("A_", tag, ".csv", sep = ""), header = FALSE)
  Z = read.csv(paste("Z_", tag, ".csv", sep = ""), header = FALSE)
  Y = read.csv(paste("Y_", tag, ".csv", sep = ""), header = FALSE)
  c 
  df = Z
  df$A = A$V1
  df$Y = Y$V1
  return(df)
}

# model A to get eta
tag = "1Dproxy"
nsamples = 1000
dataGen(1,2,2,.5,nsamples,.3,c(.5,.25),2, matrix(c(1), nrow = 1, byrow = TRUE), c(.7),3,c(4), c(2),tag)
df = dataImport(tag)

PA1 = sum(df$A == 1)/length(df$A)

# solve E[Y|Z,A=1] = alpha E[h(W)|V,A=1;eta]

df1 = df[df$A==1,]
regW1 = lm('V2~V1+1', df1)

df1$W0 = sapply(df1$V2,max,0)
regSpline1 = lm('W0~V1+1', df1)

df1$Linear = predict(regW1, df1)
df1$Spline = predict(regSpline1, df1)
regYV.1 = lm('Y~Linear*Spline', df1)

alpha = regYV.1$coefficients

df0 = df[which(df$A == 0),]
regW0 = lm('V2~V1+1', df0)
df0$W0 = sapply(df0$V2,max,0)
regSpline0 = lm('W0~V1+1', df0)

df0$Linear = predict(regW0, df0)
df0$Spline = predict(regSpline0, df0)
regYV.0 = lm('Y~Linear*Spline', df0)
gamma = regYV.0$coefficients

linear1.pred = predict(regW1, df)
spline1.pred = predict(regSpline1, df)
A1.expectation = matrix(c(rep(1,nsamples), linear1.pred, linear1.pred*spline1.pred), ncol = 3)

linear0.pred = predict(regW0, df)
spline0.pred = predict(regSpline0, df)
A0.expectation = matrix(c(rep(0,nsamples), linear0.pred, linear0.pred*spline0.pred), ncol = 3)

cate.est = rowSums((alpha[c(1,2,4)]-gamma[c(1,2,4)])*((PA1)* A1.expectation + (1-PA1)*A0.expectation))

source("computeTrue.R")

proximal.errors = c()
for (i in 1:length(df$A)){
  z = df[i,1:2]
  cate.true = trueCATE(tag, z)
  proximal.errors[i] = cate.true - cate.est[i]
}
hist(proximal.errors)
