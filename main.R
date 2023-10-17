library(lavaan)
library(cubature)
library(mvtnorm)
library(dplyr)

source("datagen.R")
source("computeTrue.R")
source("comparisonMethods.R")
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
    errors[i] = abs(cate-cate.true)
  }
  return(errors)
}

tag = "2D"
k=2
p=6

rawData = dataImport(tag)

ATE.true = trueATE(tag)

nieve.ATE = nieveEst(rawData)
linear.ATE = linearEst(rawData)
IPW.ATE = IPWest(rawData)
matching.ATE = matchingEst(rawData)

model <- '
efa("efa1")*h1 =~ V1+V2+V3
A ~ h1
A | 0*t1
A ~ 1

h1 ~ 0*1
'
model = '
efa("efa1")*h1 + efa("efa1")*h2 =~ V1+V2+V3+V4+V5+V6
A ~ h1 + h2
A | 0*t1
A ~ 1

h1 ~ 0*1
h2 ~ 0*1
'

params = fitUZA(model, rawData, k, p)
expected.df = fitExpectations(params, rawData, k, p, 1)

yModel = '
Y ~ .*.
'

betas = fitMeanModel(yModel, subset(expected.df, select = c(Y, A, expectations1)))
betas
ATE.est = ATE.est(expected.df, params, betas)

errors = CATE.errors(rawData, params, betas, tag)

mean(errors[rawData$A == 0])
mean(errors[rawData$A == 1])
