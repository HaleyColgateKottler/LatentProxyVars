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
    errors[i] = cate-cate.true
  }
  return(errors)
}

errors = c()
nieve = c()
IPW = c()
match = c()
linear = c()


# dataGen(1, 6, 2, .5, 500, -.3, c(.2, -.1, .6), 0, c(1), c(.2), .6, c(.8), c(-.2), "1Da")
# betas
# dataGen(2,6,2,.5,500,.3,c(.5,.25, .3, -.1, -.2, .4),2, matrix(c(1,-.5,-.5,1), nrow = 2, byrow = TRUE), c(.7,.3),2,c(4,1), c(1,2),"2D")


for (i in 1:10){
  print(i)
  tag = "1Db"
  azGen(tag)
  k=1
  p=3
  
  rawData = dataImport(tag)
  
  ATE.true = trueATE(tag)
  
  nieve.ATE = nieveEst(rawData)
  nieve[i] = abs(ATE.true - nieve.ATE)
  linear.ATE = linearEst(rawData)
  linear[i] = abs(ATE.true - linear.ATE)
  IPW.ATE = IPWest(rawData)
  IPW[i] = abs(ATE.true - IPW.ATE)
  matching.ATE = matchingEst(rawData)
  match[i] = abs(ATE.true - matching.ATE)
  
  model = '
  efa("efa1")*h1 + efa("efa1")*h2 =~ V1+V2+V3+V4+V5+V6
  A ~ h1 + h2
  A | 0*t1
  A ~ 1
  
  h1 ~ 0*1
  h2 ~ 0*1
  '
  
  model <- '
  efa("efa1")*h1 =~ V1+V2+V3 
  A ~ h1
  A | 0*t1
  A ~ 1
  
  h1 ~ 0*1
  '
  
  params = fitUZA(model, rawData, k, p)
  expected.df = fitExpectations(params, rawData, k, p, 1)
  
  yModel = '
  Y ~ .*.
  '
  
  betas = fitMeanModel(yModel, subset(expected.df, select = c(Y, A, expectations1)))
  
  ATEest = ATE.est(expected.df, params, betas)
  
  errors[i] = abs(ATEest-ATE.true)
}

errs.df = data.frame(errors, IPW, linear, match, nieve)
errs.df
write.table(errs.df, file = paste("errors.csv", sep = ""), sep = ",")
hist(errors-IPW)

mean(errors)
mean(IPW)
mean(linear)
mean(match)
mean(nieve)
