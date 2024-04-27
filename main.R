library(lavaan)
library(cubature)
library(mvtnorm)
library(dplyr)

source("datagen2.R")
# source("computeTrue.R")
# source("comparisonMethods.R")
source("latentMethod2.R")

dataImport = function(tag){
  A = read.csv(file.path("Data", paste("A_", tag, ".csv", sep = "")), header = FALSE)
  Z = read.csv(file.path("Data", paste("Z_", tag, ".csv", sep = "")), header = FALSE)
  Y = read.csv(file.path("Data", paste("Y_", tag, ".csv", sep = "")), header = FALSE)
  df = Z
  df$A = A$V1
  df$Y = Y$V1
  return(df)
}

# CATE.errors = function(df, params, betas, tag){
#   errors = c()
#   
#   Z = subset(df, select = -c(A, Y))
#   for(i in 1:nrow(df)){
#     z = Z[i,]
#     cate = linearCATE(t(z), betas, params)
#     cate.true = trueCATE(tag, z)
#     errors[i] = cate-cate.true
#   }
#   return(errors)
# }

tag = "1D"
savemarker = 50

for (sampleSize in c(1000)){
latent = c()
# nieve = c()
# IPW = c()
# match = c()
# linear = c()
# 2Da run of 50 trials with 1000 samples in 2D and ours majorly beats IPW
# dataGen(1, 3, 3, .5,1000, -.3, c(.2, -.1, .6), 0, c(1), c(.9), .6, c(.8), c(-.2), tag)
# betas
#dataGen(2,6,2,.5,1000,.3,c(.5,.25, .3, -.1, -.2, .4),2, matrix(c(1,-.5,-.5,1), nrow = 2, byrow = TRUE), c(.9,-.7),.8,c(.6,.5), c(.4,.2),tag)

for (j in 1:1000){
  i = j %% savemarker + 1
  
  azGen(tag, sampleSize)
  k=1
  p=3
  
  rawData = dataImport(tag)
  
  ATE.true = -.7
  
  # nieve.ATE = nieveEst(rawData)
  # nieve[i] = nieve.ATE
  # linear.ATE = linearEst(rawData)
  # linear[i] = linear.ATE
  # IPW.ATE = IPWest(rawData)
  # IPW[i] = IPW.ATE
  # matching.ATE = matchingEst(rawData)
  # match[i] = matching.ATE
  
  
  
  model = '
  efa("efa1")*h1 =~ V1+V2+V3+A
  
  h1 ~ 0*1
  '
  
  params = fitUZA(model, rawData, k, p)

  rawData$M = fitExpectations(params, rawData)
  
  AZ = subset(rawData, select = -c(Y, M))
  Z = subset(AZ, select = -c(A))
  
  yModel = '
  Y ~ M + A*M + A
  '
  betas = fitMeanModel(yModel, subset(rawData, select = c(Y, A, M)))
  
  ATEest = ATE.est(Z, params, betas)
  latent[i] = ATEest
  # 
  # latent.est = c()
  # latent.err = c()
  # cates.true = c()
  # for (i in 1:sampleSize){
  #   z = rawData[i,1:2]
  #   latent.est[i] = linearCATE(t(z), betas, params)
  #   cate.true = trueCATE(tag, z)
  #   latent.err[i] = cate.true - latent.est[i]
  #   cates.true[i] = cate.true
  # }
  # hist(latent.est)
  # plot(latent.est,cates.true)
  # abline(0,1)
  # # plot(latent.err, proximal.errors)
  # if (j %% savemarker == 0){
  #   print(sampleSize)
  #   print(j)
  #   errs.df = data.frame(latent, IPW, linear, match, nieve)
  #   write.table(errs.df, paste("errors_", tag, as.character(sampleSize), ".csv", sep = ""), sep = ",", append=TRUE,col.names=FALSE, row.names=FALSE)
  #   latent = c()
  #   nieve = c()
  #   IPW = c()
  #   match = c()
  #   linear = c()
  # }
}
}





mean(latent)
mean(IPW)
mean(linear)
mean(match)
mean(nieve)

