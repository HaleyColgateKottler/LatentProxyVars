library(lavaan)

sample_size = 500
p=6
k=2

H_mean = c(0,0)
H_covar = matrix(c(1,.2,.2,1), nrow=k)
lambda_common = matrix(c(1,1.25,1.5,2,0,0,0,.75,0,.5,-1,1), nrow=p)

epsilon_Z = MASS::mvrnorm(n = sample_size, mu = rep(0,p), Sigma = .2*diag(p))

H = MASS::mvrnorm(n=sample_size, mu = H_mean, Sigma = H_covar)
Z_set = matrix(data=NA, nrow = sample_size, ncol = p)
A_set = c()
B = matrix(c(3,2), nrow = 1)
var_A = .1

lambda_common %*% H[i,] + epsilon_Z[i,]

for (i in 1:sample_size){
  Z = lambda_common %*% H[i,] + epsilon_Z[i,]
  Z_set[i,] = t(Z)
  
  epsilon_A = rnorm(1, mean = 0, sd = sqrt(var_A))
  A_star = B%*%H[i,] + epsilon_A
  if (A_star >= 0){
    A = 1
  } else {
    A = 0
  }
  A_set[i] = A
}

meanZ = c(mean(Z_set[,1]), mean(Z_set[,2]), mean(Z_set[,3]), mean(Z_set[,4]), mean(Z_set[,5]), mean(Z_set[,6]))
meanZ
Z_set = sweep(Z_set,2,meanZ)
df = data.frame(Z=Z_set, A=A_set)

cov(df)
cor(df)

model <- '
efa("efa1")*h1 + 
efa("efa1")*h2 =~ Z.1 + Z.2 + Z.3 + Z.4 + Z.5 + Z.6

A ~ h1 + h2

A | 0*t1
'
fit <- cfa(model, data=df, ordered = "A")
summary(fit, standardized=TRUE)
inspect(fit, "partable")

model2 <- '
efa("efa")*h1 + efa("efa")*h2 =~ Z.1 + Z.2 + Z.3 + Z.4 + Z.5 + Z.6
'
fit2 <- cfa(model2, data=df)
summary(fit2, standardized=TRUE)
inspect(fit2, "partable")
