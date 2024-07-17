library(lavaan)
library(cubature)
library(mvtnorm)
library(dplyr)
library(reshape2)
library(ggplot2)
source("datagen2.R")

data.import <- function(tag, n) {
  A <- read.csv(
    file.path(
      "Data", "Samples",
      paste("A_", tag, "_", n, ".csv", sep = "")
    ),
    header = FALSE
  )
  Z <- read.csv(
    file.path(
      "Data", "Samples",
      paste("Z_", tag, "_", n, ".csv", sep = "")
    ),
    header = FALSE
  )
  Y <- read.csv(
    file.path(
      "Data", "Samples",
      paste("Y_", tag, "_", n, ".csv", sep = "")
    ),
    header = FALSE
  )
  df <- Z
  df$A <- A$V1
  df$Y <- Y$V1
  return(df)
}

tag = "2D"
p=6
k=2
constraints = matrix(1, nrow = p+1, ncol = k)

calcCATE_SE <- function(data, constraints, EM_output, a.dif, z){
  n = nrow(data)
  p = nrow(constraints)
  k = ncol(constraints)
  lambda.parms = sum(constraints==1)
  
  theta.EM = as.numeric(c(
    EM_output$mu,
    EM_output$lambda[constraints==1],
    EM_output$psi,
    EM_output$alpha,
    EM_output$gamma,
    EM_output$beta
  ))

  # estimating function
  # if weights are supplied, then it uses REM estimating equations; otherwise uses EM estimating equations
  estFUN <- function(data, constraints, a.dif, z){
    mu = apply(data[,1:p], 2, mean)
    n = nrow(data)
    p = ncol(data) - 1
    k = ncol(constraints)
    lambda = matrix(0, nrow = p, ncol = k)
    lambda.idx = which(constraints==1)
    lambda.parms = length(lambda.idx)
    d.lambda = matrix(0, nrow = n, ncol = lambda.parms)
    d.psi = matrix(0, nrow = n, ncol = p)
    d.alpha = matrix(0, nrow = n, ncol = k+1)
    d.gamma = matrix(0, nrow = n, ncol = k+1)

    function(theta){
      mu = theta[1:p]
      lambda.est = theta[(p+1):(p+lambda.parms)]
      lambda[lambda.idx] <- lambda.est
      psi = theta[(p+lambda.parms+1):(p+lambda.parms+p)]
      alpha = theta[(2*p+lambda.parms+1):(2*p+lambda.parms+k+1)]
      gamma = theta[(2*p+lambda.parms+k+2):(2*p+lambda.parms+2*k+2)]
      beta = theta[2*p+lambda.parms+2*k+3]
      # estimating equations
      Z = apply(data[,1:p], 1, function(y) y - mu)
      sigma = lambda %*% t(lambda) + diag(psi)
      inv.sigma = solve(sigma)
      m = t(lambda) %*% solve(diag(psi)) %*% lambda + diag(k)
      f.base = solve(m)%*%t(lambda)%*%solve(diag(psi))
      d.mu = t(inv.sigma %*% Z)
      for (i in 1:n){
        D = inv.sigma %*% (diag(p) - (Z[,i] %*% t(Z[,i]) %*% inv.sigma))
        d.lambda[i,] = -(D %*% lambda)[lambda.idx]
        d.psi[i,] = -0.5 * diag(D)
        f = c(1, f.base%*%Z[,i])
        d.alpha[i,] = 2*(data[i,p+1] - (alpha+gamma*data[i,p-1])%*%f) %*% t(f)
        d.gamma[i,] = data[i,p]*2*(data[i,p+1] - (alpha+gamma*data[i,p])%*%f) %*% t(f)
      }
      small.m = t(lambda[1:(p-1), 1:k]) %*% solve(diag(psi[1:(p-1)])) %*% lambda[1:(p-1), 1:k] + diag(k)
      U.mean = solve(small.m) %*% t(lambda[1:(p-1), 1:k]) %*% solve(diag(psi[1:(p-1)])) %*% (z-mu[1:(p-1)])
      d.beta = a.dif * gamma %*% c(1, U.mean) - beta

      # score function for each person and parameter
      c(d.mu, d.lambda, d.psi, d.alpha, d.gamma, d.beta)
    }
  }
  
  # use geex package to estimate standard errors
  m.results.EM <- geex::m_estimate(
    estFUN = estFUN,
    data = data,
    roots = theta.EM,
    compute_roots = FALSE,
    outer_args = list(constraints = constraints,
                      a.dif = a.dif,
                      z = z))
  
  se = sqrt(diag(geex::vcov(m.results.EM)))
  names(se) <- c(paste0('mu', 1:p),
      paste0('lambda', which(constraints==1)),
      paste0('psi', 1:p),
      paste0('alpha', 1:(k+1)),
      paste0('gamma', 1:(k+1)),
      'beta'
    )
  
  return(se)
}

calcSE <- function(data, constraints, EM_output, a.dif){
  n = nrow(data)
  p = nrow(constraints)
  k = ncol(constraints)
  lambda.parms = sum(constraints==1)
  
  theta.EM = as.numeric(c(
    EM_output$mu,
    EM_output$lambda[constraints==1],
    EM_output$psi,
    EM_output$alpha,
    EM_output$gamma,
    EM_output$beta
  ))
  
  # estimating function
  # if weights are supplied, then it uses REM estimating equations; otherwise uses EM estimating equations
  estFUN <- function(data, constraints, a.dif){
    mu = apply(data[,1:p], 2, mean)
    n = nrow(data)
    p = ncol(data) - 1
    k = ncol(constraints)
    lambda = matrix(0, nrow = p, ncol = k)
    lambda.idx = which(constraints==1)
    lambda.parms = length(lambda.idx)
    d.lambda = matrix(0, nrow = n, ncol = lambda.parms)
    d.psi = matrix(0, nrow = n, ncol = p)
    d.alpha = matrix(0, nrow = n, ncol = k+1)
    d.gamma = matrix(0, nrow = n, ncol = k+1)
    
    function(theta){
      mu = theta[1:p]
      lambda.est = theta[(p+1):(p+lambda.parms)]
      lambda[lambda.idx] <- lambda.est
      psi = theta[(p+lambda.parms+1):(p+lambda.parms+p)]
      alpha = theta[(2*p+lambda.parms+1):(2*p+lambda.parms+k+1)]
      gamma = theta[(2*p+lambda.parms+k+2):(2*p+lambda.parms+2*k+2)]
      beta = theta[2*p+lambda.parms+2*k+3]
      # estimating equations
      Z = apply(data[,1:p], 1, function(y) y - mu)
      sigma = lambda %*% t(lambda) + diag(psi)
      inv.sigma = solve(sigma)
      m = t(lambda) %*% solve(diag(psi)) %*% lambda + diag(k)
      f.base = solve(m)%*%t(lambda)%*%solve(diag(psi))
      d.mu = t(inv.sigma %*% Z)
      for (i in 1:n){
        D = inv.sigma %*% (diag(p) - (Z[,i] %*% t(Z[,i]) %*% inv.sigma))
        d.lambda[i,] = -(D %*% lambda)[lambda.idx]
        d.psi[i,] = -0.5 * diag(D)
        f = c(1, f.base%*%Z[,i])
        d.alpha[i,] = 2*(data[i,p+1] - (alpha+gamma*data[i,p-1])%*%f) %*% t(f)
        d.gamma[i,] = data[i,p]*2*(data[i,p+1] - (alpha+gamma*data[i,p])%*%f) %*% t(f)
      }
      d.beta = a.dif * gamma[1] - beta
      
      # score function for each person and parameter
      c(d.mu, d.lambda, d.psi, d.alpha, d.gamma, d.beta)
    }
  }
  
  # use geex package to estimate standard errors
  m.results.EM <- geex::m_estimate(
    estFUN = estFUN,
    data = data,
    roots = theta.EM,
    compute_roots = FALSE,
    outer_args = list(constraints = constraints,
                      a.dif = a.dif))
  
  se = sqrt(diag(geex::vcov(m.results.EM)))
  names(se) <- c(paste0('mu', 1:p),
                 paste0('lambda', which(constraints==1)),
                 paste0('psi', 1:p),
                 paste0('alpha', 1:(k+1)),
                 paste0('gamma', 1:(k+1)),
                 'beta'
  )
  
  return(se)
}

load(file.path("Data", "Parameters", paste("param_", tag, ".RData", sep = "")))
se.est = calcSE(raw.data, constraints, EM_output, 1)

ests = matrix(0, nrow = 10, ncol = 35)
for (sample.size in 1:10){
  print(sample.size)
  azGen(tag, sample.size * 100)
  raw.data = data.import(tag, sample.size * 1000)
  EM_output = list("mu" = c(Z_intercept, A_intercept),
                   "lambda" = lambda,
                   "psi" = diag(psi),
                   "phi" = diag(2),
                   "alpha" = alpha,
                   "gamma" = gamma,
                   "beta" = .47)
  
  se.est = calcSE(raw.data, constraints, EM_output, 1, rep(1, p))
  ests[sample.size,] = se.est
}

df = melt(ests[,35])
df$Var1 = factor(df$Var1)
ggplot(df, aes(x = Var1, y = value, group = Var2, color = Var2)) + geom_line()
