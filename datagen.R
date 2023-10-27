sigmaSim = function(k, p, communality){
  # p number of variables
  # k number of factors
  # communality 1 low, 2 wide, 3 high
  
  A = matrix(0, p, k)
  
  # need row entries to sum to k-1
  A[,1] = sample(0:(k-1), p, replace = TRUE)
  
  for (i in 1:p){
    if (k-2<1) break
    for (j in 2:(k-2)){
      current_sum = sum(A[i,])
      if (current_sum == k-1){
        A[i,j] = 0
      } else {
        A[i,j] = sample(0:(k-1-current_sum), 1, replace = TRUE)
      }
    }
  }
  A[,k] = (k-1)*matrix(1,1,p) - rowSums(A)
  
  # Add normal deviation
  c = .1*sample(7:9, k)
  x = rnorm(p*k)
  x1 = matrix(x, nrow = p)
  d = c(matrix(1,ncol=p, nrow = 1)/wordspace::rowNorms(x1))
  
  Y = A*c + d*x1*sqrt(1-c**2)
  
  # Apply skewing function
  Y2 = Y + abs(Y) + 0.2
  Y3 = abs(Y) + 0.2
  Z = (1.2/2.2) * (Y*Y2) / Y3
  
  g = rep(1,p)/wordspace::rowNorms(Z)
  
  # Scale to set communality
  if (communality == 1){
    B1 = diag(.1*sample(2:4, p, replace=TRUE), p)
  } else if (communality == 2){
    B1 = diag(.1*sample(2:8, p, replace=TRUE), p)
  } else if (communality == 3){
    B1 = diag(.1*sample(6:8, p, replace=TRUE), p)
  } else {
    B1 = matrix(0,p,p)
  }
  
  B2 = diag(p) - B1
  
  # Final factor loading matrix
  lambda = c(sqrt(B1)%*%g)*Z
  psi = sqrt(B2)
  
  sigma = lambda%*%t(lambda) + psi%*%t(psi)
  
  return(list(sigma, lambda, psi))
}

saveParams = function(k, p, communality, var_Y, n, A_intercept, Z_intercept, Y_intercept, H_covar, B, C0, C1, C2, tag){
  
  output = sigmaSim(k, p, communality)
  sigma = output[[1]]
  lambda = output[[2]]
  psi = output[[3]]
  save(k, p, communality, var_Y, n, sigma, lambda, psi, H_covar, B, C0, C1, C2, A_intercept, Z_intercept, Y_intercept, file = paste("param_", tag, ".RData", sep = ""))
}

azGen = function(tag){
  load(paste("param_", tag, ".RData", sep = ""))
  H_set = matrix(, nrow = n, ncol = k)
  Z_set = matrix(, nrow = n, ncol = p)
  A_set = c()
  Y_set = c()
  
  for (i in 1:n){
    epsilon_Z = MASS::mvrnorm(n = 1, mu = rep(0,p), Sigma = psi)
    epsilon_A = rnorm(1, mean = 0, sd = 1)
    epsilon_Y = rnorm(1, mean = 0, sd = sqrt(var_Y))
    
    H = MASS::mvrnorm(n=1, mu = rep(0,k), Sigma = H_covar)
    Z = lambda %*% H + Z_intercept + epsilon_Z
    A_star = B%*%H + A_intercept + epsilon_A
    if (A_star >= 0){
      A = 1
    } else {
      A = 0
    }
    Y = A * (C1 %*% H) + (C2 %*% H) + C0 * A + Y_intercept + epsilon_Y
    
    H_set[i,] = t(H)
    Z_set[i,] = Z
    A_set[i] = A
    Y_set[i] = Y
  }
  
  write.table(H_set, file = paste("H_", tag, ".csv", sep = ""), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Z_set, file = paste("Z_", tag, ".csv", sep = ""), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(A_set, file = paste("A_", tag, ".csv", sep = ""), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Y_set, file = paste("Y_", tag, ".csv", sep = ""), row.names = FALSE, col.names = FALSE, sep = ",")
}

dataGen = function(k, p, communality, var_Y, n, A_intercept, Z_intercept, Y_intercept, H_covar, B, C0, C1, C2, tag){
  saveParams(k, p, communality, var_Y, n, A_intercept, Z_intercept, Y_intercept, H_covar, B, C0, C1, C2, tag)
  azGen(tag)
}

dataGen(2,6,2,.5,1000,.3,c(.5,.25, .3, -.1, -.2, .4),2, matrix(c(1,.5,.5,1), nrow = 2, byrow = TRUE), c(.7,.3),2,c(4,1), c(1,2),"2D")

dataGen(1,6,3, .25, 500, -.2, c(.2, -.1, .4), 2, c(1), c(.9), -.7, c(-.2), c(.8), "1Db")
