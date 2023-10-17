# comparison methods

nieveEst = function(df){
  A1 = which(df$A == 1)
  A0 = which(df$A == 0)
  
  Y1 = mean(unlist(df$Y[A1]))
  Y0 = mean(unlist(df$Y[A0]))
  
  nieve.ATE = Y1-Y0
  return(nieve.ATE)
}

linearEst = function(df){
  m1 = lm(Y~.,df)
  linear.ATE = m1$coefficients["A"]
  
  return(linear.ATE)
}

IPWest = function(df){
  AZ = subset(df, select = -c(Y))
  Z = subset(AZ, select = -c(A))
  
  m1 = glm(formula = A ~ ., family = binomial(link = "logit"), data = AZ)
  x = predict(m1, newdata = Z, type="response")
  
  EY1 = mean(unlist(df$A*df$Y/x))
  EY0 = mean(unlist((1-df$A)*df$Y/(1-x)))
  
  IPW.ATE = EY1-EY0
  
  return(IPW.ATE)
}

matchingEst = function(df){
  Z = subset(df, select = -c(A,Y))
  z.mat = as.matrix(Z)
  
  A1_indices = which(df$A==1)
  A0_indices = which(df$A==0)
  
  total.diff1 = 0
  
  for (i in A1_indices){
    diffs = wordspace::rowNorms(sweep(z.mat[A0_indices,], 2, z.mat[i,]))
    match = which.min(diffs)
    total.diff1 = total.diff1 + df$Y[i] - df$Y[A0_indices[match]]
  }
  
  est1 = total.diff1 / length(A1_indices)
  
  total.diff2 = 0
  
  for (i in A0_indices){
    diffs = wordspace::rowNorms(sweep(z.mat[A1_indices,], 2, z.mat[i,]))
    match = which.min(diffs)
    total.diff2 = total.diff2 - df$Y[i] + df$Y[A1_indices[match]]
  }
  
  est2 = total.diff2 / length(A0_indices)
  
  matching.ATE = (total.diff1 + total.diff2) / length(df$A)
  
  return(matching.ATE)
}
