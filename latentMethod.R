fitUZA = function(model, df, k, p){
  smallerDF = subset(df, select = -c(Y))
  smallerDF$A = ordered(smallerDF$A, levels = c(0,1))
  
  fit <- sem(model, data=smallerDF, ordered = "A", fixed.x = c(h1 ~~ 0*h2), rotation = "varimax")
  
  parameters = list(
    lambda.est = inspect(fit,what="est")$lambda[1:p,1:k],
    psi.est = inspect(fit,what="est")$theta[1:p,1:p],
    sigma.est = inspect(fit,what="est")$psi[1:k,1:k],
    b.est = inspect(fit,what="est")$beta[k+1,1:k],
    nu.est = inspect(fit,what="est")$nu[1:p,1],
    c.est = inspect(fit,what="est")$alpha[k+1]
  )
  return(parameters)
}

fitExpectations = function(params, df, k, p, power){
  psi.est.inv = solve(params$psi.est)
  M.est = t(params$lambda.est) %*% psi.est.inv %*% params$lambda.est + solve(params$sigma.est)
  M.inv.est = solve(M.est)
  
  uCDF.generate = function(M.inv.est, d.est, b, c, power){
    uCDF = function(u){
      return(u^power*c(pnorm(u%*%b+c))*dmvnorm(u, M.inv.est%*%d.est, matrix(M.inv.est, nrow = k)))
    }
    return(uCDF)
  }
  
  constant.denom = sqrt(t(params$b.est)%*%M.inv.est%*%params$b.est+1)
  
  Z = subset(df, select = -c(Y,A))
  A = subset(df, select = c(A))
  
  expectations = matrix(, nrow = nrow(df), ncol = k)
  for(i in 1:nrow(df)){
    z <- Z[i,]
    a <- A[i,]
    d.est = t(params$lambda.est) %*% psi.est.inv %*% (t(z - params$nu.est))
    meanF = uCDF.generate(M.inv.est, d.est, params$b.est, params$c.est, power)
    inner.constant = (t(params$b.est)%*%M.inv.est%*%d.est+params$c.est)/constant.denom
    constant = c(a/pnorm(inner.constant) + (1-a)/(1-pnorm(inner.constant)))
    expectations[i,] = constant*cubintegrate(meanF, rep(-Inf,k), rep(Inf,k), fDim=k)$integral
  }
  
  existing_names = names(df)
  new_col_name = paste0("expectations", power)
  while (new_col_name %in% existing_names) {
    if (grepl("[0-9]$", new_col_name)) {
      number <- as.numeric(gsub(".*([0-9]+)$", "\\1", new_col_name)) + 1
      new_col_name <- gsub("[0-9]$", "", new_col_name)
      new_col_name <- paste0(new_col_name, number)
    } else {
      new_col_name <- paste0(new_col_name, "1")
    }
  }
  
  df <- df %>%
    mutate(!!new_col_name := expectations)
  
  return(df)
}

fitMeanModel = function(model, df){
  reg = lm(model, df)
  
  betas = coef(reg)
  return(betas)
}

linearCATE = function(z, coefs, params){
  k = sqrt(length(params$sigma.est))
  d = t(params$lambda) %*% solve(params$psi) %*% (z - params$nu)
  M = t(params$lambda) %*% solve(params$psi) %*% params$lambda + solve(params$sigma.est)
  lower.index = 3+k
  higher.index = 2+2*k
  return(c(coefs[lower.index:higher.index])%*%solve(M)%*%d + coefs[2])
}

ATE.est = function(df, params, coefs, method = linearCATE){
  ATE.est = 0
  Z = subset(df, select = -c(A, Y, expectations1))
  for(i in 1:nrow(df)){
    z = Z[i,]
    cate = method(t(z), coefs, params)
    ATE.est = ATE.est + cate
  }
  ATE.est = ATE.est/nrow(df)
  return(ATE.est)
}
