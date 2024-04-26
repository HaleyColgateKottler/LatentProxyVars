fitUZA = function(model, df, k, p){
  smallerDF = subset(df, select = -c(Y))

  fit <- sem(model, data=smallerDF, rotation = "varimax")
  
  parameters = list(
    lambda.est = inspect(fit,what="est")$lambda,
    psi.est = inspect(fit,what="est")$theta,
    sigma.est = inspect(fit,what="est")$psi,
    nu.est = inspect(fit,what="est")$nu
  )
  return(parameters)
}

fitExpectations = function(params, df){
  psi.est.inv = solve(params$psi.est)
  lambda.lambda = params$lambda.est %*% t(params$lambda.est)
  M.inv = (diag(k)-t(params$lambda.est) %*% solve(psi.est.inv + lambda.lambda) %*% params$lambda.est)
  
  ZA = subset(df, select = -c(Y))
  d = apply(ZA, 1, function(vec.ZA){ t(params$lambda.est) %*% psi.est.inv %*% (vec.ZA - params$nu.est)})
  Md = sapply(d, function(d.vec){M.inv %*% d.vec})
  return(Md)
}

fitMeanModel = function(model, reg.df){
  reg = lm(model, reg.df)
  
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
