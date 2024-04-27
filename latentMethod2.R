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
  ZA = subset(df, select = -c(Y))
  
  psi.inv.centered = t(apply(ZA, 1, function(za.vec){solve(params$psi.est, za.vec - params$nu.est)}))
  psi.inv.lambda = solve(params$psi.est, params$lambda.est)
  k = ncol(params$lambda.est)
  getM = function(psi.inv.za){
    solve(t(params$lambda.est)%*%psi.inv.lambda+diag(k), t(params$lambda.est)%*%psi.inv.za)
  }
  
  M = apply(psi.inv.centered, 1, getM)
  return(M)
}


# regress Y onto tau (lambda'Z + gamma A)9

fitMeanModel = function(model, reg.df){
  reg = lm(model, reg.df)
  
  betas = coef(reg)
  return(betas)
}

linearCATE = function(z, coefs, params, a1, a2){
  gamma = c(betas[["M:A"]], betas[["A"]])
  p = length(z)
  k = ncol(params$lambda.est)
  psi.inv.centered = solve(params$psi.est[1:p,1:p], z - params$nu.est[1:p])
  psi.inv.lambda = solve(params$psi.est[1:p,1:p], params$lambda.est[1:p,])
  M = solve(t(params$lambda.est[1:p,])%*%psi.inv.lambda+diag(k), t(params$lambda.est[1:p,])%*%psi.inv.centered)
  
  cate = t(gamma) %*% c(M,1)*(a1-a2)
  return(cate)
}

ATE.est = function(Z, params, coefs, method = linearCATE){
  ATE.est = 0
  for(i in 1:nrow(Z)){
    z = Z[i,]
    cate = method(z, coefs, params, 1, 0)
    ATE.est = ATE.est + cate
  }
  ATE.est = ATE.est/nrow(Z)
  return(ATE.est)
}
