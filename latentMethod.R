fitUZA <- function(model, df, k, p) {
  smallerDF <- subset(df, select = -c(Y))
  
  fit <- sem(model, data = smallerDF, rotation = "varimax")
  
  parameters <- list(
    lambda.est = inspect(fit, what = "est")$lambda,
    psi.est = inspect(fit, what = "est")$theta,
    sigma.est = inspect(fit, what = "est")$psi,
    nu.est = inspect(fit, what = "est")$nu
  )
  return(parameters)
}

fitExpectations <- function(params, df) {
  ZA <- subset(df, select = -c(Y))
  psi.inv.centered <- t(apply(ZA, 1, function(za.vec) {
    solve(params$psi.est, za.vec - params$nu.est)
  }))
  psi.inv.lambda <- solve(params$psi.est, params$lambda.est)
  k <- ncol(params$lambda.est)
  getM <- function(psi.inv.za) {
    solve(t(params$lambda.est) %*% psi.inv.lambda + diag(k), t(params$lambda.est) %*% psi.inv.za)
  }
  
  M <- apply(psi.inv.centered, 1, getM)
  return(M)
}

fitExpectationsNCE <- function(params, df) {
  ZA <- subset(df, select = -c(Y))
  psi.inv.centered <- t(apply(ZA, 1, function(za.vec) {
    solve(params$psi.est, za.vec - params$nu.est)
  }))
  psi.inv.lambda <- solve(params$psi.est, params$lambda.est)
  k <- ncol(params$lambda.est)
  getM <- function(psi.inv.za) {
    solve(t(params$lambda.est) %*% psi.inv.lambda + diag(k), t(params$lambda.est) %*% psi.inv.za)
  }
  
  M <- apply(psi.inv.centered, 1, getM)
  return(M)
}


# regress Y onto tau (lambda'Z + gamma A)9

fitMeanModel <- function(model, reg.df) {
  reg <- lm(model, reg.df)
  betas <- coef(reg)
  return(betas)
}

linearCATE <- function(z, betas, params, a1, a2) {
  p <- length(z)
  k <- ncol(params$lambda.est)
  
  z.vars <- paste0("V", 1:p)
  M.vars <- paste0("M", 1:k)
  # gamma = c(betas[which(c(rep(0,k+1), 0, rep(1,k))==1)], betas[which(c(rep(0,k+1),1,rep(0,k))==1)])
  gamma <- c(betas[paste0(M.vars, ":A")], betas["A"])
  
  psi.inv.centered <- solve(params$psi.est[1:p, 1:p], z - params$nu.est[1:p])
  psi.inv.lambda <- solve(params$psi.est[1:p, 1:p], params$lambda.est[1:p, ])
  M <- solve(t(params$lambda.est[1:p, ]) %*% psi.inv.lambda + diag(k), t(params$lambda.est[1:p, ]) %*% psi.inv.centered)
  
  cate <- (a1 - a2) * gamma %*% c(M, 1)
}

ATE.est <- function(Z, params, coefs, method = linearCATE) {
  ATE.est <- 0
  for (i in 1:nrow(Z)) {
    z <- Z[i, ]
    cate <- method(z, coefs, params, 1, 0)
    ATE.est <- ATE.est + cate
  }
  ATE.est <- ATE.est / nrow(Z)
  return(ATE.est)
}


latent.ATE <- function(rawData, k, p) {
  z.vars <- names(rawData)[grep("^V", names(rawData))]
  h.vars <- paste0("efa('efa1')*h", 1:k)
  
  model <- paste(
    paste(h.vars, collapse = " + "),
    " =~ ",
    paste(z.vars, collapse = " + "),
    "+ A
                ",
    paste0("h", 1:k, sep = "", collapse = "~ 0*1\n"),
    "~0*1"
  )
  
  params <- fitUZA(model, rawData, k, p)
  
  M <- fitExpectations(params, rawData)
  
  AZ <- subset(rawData, select = -c(Y))
  Z <- subset(AZ, select = -c(A))
  
  if (k > 1) {
    M.df <- t(M)
  } else {
    M.df <- M
  }
  prior.names <- colnames(rawData)
  M.vars <- paste0("M", 1:k)
  rawData <- cbind(rawData, M.df)
  colnames(rawData) <- c(prior.names, M.vars)
  yModel <- paste(
    "Y ~",
    paste(M.vars, collapse = " + "),
    "+", paste0(M.vars, "*A", collapse = " + "),
    "+ A"
  )
  betas <- fitMeanModel(yModel, rawData)
  ATEest <- ATE.est(Z, params, betas)
}

latent.ATE.NCE <- function(rawData, k, p) {
  z.vars <- names(rawData)[grep("^V", names(rawData))]
  h.vars <- paste0("efa('efa1')*h", 1:k)
  
  model <- paste(
    paste(h.vars, collapse = " + "),
    " =~ ",
    paste(z.vars, collapse = " + "),
    "+ A
                ",
    "A ~", paste(z.vars, collapse = " + "),
    "
                ",
    paste0("h", 1:k, sep = "", collapse = "~ 0*1\n"),
    "~0*1"
  )
  
  params <- fitUZA(model, rawData, k, p)
  
  M <- fitExpectationsNCE(params, rawData)
  
  AZ <- subset(rawData, select = -c(Y))
  Z <- subset(AZ, select = -c(A))
  
  if (k > 1) {
    M.df <- t(M)
  } else {
    M.df <- M
  }
  prior.names <- colnames(rawData)
  M.vars <- paste0("M", 1:k)
  rawData <- cbind(rawData, M.df)
  colnames(rawData) <- c(prior.names, M.vars)
  yModel <- paste(
    "Y ~",
    paste(M.vars, collapse = " + "),
    "+", paste0(M.vars, "*A", collapse = " + "),
    "+ A"
  )
  betas <- fitMeanModel(yModel, rawData)
  ATEest <- ATE.est(Z, params, betas)
}
