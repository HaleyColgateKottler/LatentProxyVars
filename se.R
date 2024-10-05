library(lavaan)
library(cubature)
library(mvtnorm)
library(dplyr)
library(reshape2)
library(ggplot2)

# tag = "1D"
# p=9
# k=1

calcCATE_SE <- function(data, constraints, output, a.dif, z) {
  n <- nrow(data)
  p <- nrow(constraints)
  k <- ncol(constraints)
  lambda.parms <- sum(constraints == 1)

  theta <- as.numeric(c(
    output$nu,
    output$lambda[constraints == 1],
    output$psi,
    output$alpha,
    output$gamma,
    output$beta
  ))

  # estimating function
  # if weights are supplied, then it uses REM estimating equations; otherwise uses EM estimating equations
  estFUN <- function(data, constraints, a.dif, z) {
    nu <- apply(data[, 1:p], 2, mean)
    n <- nrow(data)
    p <- ncol(data) - 1
    k <- ncol(constraints)
    lambda <- matrix(0, nrow = p, ncol = k)
    lambda.idx <- which(constraints == 1)
    lambda.parms <- length(lambda.idx)
    d.lambda <- matrix(0, nrow = n, ncol = lambda.parms)
    d.psi <- matrix(0, nrow = n, ncol = p)
    d.alpha <- matrix(0, nrow = n, ncol = k + 1)
    d.gamma <- matrix(0, nrow = n, ncol = k + 1)

    function(theta) {
      nu <- theta[1:p]
      lambda.est <- theta[(p + 1):(p + lambda.parms)]
      lambda[lambda.idx] <- lambda.est
      psi <- theta[(p + lambda.parms + 1):(p + lambda.parms + p)]
      alpha <- theta[(2 * p + lambda.parms + 1):(2 * p + lambda.parms + k + 1)]
      gamma <- theta[(2 * p + lambda.parms + k + 2):(2 * p + lambda.parms + 2 * k + 2)]
      beta <- theta[2 * p + lambda.parms + 2 * k + 3]
      # estimating equations
      Z <- apply(data[, 1:p], 1, function(y) y - nu)
      sigma <- lambda %*% t(lambda) + diag(psi)
      inv.sigma <- solve(sigma)
      m <- t(lambda) %*% solve(diag(psi)) %*% lambda + diag(k)
      f.base <- solve(m) %*% t(lambda) %*% solve(diag(psi))
      d.nu <- t(inv.sigma %*% Z)
      for (i in 1:n) {
        D <- inv.sigma %*% (diag(p) - (Z[, i] %*% t(Z[, i]) %*% inv.sigma))
        d.lambda[i, ] <- -(D %*% lambda)[lambda.idx]
        d.psi[i, ] <- -0.5 * diag(D)
        f <- c(1, f.base %*% Z[, i])
        d.alpha[i, ] <- 2 * (data[i, p + 1] - (alpha + gamma * data[i, p - 1]) %*% f) %*% t(f)
        d.gamma[i, ] <- data[i, p] * 2 * (data[i, p + 1] - (alpha + gamma * data[i, p]) %*% f) %*% t(f)
      }
      small.m <- t(lambda[1:(p - 1), 1:k]) %*% solve(diag(psi[1:(p - 1)])) %*% lambda[1:(p - 1), 1:k] + diag(k)
      U.mean <- solve(small.m) %*% t(lambda[1:(p - 1), 1:k]) %*% solve(diag(psi[1:(p - 1)])) %*% (z - nu[1:(p - 1)])
      d.beta <- a.dif * gamma %*% c(1, U.mean) - beta

      # score function for each person and parameter
      c(d.nu, d.lambda, d.psi, d.alpha, d.gamma, d.beta)
    }
  }

  # use geex package to estimate standard errors
  m.results.EM <- geex::m_estimate(
    estFUN = estFUN,
    data = data,
    roots = theta,
    compute_roots = FALSE,
    outer_args = list(
      constraints = constraints,
      a.dif = a.dif,
      z = z
    )
  )

  se <- sqrt(diag(geex::vcov(m.results.EM)))
  names(se) <- c(
    paste0("nu", 1:p),
    paste0("lambda", which(constraints == 1)),
    paste0("psi", 1:p),
    paste0("alpha", 1:(k + 1)),
    paste0("gamma", 1:(k + 1)),
    "beta"
  )

  return(se)
}

calcSE <- function(data, p, k, output, a.dif) {
  n <- nrow(data)
  lambda.parms <- (p + 1) * k

  theta <- as.numeric(c(
    output$nu,
    output$lambda,
    diag(output$psi),
    output$alpha,
    output$gamma,
    output$beta
  ))

  # estimating function
  estFUN <- function(data, p, k, a.dif) {
    nu <- apply(data[, 1:(p + 1)], 2, mean)
    n <- nrow(data)
    lambda <- matrix(0, nrow = p + 1, ncol = k)
    lambda.idx <- 1:((p + 1) * k)
    lambda.parms <- (p + 1) * k
    d.lambda <- matrix(0, nrow = n, ncol = lambda.parms)
    d.psi <- matrix(0, nrow = n, ncol = p + 1)
    d.alpha <- matrix(0, nrow = n, ncol = k + 1)
    d.gamma <- matrix(0, nrow = n, ncol = k + 1)

    function(theta) {
      nu <- theta[1:(p + 1)]
      lambda.est <- theta[(p + 2):(p + 1 + lambda.parms)]
      lambda[lambda.idx] <- lambda.est
      psi <- theta[(p + lambda.parms + 2):(p + lambda.parms + p + 2)]
      alpha <- theta[(2 * p + lambda.parms + 3):(2 * p + lambda.parms + k + 3)]
      gamma <- theta[(2 * p + lambda.parms + k + 4):(2 * p + lambda.parms + 2 * k + 4)]
      beta <- theta[2 * p + lambda.parms + 2 * k + 5]
      # estimating equations
      Z <- apply(data[, 1:(p + 1)], 1, function(y) y - nu)
      sigma <- lambda %*% t(lambda) + diag(psi)
      inv.sigma <- solve(sigma)
      m <- t(lambda) %*% solve(diag(psi)) %*% lambda + diag(k)
      f.base <- solve(m) %*% t(lambda) %*% solve(diag(psi))
      d.nu <- t(inv.sigma %*% Z)
      for (i in 1:n) {
        D <- inv.sigma %*% (diag(p + 1) - (Z[, i] %*% t(Z[, i]) %*% inv.sigma))
        d.lambda[i, ] <- -(D %*% lambda)[lambda.idx]
        d.psi[i, ] <- -0.5 * diag(D)
        f <- c(1, f.base %*% Z[, i])
        d.alpha[i, ] <- 2 * (data[i, p + 2] - (alpha + gamma * data[i, p]) %*% f) %*% t(f)
        d.gamma[i, ] <- data[i, p + 1] * 2 * (data[i, p + 2] - (alpha + gamma * data[i, p + 1]) %*% f) %*% t(f)
      }
      d.beta <- a.dif * gamma[1] - beta

      # score function for each person and parameter
      c(d.nu, d.lambda, d.psi, d.alpha, d.gamma, d.beta)
    }
  }

  # use geex package to estimate standard errors
  m.results.EM <- geex::m_estimate(
    estFUN = estFUN,
    data = data,
    roots = theta,
    compute_roots = FALSE,
    outer_args = list(
      p = p,
      k = k,
      a.dif = a.dif
    )
  )

  se <- sqrt(diag(geex::vcov(m.results.EM)))
  names(se) <- c(
    paste0("nu", 1:(p + 1)),
    paste0("lambda", 1:((p + 1) * k)),
    paste0("psi", 1:(p + 1)),
    paste0("alpha", 1:(k + 1)),
    paste0("gamma", 1:(k + 1)),
    "beta"
  )

  return(se)
}

# azGen(tag, 1000)
# raw.data = data.import(tag, 1000)
# ate.est <- latent.ATE(raw.data, 1:3, p)
# load(file.path("Data", "Parameters", paste("param_", tag, ".RData", sep = "")))
# output = list("nu" = c(Z_intercept, A_intercept),
#               "lambda" = lambda,
#               "psi" = diag(psi),
#               "alpha" = alpha,
#               "gamma" = gamma,
#               "beta" = ate.est)
#
# se.est = calcSE(raw.data, p, k, output, 1)
#
# raw.data <- data.import(tag, 1000)
# se.est = calcSE(raw.data, p, k, output, 1)
#
# ests = matrix(0, nrow = 1, ncol = 100)
# se.ests = matrix(0, nrow = 1, ncol = 100)
#
# for (sample.size in 1:1){
#   print(sample.size)
#   for (trial.num in 1:100){
#     print(trial.num)
#     azGen(tag, sample.size * 1000)
#     raw.data = data.import(tag, sample.size * 1000)
#     ate.est <- latent.ATE(raw.data, 1:3, p)
#     se.ests[sample.size,trial.num] = ate.est$se[35]
#     ests[sample.size,trial.num] = ate.est$estimate
#   }
# }
#
# samp.df <- data.frame('estimates' = ests[1,],
#                       'se' = se.ests[1,],
#                       'lowerCI' = 0,
#                       'upperCI' = 0,
#                       'cover' = 0)
# covered <- 0
# for (i in 1:100){
#   if (!(is.na(ests[1,i]) | is.na(se.ests[1,i]))){
#     lower.ci <- ests[1,i] - 1*se.ests[1,i]
#     upper.ci <- ests[1,i] + 1*se.ests[1,i]
#     samp.df[i,'lowerCI'] <- lower.ci
#     samp.df[i,'upperCI'] <- upper.ci
#     if (lower.ci < .5 && upper.ci > .5){
#       covered <- covered + 1
#       samp.df[i,'cover'] <- 1
#     }
#   }
# }
# samp.df$cover <- factor(samp.df$cover)
#
# write.csv(samp.df, paste("coverage_test_", tag, ".csv", sep = ""))
#
# ggplot(samp.df) + geom_hline(yintercept = .5) +
#   geom_point(aes(x=1:100, y = estimates, color = cover)) +
#   geom_errorbar(aes(x = 1:100, ymin = lowerCI, ymax = upperCI,
#                     color = cover))
# ggsave(paste("coverage_plot_", tag, ".png", sep = ""))
