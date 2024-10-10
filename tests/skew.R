library(sn)

skew_param_gen <- function(k, p, alpha, gamma, slant.alpha, tag) {
  factor.loadings <- sigmaSim(k, p, 1)
  lambda <- factor.loadings[[1]]
  psi <- factor.loadings[[2]]
  psi[p + 1, p + 1] <- 1

  var_Y <- .1
  
  omega <- sqrt(1/(1-2*slant.alpha^2/(pi*(1+slant.alpha^2))))
  xi <- -omega*slant.alpha/(sqrt(1+slant.alpha^2))*sqrt(2/pi)

  save(k, p, var_Y, lambda, psi, alpha, gamma, slant.alpha,
       xi, omega,
    file = file.path(
      "Data", "Parameters",
      paste('param_', tag, sample.size,
        ".RData",
        sep = ""
      )
    )
  )
}

skew_az_gen <- function(tag, sample.size) {
  load(file.path(
    "Data", "Parameters",
    paste("param_", tag, sample.size, ".RData", sep = "")
  ))

  H <- rsn(sample.size, xi, omega, slant.alpha)
  Z <- t(lambda[1:p, ] %*% t(H)) + rmvnorm(sample.size, sigma = psi[1:p, 1:p])
  colnames(Z) <- paste0("V", 1:p)
  A <- t(lambda[p + 1, ] %*% t(H)) + rnorm(sample.size)
  Y <- rowSums((alpha + A %*% gamma) * cbind(1, H)) + rnorm(sample.size, sd = var_Y)

  write.table(H, file = file.path("Data", "Samples", paste("H_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Z, file = file.path("Data", "Samples", paste("Z_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(A, file = file.path("Data", "Samples", paste("A_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Y, file = file.path("Data", "Samples", paste("Y_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
}

skew_test <- function(kvals, p, sample.size, reps, tag, savemarker = 100) {
  est.df <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(est.df) <- c("latent", "linear", "IPW", "IV", "proximal")
  latent <- c()
  ipw <- c()
  linear <- c()
  iv <- c()
  proximal <- c()

  for (j in 1:trials) {
    i <- j %% savemarker

    skew_az_gen(tag, sample.size)
    raw.data <- data.import(tag, sample.size)

    linear.ate <- linearEst(raw.data)
    linear[i] <- linear.ate
    ipw.ate <- continuousIPWest(raw.data)
    ipw[i] <- ipw.ate
    iv.ate <- IVest(raw.data)
    iv[i] <- iv.ate
    proximal.ate <- proximal_causal(raw.data, paste("V", 1:floor(p/2), sep = ""),
                                    paste("V", (floor(p/2)+1):p, sep = ""))
    proximal[i] <- proximal.ate

    ate.est <- latent.ATE(raw.data, kvals, p)
    latent[i] <- ate.est$estimate

    if (j %% savemarker == 0 | j == trials) {
      print(sample.size)
      print(j)
      est.df <- rbind(est.df, cbind(latent, linear, ipw, iv, proximal))
      write.table(est.df,
        file.path(
          "Data", "Estimates",
          paste("ests_", tag, as.character(sample.size),
            ".csv",
            sep = ""
          )
        ),
        sep = ",",
        row.names = FALSE
      )
      latent <- c()
      ipw <- c()
      linear <- c()
      iv <- c()
      proximal <- c()
    }
  }
}
