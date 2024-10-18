iv_param_gen <- function(k, p, alpha, gamma, tag) {
  factor.loadings <- sigmaSim(k + p, 0, 2)
  lambda <- factor.loadings[[1]]
  psi <- factor.loadings[[2]]
  var_Y <- .1
  save(k, p, var_Y, lambda, alpha, gamma,
    file = file.path(
      "Data", "Parameters",
      paste("param_", tag,
        ".RData",
        sep = ""
      )
    )
  )
}

iv_az_gen <- function(tag, sample.size) {
  load(file.path(
    "Data", "Parameters",
    paste("param_", tag, ".RData", sep = "")
  ))

  Z <- rmvnorm(sample.size, sigma = diag(p))
  colnames(Z) <- paste0("V", 1:p)
  H <- rmvnorm(sample.size, mean = rep(0, k))
  A <- t(lambda %*% t(cbind(Z, H))) + rnorm(sample.size)
  Y <- rowSums((alpha + A %*% gamma) * cbind(1, H)) + rnorm(sample.size, sd = var_Y)

  write.table(H, file = file.path("Data", "Samples", paste("H_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Z, file = file.path("Data", "Samples", paste("Z_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(A, file = file.path("Data", "Samples", paste("A_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Y, file = file.path("Data", "Samples", paste("Y_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
}

iv_test <- function(kvals, p, sample.size, reps, tag, savemarker = 100) {
  est.df <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(est.df) <- c("latent", "linear", "IPW", "IV", "proximal")
  latent <- c()
  ipw <- c()
  linear <- c()
  iv <- c()
  proximal <- c()

  for (j in 1:reps) {
    i <- j %% savemarker
    iv_az_gen(tag, sample.size)
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
    if (!is.null(ate.est)) {
      est <- ate.est$estimate
    } else {
      est <- NaN
    }
    latent[i] <- est

    if (j %% savemarker == 0 | j == reps) {
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
