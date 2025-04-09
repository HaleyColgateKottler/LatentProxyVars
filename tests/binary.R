binary_param_gen <- function(k, p, alpha, gamma, tag) {
  factor.loadings <- sigmaSim(k, p, 3)
  lambda <- factor.loadings[[1]]
  psi <- factor.loadings[[2]]
  psi[p + 1, p + 1] <- 1

  var_Y <- .1

  save(k, p, var_Y, lambda, psi, alpha, gamma,
    file = file.path(
      "Data", "Parameters",
      paste("param_", tag,
        ".RData",
        sep = ""
      )
    )
  )
}

binary_az_gen <- function(tag, sample.size) {
  load(file.path(
    "Data", "Parameters",
    paste("param_", tag, ".RData", sep = "")
  ))

  H <- rmvnorm(sample.size, mean = rep(0, k))
  Z <- t(lambda[1:p, ] %*% t(H)) + rmvnorm(sample.size, sigma = psi[1:p, 1:p])
  colnames(Z) <- paste0("V", 1:p)
  A.temp <- t(lambda[p + 1, ] %*% t(H)) + rnorm(sample.size)
  A <- ifelse(A.temp > 0, 1, 0)
  Y <- rowSums((alpha + A %*% gamma) * cbind(1, H)) + rnorm(sample.size, sd = var_Y)

  write.table(H, file = file.path("Data", "Samples", paste("H_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Z, file = file.path("Data", "Samples", paste("Z_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(A, file = file.path("Data", "Samples", paste("A_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Y, file = file.path("Data", "Samples", paste("Y_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
}

binary_test <- function(kvals, p, sample.size, reps, tag, savemarker = 100) {
  file.name <- file.path(
    "Data", "Estimates",
    paste("ests_", tag, as.character(sample.size),
          ".csv",
          sep = ""
    )
  )
  if (file.exists(file.name)){
    est.df <- read.csv(file.name)
  } else {
    est.df <- data.frame(matrix(nrow = 0, ncol = 7))
    colnames(est.df) <- c("latent", "low", "high", "linear", "IPW", "IV", "proximal")
  }
  latent <- c()
  low <- c()
  high <- c()
  ipw <- c()
  linear <- c()
  iv <- c()
  proximal <- c()

  for (j in 1:trials) {
    i <- j %% savemarker + 1
    binary_az_gen(tag, sample.size)
    raw.data <- data.import(tag, sample.size)

    linear.ate <- linearEst(raw.data)
    linear[i] <- linear.ate
    ipw.ate <- binaryIPWest(raw.data)
    ipw[i] <- ipw.ate
    iv.ate <- IVest(raw.data)
    iv[i] <- iv.ate
    proximal.ate <- proximal_causal(raw.data, paste("V", 1:floor(p/2), sep = ""),
                                    paste("V", (floor(p/2)+1):p, sep = ""))
    proximal[i] <- proximal.ate

    ate.est <- bootstrap.latent(raw.data, kvals, p, 6)
    
    latent[i] <- ate.est[1]
    low[i] <- ate.est[2]
    high[i] <- ate.est[3]

    if (j %% savemarker == 0 | j == trials) {
      print(sample.size)
      print(j)
      est.df <- rbind(est.df, cbind(latent, low, high, linear, ipw, iv, proximal))
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
      low <- c()
      high <- c()
      ipw <- c()
      linear <- c()
      iv <- c()
      proximal <- c()
    }
  }
}
